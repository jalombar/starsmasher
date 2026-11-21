SPH on the GPUs
===============

``ngpusph`` in ``sph.input`` puts parts of the SPH on the GPUs.  Its value is a bit mask; the default, 0,
is the code as it was.  Bit 1 (``ngpusph = 1``) solves the smoothing lengths on the GPU and makes the
densities, the other sums over the neighbours (``omega``, ``0mega``, ``psi`` and ``wn`` of the code) and
the neighbour lists there; everything else of the step (the pair forces, the pressure, the integrator)
stays where it was.  Bit 2 (``ngpusph = 3``, always together with bit 1) also makes the pair sums of the
forces on the GPU (below).

With ``ngpusph = 1`` the results are the same bits as with ``ngpusph = 0`` (for the same number of MPI ranks, and for any number of threads and GPUs).
With ``ngpusph = 3`` and ``7`` the pair sums of the forces are made on the GPU and the results are not the bits of the CPUs: they differ at the
level of the rounding of the forces (see the paragraph on the forces below), and ``ngpusph = 7`` gives the bits of ``ngpusph = 3``.

The kernel does the arithmetic of ``rho_and_h`` one operation at a time in double precision, without fused
multiply-add, in the same order (the sums run over the neighbours in the order of their index), with the
same tables of the kernel and the same steps of the root finder: the same sequence of trial smoothing lengths
and the same final sums.  A neighbour list is every particle at a squared distance below :math:`4h^2`,
in the order of the index, which is what the CPU makes.

How it works
------------

For each call the GPU makes a tree of boxes over all the particles (a Morton curve, leaves of 8 particles),
walks it with one thread for each particle of the rank, tests the points in single precision with a margin
(it can only let too many pass) and then in double precision, and keeps the candidates in the order of the
index.  A rank with several GPUs (``ngpuperrank``) divides its particles among them in fixed parts.  The
particles are done in chunks and the results of a chunk are copied to the CPU while the next one is computed.

A particle with more than 128 candidates is left to the CPU, which solves it with the code that always did
it; the CPU tree is then made, for those particles only.

When the CPU does the step
--------------------------

The criteria are fixed by the build and the input, never by how long something takes.  The CPU does the
whole step, and ``log0.sph`` says why (once, and again when the reason changes), in these cases: the build has
no GPU code (``make cpu``, the Metal build); the CPU code of the build fuses multiplications and additions
(``MARCH``, or aarch64), whose bits are not the GPU's; a periodic box; ``dynhco`` above 0; block time steps
(``nblock = 1``); a process without a GPU of its own; the first step of a run; the lists do not fit in ``list``; an error of CUDA
(the run goes on on the CPUs).  For the forces (bit 2) the reasons are the same, except that block steps and the first step
are done on the GPU, and that ``nav = 4`` (the Riemann solver) and ``nintvar = 3`` are left to the CPUs.

With the direct sum of the gravity on one GPU and ``ngravoverlap = 1`` the kernels of the gravity occupy the GPU while the smoothing
lengths are solved, which the CPUs would do meanwhile at no cost: ``ngpusph = 1`` is then slower for many particles (106 ms against
87 ms per step at 500000 particles) and faster for few (a 100000-particle star on one thread: 141 ms against 328 ms).  ``ngpusph``
is a choice, and it is made as asked.

The forces (``ngpusph = 3``)
----------------------------

Bit 2 makes ``div v`` and ``curl v`` (``getderivs``), the accelerations and ``du/dt`` (the pair loop of ``uvdots`` with
Balsara's viscosity) and ``uijmax`` on the GPU, for every particle that takes part in the step: all of them with shared
steps, the active ones and the ones refreshed for them with block steps (``nblock = 1``), the same particles as the CPU.
The pressure, the gravity, the integrator and the sums between the ranks stay where they were.

One thread for each particle forms the terms of its pairs with the operations, tables and order of the CPU
(double precision, no fused multiply-add) and keeps them; the terms are then sorted by the particle they act on and
added to it in the order the CPU adds them (the part from its own list first, then the part from the other particles
chunk by chunk, as ``pairbuf.f90`` does).  There is no atomic operation on floating-point numbers, so a run gives the same
bits every time, with any number of devices.  The results are the bits of the CPU, with one exception: the CPU takes
the root of the viscosity term with ``r2**0.5d0``, which the compiler turns into a call to ``pow``, and the GPU with the
exactly rounded ``sqrt``; the two differ in the last bit of some terms, so that where particles approach each other the
accelerations differ from the CPU's at the level of the rounding (about 1e-16 of the particle's), and the bits of a run
with ``ngpusph = 3`` are not those of ``ngpusph = 1``.  Momentum and angular momentum are conserved by the same
pairing of terms as on the CPU, to the rounding of the sums.

Checking and measuring
----------------------

* ``ngpusphcheck = 1`` makes the CPU do the work too and compares: the run stops at the first particle whose
  smoothing length, density, sums, number of neighbours or list has other bits, saying which quantity and the
  two values.  It is a diagnostic and costs the time of both.
* ``ngpusphcheck = 2`` makes the CPUs also form the pair sums of the forces after the GPU, on the same data, and writes one
  line to ``log0.sph`` for each evaluation with the differences: the largest and rms difference of the acceleration relative
  to the particle's own, the largest difference and the reference of ``du/dt``, ``uijmax`` and ``div v``, and how many
  particles have exactly the same bits.  The run goes on with the results of the CPUs.  ``ngpusphcheck = 3`` is both checks.  A rank that has no GPU
  of its own (more ranks than GPUs) or has lost it after an error of CUDA takes part in the comparison with the sums of its CPUs on both sides.
* ``SPHGPU_PROF=1`` in the environment writes to the error output of every rank the time of the phases of
  each call (upload and tree, kernel and copies, particles left to the CPU); ``SPHGPU_PROF=2`` also the
  events between the kernels.  ``SPHGPU_CHUNK`` is the number of particles of a chunk (131072).
* ``SPHGPU_FLAGMOD=K`` (a test of the code, not for a run): every Kth particle is left to the CPU as if it had too many candidates, to
  check that the particles solved by the CPU and by the GPU in the same step give the same bits as the CPU alone.
