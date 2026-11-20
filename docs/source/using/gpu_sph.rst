SPH on the GPUs
===============

``ngpusph`` in ``sph.input`` puts parts of the SPH on the GPUs.  Its value is a bit mask; the default, 0,
is the code as it was.  Bit 1 (``ngpusph = 1``) solves the smoothing lengths on the GPU and makes the
densities, the other sums over the neighbours (``omega``, ``0mega``, ``psi`` and ``wn`` of the code) and
the neighbour lists there; everything else of the step (the pair forces, the pressure, the integrator)
stays where it was.

The results are the same bits as with ``ngpusph = 0``, for any number of MPI ranks, threads and GPUs.
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
(the run goes on on the CPUs).

With the direct sum of the gravity on one GPU and ``ngravoverlap = 1`` the kernels of the gravity occupy the GPU while the smoothing
lengths are solved, which the CPUs would do meanwhile at no cost: ``ngpusph = 1`` is then slower for many particles (106 ms against
87 ms per step at 500000 particles) and faster for few (a 100000-particle star on one thread: 141 ms against 328 ms).  ``ngpusph``
is a choice, and it is made as asked.

Checking and measuring
----------------------

* ``ngpusphcheck = 1`` makes the CPU do the work too and compares: the run stops at the first particle whose
  smoothing length, density, sums, number of neighbours or list has other bits, saying which quantity and the
  two values.  It is a diagnostic and costs the time of both.
* ``SPHGPU_PROF=1`` in the environment writes to the error output of every rank the time of the phases of
  each call (upload and tree, kernel and copies, particles left to the CPU); ``SPHGPU_PROF=2`` also the
  events between the kernels.  ``SPHGPU_CHUNK`` is the number of particles of a chunk (131072).
* ``SPHGPU_FLAGMOD=K`` (a test of the code, not for a run): every Kth particle is left to the CPU as if it had too many candidates, to
  check that the particles solved by the CPU and by the GPU in the same step give the same bits as the CPU alone.
