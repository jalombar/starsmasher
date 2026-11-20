SPH on the GPUs
===============

``ngpusph`` in ``sph.input`` puts parts of the SPH on the GPUs.  Its value is a bit mask: 0 (the code as it was, which is the default), 1, 3 or 7;
-1 chooses one of 0, 1 and 7 from the configuration (below).  The SPH on the GPUs is off unless it is asked for; with a CUDA build
``ngpusph = -1``, together with ``ngravtree = -1``, is the value recommended for speed.  Bit 1 (``ngpusph = 1``) solves the smoothing lengths on the GPU and makes the
densities, the other sums over the neighbours (``omega``, ``0mega``, ``psi`` and ``wn`` of the code) and
the neighbour lists there; everything else of the step (the pair forces, the pressure, the integrator)
stays where it was.  Bit 2 (``ngpusph = 3``, always together with bit 1) also makes the pair sums of the
forces on the GPU (below).  Bit 4 (``ngpusph = 7``) keeps the candidate lists of the
smoothing-length search between steps (Verlet lists, below).

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
(it can only let too many pass), puts those that pass in the order of the index (a merge sort in the thread) and
then tests them in double precision, in that order, so that the candidates are in the order of the index.  A rank with several GPUs (``ngpuperrank``) divides its particles among them in fixed parts.  The
particles are done in chunks and the results of a chunk are copied to the CPU while the next one is computed.

The smoothing length of a particle is found by a root finder that evaluates a sum over the candidates (a function of the trial length) once for each trial.  Most
particles need two evaluations; the few that a shock has heated need tens or hundreds, one after the other, and a warp is as slow as its slowest thread.  At the collision
of two stars the few threads that make more than 8 evaluations (one in 200) took a third to a half of the time of the kernel on a GPU with slow double precision (the RTX cards, the L40S).  A thread therefore makes its first three
evaluations by itself, and when at most three threads of its warp are not done the 32 threads of the warp make the evaluations of each of them together: every thread has some of the candidates,
the terms of the sum are made at once and added in the order of the index, so the sums and the results are the same bits.  On an L40S the kernel takes 28 to 65% less time on dumps of the code (29% at the 500000-particle star, 55% at the pericentre of the collision), on an H100 9 to 49% less (33% at the pericentre).

A particle with more than 128 candidates is left to the CPU, which solves it with the code that always did
it; the CPU tree is then made, for those particles only.  With one rank the CPUs make that tree while the GPU solves the other
particles, if the step before left some to them (the tree is the one they would make afterwards from the same positions).  With several devices for the rank (``ngpuperrank``) the
neighbour lists stay on the devices in such a step too: the forces take the smoothing length, the sums and the list that the CPU made of
each of those particles to the device that owns it, so that the pair sums of a step with such particles are divided among the devices and added
in the same way as in a step without them.  A run with several devices therefore gives the same bits whichever steps have such particles
(``SPHGPU_FLAGMOD`` below changes nothing); with one device nothing changed.  The last digits of a run with several devices and such steps are not those of the
versions before this change: the sums are added in another grouping, so that they differ at the level of the rounding and the difference
grows as the run goes on.  In 200 steps of the collision of two Suns at 2 or 4 devices the dump of the new version differs from the old one by 1e-10 to 6e-9
(largest, positions) and 1e-6 to 1e-5 (largest, velocities), 3e-10 to 1e-9 in the median of the positions at the epoch with the most such steps,
and one or two of the 200 rows of ``energy0.sph`` change in the seventh digit.  Runs without such steps (no particle with more than 128 candidates:
the models of the tests, the stars at rest) and runs with one device have the bits of the old version.

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

The automatic choice (``ngpusph = -1``)
---------------------------------------

``ngpusph = -1`` chooses 0, 1 or 7 once, at the start of the run (and of every restart), from the configuration alone and never from a time measured, and writes
the choice and the reason to ``log0.sph`` (as ``ngravtree = -1`` does), as the fastest.  The rules come from measurements of the step time with ``ngpusph`` 0, 1 and 7 for
the direct sum on a GPU, the tree on a GPU and the tree on the CPUs (1, 4, 16 and 63 threads, 100000 and 500000 particles, 4 and 8 ranks and 8 GPUs;
for three or more GPUs, relaxed Suns of 9966 to 499860 particles and the collision of two Suns, 7 to 63 threads per rank, on 4 H100 and on 8 RTX PRO 6000, in the layouts of :doc:`gravity`;
relaxed Suns of 9966 to 499860 particles, 7, 16, 32 and 63 threads per rank, one rank and one GPU, on a Threadripper 7980X with a RTX 4070 Ti SUPER and on
an EPYC 9335 with a RTX PRO 6000; and, as a check, on a third node, a Xeon 8480C with a H100: 19938 to 499860 particles, 16, 32, 48 and 63 threads, mean step).  A choice is taken as the fastest if it is the fastest of the cells measured or lies within the noise of the measurement (the spread of the repetitions).

* No GPU code (CPU build, no GPU on the node), a build whose CPU code fuses multiplications and additions (``MARCH``), a periodic box, ``dynhco`` above 0: 0.
* ``nav = 4`` or ``nintvar = 3`` (their forces stay on the CPUs): 1.
* The tree walked on the GPU, the tree on the CPUs, or no gravity: 7 (tree on the GPU, 63 threads, 500000 particles: 24.5 ms against 37 with 0; 1 thread: 124
  against 1457; tree on the CPUs, 63 threads: 60 against 75).
* The direct sum on two GPUs (two ranks, or one rank with ``ngpuperrank = 2``): 7 below 450000 particles (the collision of two Suns, 199908 particles, 2 ranks of 8 threads:
  17.9 ms a step with 7, 26.1 with 1, 53.6 with 0) and 1 from there (500000 particles, one rank of 16 threads: 56.8 ms with 1, 62.9 with 7, 110.7 with 0).
* The direct sum on three or more GPUs (ranks times ``ngpuperrank``): 7 whatever the number of particles (4 H100, one rank of 63 threads: 3.77 ms a step with 7 against 4.92 with 1 at 9966 particles, 4.02 against 5.24 at 30090;
  at 499860 particles 1 is ahead by 2% with 63 threads and by 5% with 32; 8 RTX PRO 6000, 16 threads, 499860 particles: 30.34 ms with 7 against 52.01 with 1).  With ``ngpuperrank = 0`` (all the devices of the node divided among the processes),
  not measured: 1.
* The direct sum on one GPU, which the kernels of the gravity keep busy, depends on the threads per rank and on the particles, and with 32 threads or more on the class of the GPU: with at most 8 threads 7 (7 threads, 100000 particles: 25 ms against 60 with 0;
  1 thread, 500000 particles: 180 ms against 1425); with 9 to 31 threads 7 below 450000 particles (16 threads, 100000 particles: 24 ms against 30 on the Threadripper, 11 against 24 on the EPYC; 399816 particles: 70.7 ms against 91.3 on the EPYC)
  and 0 from there (499860 particles: 100 ms against 105); with 32 threads or more the limits are in the table below.  Below 9966 particles, the smallest star measured, the choice stays what it was in the earlier versions
  (at most 8 threads 7, 32 or more 0, in between 7).

**The class of the GPU.**  With 32 threads or more the CPUs take the SPH from the GPU earlier when the GPU is not a strong one, so the limits of the direct sum on one GPU depend on its class.  The class is read from the GPU: the number of multiprocessors
it says it has (``multiProcessorCount`` of CUDA; the smallest of the GPUs the process sees; rank 0 reads it and its choice is the one of every rank).  Below 100 the class is *mid*, from 100 up, or if the number cannot be
read (the CPU build, the Metal build, an error of CUDA), *strong*.  The number of multiprocessors is a coarse hint of the speed of the GPU: a RTX 4070 Ti SUPER has 66, a H100 132, a RTX PRO 6000 Blackwell 188.
``ngpuclass`` in ``sph.input`` writes the class by hand: 1 mid, 2 strong, -1 (the default) from the GPU; it serves to get the same choice on machines of different kinds, or to try the other class.  Any other value stops the run.
When the class decides the choice, the line of ``log0.sph`` says it and the number read, for example ``SPH on the GPU: ngpusph=-1 (auto): direct sum, 63 threads per rank and N=83000 >= 40000 on one GPU of class mid (66 SMs) -> ngpusph=0``,
or ``(class strong, written in ngpuclass)``.  ``ngpusph`` 7 is chosen below the number of particles of the table and 0 from there:

=================  ===================  ==============
threads per rank   strong class         mid class
=================  ===================  ==============
32 to 47           150000               75000
48 or more         75000                40000
=================  ===================  ==============

The cells that decide (step times, ms; the H100 ones are the mean step, the others the median of the repetitions).  Strong class, 32 to 47 threads: 99954 particles 10.6 with 7 against 14.3 with 0 on the EPYC, 16.1 against 20.1 on the H100;
199908 particles 37.4 with 0 against 42.0 with 7 on the H100; 499860 particles 157.6 against 171.0 on the H100 (0); any limit from 99955 to 199908 loses nothing in these cells.  Strong, 48 threads or more: 50160 particles 5.4 with 7
against 8.0 with 0 on the H100 with 63 threads, 5.35 against 9.2 with 48; 99954 particles 9.6 with 0 against 10.8 with 7 on the EPYC, 15.8 against 15.9 with 63 threads and 16.6 against 15.6 with 48 on the H100: from 50161 to 99954 the worst
cell is the H100 at 48 threads (0 chosen, 6.2% slower than 7), from 99955 on the EPYC loses 12.5%.  Mid class (the Threadripper with the RTX 4070 Ti SUPER, the only machine of this class that was measured), 32 threads: 50160 particles 7.9 with 7
against 9.3 with 0, 99954 particles 17.6 with 0 against 22.8 with 7, 199908 particles 45.8 against 65.6; 63 threads: 30090 particles 4.3 with 7 against 4.6 with 0, 50160 particles 6.8 with 0 against 8.0 with 7, 99954 particles 16.6 against 22.2;
48 threads (the mean of three repetitions): 19938 particles 3.3 with 7 against 3.7 with 0, 50160 particles 7.5 with 0 against 8.2 with 7, 99954 particles 15.9 against 23.0.  The limits of the two classes
are the numbers that lose least in the worst cell of the machines of the class, in the middle of the interval between the last cell where 7 wins and the first where 0 wins.

With ``ngravtree = -1`` as well (the default of ``ngravtreeauto`` is -1) the gravity and the SPH are chosen together, as the fastest pair: the tree walked on the GPU with the SPH on the GPU (7) from 9966 particles
on (see :doc:`gravity`).

The rule does not know how fast your GPU is, and the nodes do not always agree.  To find the fastest pair on your machine, for your input, measure it: :doc:`autotune` runs the
pairs for a few tens of steps and says which values to write.

Where the choice is 3 or 7 the results differ from the CPU's at the rounding level (the forces are formed once on the GPU, the root of the viscosity term differs), as
in the paragraph below; ``ngpusph = 0`` gives back the CPU code of earlier versions.  The measurements are for three kinds of node (a Threadripper with a RTX 4070 Ti SUPER,
an EPYC with RTX PRO 6000 Blackwell GPUs, a Xeon 8480C with a H100, on which the limits of the direct sum on one GPU were checked and moved): one is of the mid class and two of the strong one,
and the number of multiprocessors says only roughly how fast a GPU is against the CPUs of its node, so other hardware may need other rules (``ngpuclass`` to try the other class).

The forces (``ngpusph = 3``)
----------------------------

Bit 2 makes ``div v`` and ``curl v`` (``getderivs``), the accelerations and ``du/dt`` (the pair loop of ``uvdots`` with
Balsara's viscosity) and ``uijmax`` on the GPU, for every particle that takes part in the step: all of them with shared
steps, the active ones and the ones refreshed for them with block steps (``nblock = 1``), the same particles as the CPU.
The pressure, the gravity, the integrator and the sums between the ranks stay where they were.  With one rank the CPUs make, while the
GPUs do the pair sums, the sums of the energy line of the log that do not need the potential (the box, the extremes of the neighbour numbers,
density, u, h and pressure, the internal and kinetic energies); the sum of the potential is made when the gravity is collected.

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

Verlet lists (``ngpusph = 7``)
------------------------------

The search of the neighbours is the dearest part of the smoothing-length kernel (the visit of the tree and the sorting of the candidates by index).
With bit 4 the candidates of the first search of each particle, which are every particle within a radius larger by ``ngpusphskin`` (0.3 by default)
than the one asked for, are kept on the GPU from the call that made them.  Later calls read them from there, in the order of the index and with the
squared distances of the CPU, as long as the lists are complete: the radius they cover from the positions of the build, less twice the largest
displacement of a particle since then, must reach the radius the search needs.  A particle whose list does not reach goes to the tree as before, and
when more than 2% of the lists do not reach the next call builds them again.  The criterion is geometric, never a time.  The candidates are the same,
in the same order, and the final filter on the squared distance is the same, so the results are the bits of ``ngpusph = 3``, for any value of
``ngpusphskin``.  It works with one GPU for the rank; with several the search is as without it.  With the forces on the GPU (``ngpusph = 7``, shared steps)
the neighbour lists then stay on the GPU, where the forces read them, and are not copied to the CPU; they are copied only if the CPU needs them (particles it
solves itself, an error of CUDA in the forces, or a time step that uses the smoothed accelerations: ``nrelax`` other than 1 with ``cn2`` or ``cn4`` below
``1d29``, the default is ``1d30``, which makes ``tstep`` sum the accelerations over the lists on the CPU).  A particle whose list would hold more than 128
candidates has none.  ``SPHGPU_PROF=1`` writes whether the lists were built or used, the largest displacement and the number of lists that did not reach.

The host side of a step (one device)
------------------------------------

The step of the collision of two Suns (199908 particles) on an H100 with 63 threads takes 9 to 12 ms, and the GPU was idle for 2.1 to 2.7 ms of it: the CPU makes the arrays the next kernel needs, and the results
of one kernel make their way back before the next can start.  The code now shortens that path, without changing a bit of any result (``tests/run_tests.py gpu_hostpath``):

* With one device and the lists kept on it (``ngpusph = 7``, up to 8 chunks, that is up to about a million particles) the copies of the results of the smoothing-length kernel are queued behind the kernels before the CPUs start their own work
  (the tree of the gravity), and a kernel puts the records in the order of the particles first, so that what comes back is what the CPU reads (hp, rho, the number of neighbours and whether the particle was left to the CPU:
  24 bytes for each particle instead of 64) and the CPU stores it in order.  The rest of the records (omega, 0mega, psi, wn and where the lists are) is read only by the pair sums of the CPU, which are not made in such a run; it is
  brought and put in the arrays of the CPU if they are needed (``gpu_sph_getlists``, or forces made from the arrays of the CPU).  ``SPHGPU_SLOWDRAIN=1`` gives the earlier way (a test of the code, with the same bits).
* The velocities and u, which the predictor made, are sent to the device and packed with the positions while it solves the smoothing lengths, and the forces send only the pressure term (in the steps of the integration, with shared steps and ``nrelax = 0``;
  ``SPHGPU_NOEARLY=1`` gives the earlier way;
  ``SPHGPU_FECHECK=1`` makes the forces compare the arrays of the CPU with what was sent and stop at a difference, a test of the code).
* The positions are copied first and the masses, u, h and div v on a second stream while the box and the sort are made.
* The particles left to the CPU (more than 128 candidates) are listed once, by the threads, and only those are solved on the CPU; the count and the list take one parallel region, and the loop that used to scan all the particles with chunks of 64 is gone.
* The sums that come back from the forces, the results of the gravity of the tree, and the accelerations are copied and added in fewer passes and regions, each addition as it was.

Measured on the H100 of the DGX with 63 threads and one device (mean step, 5 repetitions alternated with the code before the change): the collision of two Suns at its start 9.49 -> 8.34 ms (-12.1%), at t=8 12.01 -> 10.42 (-13.2%),
at the pericentre 11.05 -> 9.57 (-13.3%) and late 11.48 -> 10.18 (-11.3%), the star of 500000 particles 23.04 -> 20.19 (-12.4%) and that of 2 million particles 95.74 -> 91.70 (-4.2%), which has 16 chunks and keeps the earlier copies of the results.
The GPU was idle 2.1 to 2.7 ms of a step of the collision, and is now idle 1.2 to 1.7 ms.  Variables of the environment
(not in the code): ``OMP_WAIT_POLICY=active`` or a large ``GOMP_SPINCOUNT`` keep the threads awake during the waits for the GPU and gave 1.4 to 1.7% on the same cells;
``OMP_WAIT_POLICY=passive`` costs 36 to 79%.

Single precision for the forces (``ngpusphprec = 1``)
-----------------------------------------------------

Never the default.  The terms of the pairs of the forces are made in single precision: the differences of positions and velocities are made in double
precision and rounded to float, the term is float (the kernel tables too), and what is added up (accelerations, du/dt, ``uijmax``) is double.  ``div v`` and
``curl v``, the smoothing lengths and the densities stay in double precision.  The results are not the bits of ``ngpusphprec = 0`` but differ from them, and from
the CPU's, by about 1e-6 of the acceleration of a particle (``ngpusphcheck = 2`` writes the differences); the run is as repeatable as the double precision one
(the same bits every time).  The tolerances it is qualified with are those of the test ``gpu_mixed`` (see :doc:`test_suite`).

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
* ``SPHGPU_NDEV=K`` (a test of the code, not for a run): the SPH and its forces use K logical devices, dealt out round the ``ngpuperrank`` physical devices of
  the rank (logical device k is on physical device ``first + k % ngpuperrank``), each with its own stream and buffers; the gravity keeps the physical ones.  On one
  GPU it gives the code paths, the division of the work and the bits of K devices, which is what the tests and the checks of the bits need; the times are not those
  of K GPUs.
* ``SPHGPU_FLAGMOD=K`` (a test of the code, not for a run): every Kth particle is left to the CPU as if it had too many candidates, to
  check that the particles solved by the CPU and by the GPU in the same step give the same bits as the CPU alone.
