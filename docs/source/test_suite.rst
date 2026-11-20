Test suite
==========

::

   python3 tests/run_tests.py

The suite runs short calculations and checks them against things that must be
true, rather than against stored output.  A golden-file test tells you that
something changed.  These tell you whether the answer is right.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Test
     - Checks
   * - ``virial``
     - A relaxed :math:`n=1.5` polytrope satisfies
       :math:`W = -\frac{3}{5-n}\frac{GM^2}{R}`, and :math:`U = -W/2`.
   * - ``energy``
     - Total energy is conserved through the run.
   * - ``quiet``
     - A relaxed star stays put: kinetic energy remains negligible.
   * - ``ranks``
     - The answer does not depend on how many MPI ranks are used.
   * - ``cpu_gpu``
     - The CPU and GPU builds agree.
   * - ``block_cpu_gpu``
     - With block timesteps, where gravity is computed on the active particles
       only, the CPU and GPU builds agree.
   * - ``gpu_hilo``
     - GPU gravity with and without the float remainders of the positions
       (``ngravhilo``) agrees for a single star.
   * - ``gpu_tree``
     - The tree gravity on the GPU (``ngravtree = 2``) agrees with the tree on
       the CPUs (``ngravtree = 1``) for the same opening angle and leaf size:
       energies, accelerations and potentials of a relaxed star; a second GPU run
       has the threshold of the leaves with a very long list so low that most
       leaves are divided in chunks.
   * - ``gravity_auto``
     - ``ngravtree = -1`` chooses the gravity from the number of particles, the
       threads, the ranks and the SPH (with a GPU: the direct sum below
       ``ngravtreeauto`` and the tree on the GPU from it up when it is written,
       the limits measured for the SPH that goes with the gravity by default;
       the tree on the CPUs in the CPU build), writes the choice in
       ``log0.sph``, and gives exactly the energies and dumps of the value it
       chose, set by hand.
   * - ``gpu_sph``
     - ``ngpusph = 1`` (the smoothing lengths, densities and neighbour lists on the GPU) gives exactly the bits of
       ``ngpusph = 0``: the energies of every row and the bytes of every dump, on
       one rank and on several; a run with ``ngpusphcheck = 1`` finishes (the CPUs
       and the GPU agree particle by particle); the CPU build says in ``log0.sph``
       that it has no GPU code and gives the same bits.
   * - ``gpu_forces``
     - ``ngpusph = 3`` (the pair sums of the forces on the GPU too) agrees with
       ``ngpusph = 1`` to the rounding of the sums, with shared and with block
       steps, on one rank and on several: the energies of every row to one unit of the last digit that
       ``energy0.sph`` prints (7 digits), and with ``ngpusphcheck = 2`` the accelerations to 1e-12 and du/dt to 1e-10 of
       the largest; the CPU build says in ``log0.sph`` that it has no GPU code and
       gives the same bits.
   * - ``gpu_ranks``
     - ``ngpusphcheck = 2`` and ``3`` with more ranks than GPUs
       (``CUDA_VISIBLE_DEVICES=0`` with 2 and with 3 ranks): the ranks without a
       GPU of their own take part in the comparison of the forces, the run
       finishes (it used to hang) and agrees with ``ngpusph = 1`` to 1e-9.
   * - ``gpu_resident``
     - ``ngpusph = 7`` (the candidates of the smoothing-length search kept on the
       GPU between steps, with a small margin so that they are built again in the
       run) gives exactly the bits of ``ngpusph = 3``: the energies of every row
       and the bytes of every dump, on one rank and on several, also with
       ``ngpusphcheck = 3``; the CPU build says in ``log0.sph`` that it has no GPU
       code and gives the same bits.
   * - ``gpu_mixed``
     - ``ngpusphprec = 1`` (the terms of the pairs of the forces on the GPU in
       single precision, sums in double) stays within the tolerances:
       with ``ngpusphcheck = 2`` the accelerations differ from the CPU's by less
       than 1e-4 and du/dt by less than 1e-3 of the largest at every step, and the
       energies are within 1e-4 of those in double precision; the CPU build says in
       ``log0.sph`` that it has no GPU code and gives the same bits.
   * - ``sph_auto``
     - ``ngpusph = -1`` chooses from the configuration, writes the
       choice and the reason in ``log0.sph`` (the CPU build chooses 0, a GPU build
       with the tree 7, with the direct sum by threads per rank, particles and
       ranks), and gives exactly the bits of a run with the value it chose set by
       hand, on one rank and on several; each limit of the rule is tried with a
       star below it and a star above it, for the mid and the strong class of the
       GPU (``ngpuclass``) where the limit depends on it; with ``ngpuclass = -1``
       the class that ``log0.sph`` says must agree with the number of
       multiprocessors it says it read, and ``ngpuclass = 3`` must be refused
       (in a build without OpenMP the cases that need more than one thread per rank
       are skipped).
   * - ``vdotsm``
     - With ``cn2`` and ``cn4`` finite (so that ``tstep`` uses the smoothed
       accelerations of ``vdotsm``, which sums over the neighbour lists on the
       CPU) ``ngpusph = 1`` gives exactly the bits of ``ngpusph = 0`` and
       ``ngpusph = 7`` those of ``ngpusph = 3``: the energies of every row, the
       bytes of every dump and the number of steps, on one rank and on several;
       the same run with the default ``cn2`` and ``cn4`` comes out different.
   * - ``gpu_multidev``
     - The SPH and the forces with several devices for the rank, tried on one GPU
       with ``SPHGPU_NDEV`` (logical devices dealt out round the GPUs): with 1, 2
       and 4 logical devices the run reaches the end, with ``ngpusphcheck = 2`` the
       accelerations differ from the CPU's by less than 1e-12 and du/dt by less than
       1e-10 of the largest, and the energies are those of one device to 1e-9; with
       3 logical devices the dump and the energies do not depend on
       ``SPHGPU_FLAGMOD``.
   * - ``guard_ngravprocs``
     - Asking for more gravity processes than there are GPUs is clamped rather
       than fatal.
   * - ``startfiles``
     - Bodies are read from the file names given in ``sph.input``, not from
       names compiled into the code.
   * - ``periodic``
     - A fully periodic box conserves momentum to round-off and energy closely,
       keeps every particle inside the box, and shows no density seam at its
       faces.
   * - ``periodic_block``
     - Block timesteps in a periodic box, with a hot half on shorter steps,
       give the same densities as shared timesteps to within 1% and keep
       every particle inside the box.
   * - ``guard_periodic``
     - A periodic box too small for the kernel (2h beyond half a period) stops
       the run with a message rather than missing neighbours.
   * - ``hdf5_dump``
     - With ``nhdf5=1`` the ``out*.h5`` snapshots hold exactly the numbers of
       the ``out*.sph`` ones, bit for bit, in every column and every entry of
       the header; with ``nhdf5=2`` the same files appear without
       ``out*.sph``; ``restartrad.sph`` stays native.  Skipped unless the
       build has HDF5 and ``h5py`` is installed.

Options::

   python3 tests/run_tests.py --list            # what exists
   python3 tests/run_tests.py virial energy     # named tests only
   python3 tests/run_tests.py --np 2            # fewer ranks
   python3 tests/run_tests.py --exe /path/to/..._sph   # skip the build
   python3 tests/run_tests.py --with-cpu        # also build the CPU version
   python3 tests/run_tests.py --quick           # smaller star, same assertions
   python3 tests/run_tests.py --hdf5            # build with make HDF5=1, so that hdf5_dump runs

Each build happens in a copy of the source under the suite's own scratch
directory, so a test run leaves the repository untouched.

Adding a test
-------------

Write a function taking a context and returning a one-line description of what
it verified, raise ``Fail`` with the numbers when it does not hold, and add it to
the ``TESTS`` dictionary in ``tests/run_tests.py``.  Prefer an invariant that has
to be true over a value that happened to be produced.
