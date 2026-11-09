Installation
============

StarSmasher is compiled against whatever MPI and CUDA your machine provides, so
it is built from source rather than installed from a package.

What you need
-------------

* Linux, or macOS on Apple Silicon
* An MPI implementation providing ``mpif90``.  OpenMPI and MPICH both work
* A Fortran compiler.  Tested on gfortran 10 and later, and on ifort
* Optionally, to write the snapshots in HDF5 as well (see :doc:`using/hdf5`),
  the HDF5 library with its headers (``libhdf5-dev`` on Debian and Ubuntu) and
  an MPI C compiler, ``mpicc``
* On Linux, optionally an NVIDIA card and CUDA.  On Apple Silicon, the
  integrated GPU is supported through Metal.  Without a GPU backend,
  ``make cpu`` builds a working executable

Apple Silicon users should follow the `Apple Silicon installation guide
<https://github.com/jalombar/starsmasher/blob/master/documentation/installation_applesilicon.md>`_.
The build detects ``Darwin/arm64`` and selects Metal automatically; CUDA is not
required and no source or Makefile edits are needed.

The repository's own `documentation/installation.md
<https://github.com/jalombar/starsmasher/blob/master/documentation/installation.md>`_ covers setting
these up, including installing CUDA and OpenMPI and dealing with ``module``
environments on clusters.  What follows is the short version.

Getting the code
----------------

.. code-block:: console

   $ git clone https://github.com/jalombar/starsmasher.git

Clone rather than downloading an archive: the code is compiled against your own
libraries, and a clone can be updated later.

Checking what the build will use
--------------------------------

From ``parallel_bleeding_edge/src``::

   $ make config

which prints something like::

   Platform
     OS               Linux
     architecture     x86_64
     GPU backend      cuda
   MPI
     mpif90          /usr/lib64/openmpi/bin/mpif90
     wrapping        GNU
     fixed-form flag -ffixed-line-length-132 -fallow-argument-mismatch
   CUDA (used only when GPU backend is cuda)
     nvcc            /usr/local/cuda/bin/nvcc
     CUDAPATH        /usr/local/cuda
     runtime libdir  /usr/local/cuda/lib64
   Compilation
     FFLAGS          -O4 -mcmodel=medium
     code model      -mcmodel=medium
   Output
     GPU executable  ../parallel_bleeding_edge_gpu_sph
     CPU executable  ../parallel_bleeding_edge_cpu_sph

Everything there is detected, not hardcoded.  If a line is wrong, override it on
the command line rather than editing the Makefile:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - If this is wrong
     - Set
   * - the ``mpif90`` path, or the compiler behind it
     - ``make MPIF90=/full/path/to/mpif90``
   * - ``CUDAPATH``, or the runtime library directory
     - ``make CUDAPATH=/opt/cuda``
   * - the optimisation flags
     - ``make OLEVEL=-O2``
   * - the compute capability, if ``nvidia-smi`` cannot report it
     - ``make COMPUTE_CAPABILITY=120``
   * - where HDF5 is, for ``make HDF5=1``
     - ``make HDF5=1 HDF5ROOT=/path/to/hdf5``, or ``HDF5_CFLAGS`` and ``HDF5_LIBS``

Building
--------

One command, from ``parallel_bleeding_edge/src``::

   $ make

That builds the gravity library too, so there is no need to build it
separately.  To rebuild only the library, which is occasionally useful when
chasing a CUDA problem, use ``make -C SPHgrav_lib2``.

On Apple Silicon the same command builds ``SPHgrav_lib_metal`` instead.  The
CPU-only target remains available with ``make cpu``, but the Metal build is
recommended for normal runs.  Current M1 Pro measurements are about five times
faster than CPU-only runs for the tested simulations.  As a rough estimate for
gravity-dominated work, Pro-class chips should be about five times faster with
Metal than with CPU alone, and Max-class chips about ten times faster.

The library works out your card's compute capability for itself, so there is
normally nothing to set.

.. dropdown:: Setting the compute capability by hand
   :icon: gear

   Where the machine you compile on differs from the machine you run on, such
   as a login node with no card in it or a cluster with mixed hardware, give the
   value explicitly::

       $ make COMPUTE_CAPABILITY=120

   Write it as the capability times ten, with no decimal point: 9.0 is ``90``,
   12.0 is ``120``.  The setting carries from ``src`` down into the gravity
   library, so it works on either ``make`` line.

   If you are unsure of the number, ``__nvcc_device_query``, which ships with
   CUDA, prints it.  If you are unsure which of several cards to target, choose
   the lowest: code built for a lower capability still runs on a higher card,
   but not the other way round.

A successful build ends with::

   ***MADE VERSION THAT USES GPUS***  ->  ../parallel_bleeding_edge_gpu_sph

The executable is named after the directory containing ``src``, and is copied one
level up.  ``make cpu`` produces ``..._cpu_sph`` in the same place.

.. note::

   ``nmax``, ``nnmax`` and ``ntab`` are generated by the Makefile.  Their
   defaults can be overridden without editing source files, for example with
   ``make NMAX=200000 NNMAX=128 NTAB=400000``.  Apple Silicon uses smaller
   defaults for ``nmax`` and ``nnmax``, suitable for unified-memory Macs; other
   platforms keep the original values.  ``ntab`` is the number of entries of
   the kernel tables, which are read with linear interpolation: the default of
   200001 is more accurate than the former 10\ :sup:`7` entries read without
   interpolation, except for the cubic spline within 0.03\ *h* of a particle.
   ``kdm`` remains a compile-time constant in ``starsmasher.h``.

HDF5 snapshots (optional)
-------------------------

The snapshots are ``out*.sph`` files, a format of the code's own.  To have
``out*.h5`` files as well, which Python reads directly (see :doc:`using/hdf5`),
build with ``HDF5=1``::

   $ make HDF5=1            # the GPU version
   $ make cpu HDF5=1        # the CPU version

and set ``nhdf5=1`` in ``sph.input``.  Without ``HDF5=1`` the build is the same
as before and needs nothing new, not even a C compiler.  With it you need the
HDF5 library with its headers, and ``mpicc``, which compiles the one C file.

**Installing HDF5.**  Any HDF5 from 1.10 on will do, in its serial form or with
MPI.  The writer calls only the C interface of the library, not the Fortran
modules, whose ``.mod`` files depend on the version of the Fortran compiler, so
the HDF5 on your machine does not have to have been built with yours.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Where
     - How
   * - Debian, Ubuntu
     - ``sudo apt install libhdf5-dev pkg-config`` (or ``libhdf5-openmpi-dev``)
   * - macOS, Homebrew
     - ``brew install hdf5 pkg-config``
   * - A cluster
     - ``module load hdf5``, or ``spack load hdf5``, if your site has one; then
       give the directory with ``HDF5ROOT``, below
   * - Python, to read the files
     - ``pip install h5py numpy matplotlib`` (pandas too, for ``to_dataframe``)

**Where the build looks.**  In this order:

#. ``HDF5ROOT``: the directory that holds ``include`` and ``lib``, as in
   ``make HDF5=1 HDF5ROOT=/opt/hdf5``.
#. ``HDF5_CFLAGS`` and ``HDF5_LIBS``, if you give them yourself::

      $ make HDF5=1 HDF5_CFLAGS='-I/opt/hdf5/include' HDF5_LIBS='-L/opt/hdf5/lib -lhdf5'

#. Otherwise ``pkg-config``, asked for ``hdf5-serial``, then ``hdf5``, then
   ``hdf5-openmpi``, and the first that exists is used.  Without ``pkg-config``
   the link line gets a plain ``-lhdf5``.

**Checking it.**  ``make config HDF5=1`` prints ``mpicc``, ``HDF5_CFLAGS`` and
``HDF5_LIBS``, which is what the build is going to use; if one is wrong, set it
as above.  The link line of the build then ends with the HDF5 library.  A run
that asks for ``nhdf5`` above 0 with a build made without ``HDF5=1`` stops at
the start with ``init: nhdf5>0 needs a build with make HDF5=1``, so a build that
lacks the writer cannot be mistaken for one that has it.  The test suite has a
test for the writer, ``python3 tests/run_tests.py --hdf5 hdf5_dump``, which needs
``h5py``; without HDF5 in the build, or without ``h5py``, it is skipped with a
message.

Checking that it works
----------------------

From the top of the repository, not from ``src``::

   $ cd ../..
   $ python3 tests/run_tests.py

The suite builds its own copy of the code and runs seven short calculations,
each checking something that has to be true rather than comparing against stored
output.  It relaxes a polytrope and checks it against the virial theorem,
verifies energy conservation, confirms the star stays put once relaxed, checks
that the answer does not depend on the number of MPI ranks or on whether the GPU
or CPU build is used, and confirms that two input guards still work.

Expect about seven minutes on a GPU, and longer without one.  For a quicker
check, ``--quick`` uses a smaller star and finishes in well under a minute::

   $ python3 tests/run_tests.py --quick

Nothing about what is asserted changes, only the size of the model.  Use
``--list`` to see the individual tests, and name one to run it alone.
