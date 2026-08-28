Installation
============

StarSmasher is compiled against whatever MPI and CUDA your machine provides, so
it is built from source rather than installed from a package.

What you need
-------------

* Linux
* An MPI implementation providing ``mpif90``.  OpenMPI and MPICH both work
* A Fortran compiler.  Tested on gfortran 10 and later, and on ifort
* Optionally, an NVIDIA card and CUDA.  Without one, ``make cpu`` builds a
  working executable

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

   MPI
     mpif90          /usr/lib64/openmpi/bin/mpif90
     wrapping        GNU
     fixed-form flag -ffixed-line-length-132 -fallow-argument-mismatch
   CUDA
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

Building
--------

One command, from ``parallel_bleeding_edge/src``::

   $ make

That builds the gravity library too, so there is no need to build it
separately.  To rebuild only the library, which is occasionally useful when
chasing a CUDA problem, use ``make -C SPHgrav_lib2``.

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

   Two settings are compile-time constants in ``starsmasher.h`` rather than
   inputs: ``nmax``, the maximum number of particles, and ``kdm``, the maximum
   number of zones in a stellar-evolution profile.  They ship as 900000 and
   5000.  Raising either means editing that file and rebuilding.

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
