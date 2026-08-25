Installation
============

StarSmasher is compiled against whatever MPI and CUDA your machine provides, so
it is built from source rather than installed from a package.

What you need
-------------

* Linux
* An MPI implementation providing ``mpif90`` -- OpenMPI and MPICH both work
* A Fortran compiler.  gfortran 10 or later, or ifort
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

Building
--------

The gravity library first, then the code::

   $ make -C SPHgrav_lib2
   $ make

The gravity library reads your card's compute capability from ``nvidia-smi``, so
there is nothing to edit.  On a machine without a GPU, pass the two digits
yourself with ``make COMPUTE_CAPABILITY=75``.  Do not override ``NVCCFLAGS``
directly: a value given on the command line replaces the whole variable and
quietly drops the optimisation level and include paths with it.

A successful build ends with::

   ***MADE VERSION THAT USES GPUS***  ->  ../parallel_bleeding_edge_gpu_sph

The executable is named after the directory containing ``src``, and is copied one
level up.  ``make cpu`` produces ``..._cpu_sph`` in the same place.

.. note::

   Two settings are compile-time constants in ``starsmasher.h`` rather than
   inputs: ``nmax``, the maximum number of particles, and ``kdm``, the maximum
   number of zones in a stellar-evolution profile.  Raising either means editing
   that file and rebuilding.

   The static arrays are already close to a limit.  At the default ``nmax`` the
   block comes to about 2040 MB against the 2048 MB that ``-mcmodel=small``
   allows on x86-64, which is why the Makefile uses ``-mcmodel=medium`` where it
   is supported.  Raising the limits without it gives a relocation overflow at
   link time and an error that does not say what went wrong.

Checking that it works
----------------------

::

   $ python3 tests/run_tests.py

runs the test suite, which relaxes a polytrope and checks it against the virial
theorem, verifies energy conservation, and confirms that the answer does not
depend on the number of MPI ranks or on whether the GPU or CPU build is used.
