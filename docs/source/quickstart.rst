Quickstart
==========

The fastest path from a clone to a running calculation.  If anything here fails,
:doc:`installation` covers the same ground in more detail.

You will need a Linux machine with MPI.  An NVIDIA card is worth having but is
not required.

.. code-block:: console

   $ git clone https://github.com/jalombar/starsmasher.git
   $ cd starsmasher

Check what the build will use::

   $ cd parallel_bleeding_edge/src
   $ make config

That prints the ``mpif90`` it found, the compiler behind it, the CUDA
installation and where the executables will go.  If a line is wrong, override it
on the command line rather than editing the Makefile.  See :doc:`installation`.

Then build::

   $ make

That compiles the gravity library as well, so it is the only build command you
need.  The last line should be::

   ***MADE VERSION THAT USES GPUS***  ->  ../parallel_bleeding_edge_gpu_sph

The executable is named after the directory holding ``src``, and is copied one
level up.

Without an NVIDIA card, use ``make cpu`` instead, which produces
``..._cpu_sph``.  See :doc:`using/running` for more on the cpu version of the code.

Now relax a star.  Work in a copy, so the pristine source stays clean.  The copy
brings the executable with it::

   $ cd ../..
   $ cp -r parallel_bleeding_edge my_first_star
   $ cd my_first_star
   $ cp ../example_input/relaxation_preMS/sph.in* .

Edit ``sph.input`` to something small, say ``n=10000`` and ``tf=30``, then run::

   $ mpirun -np 4 ./parallel_bleeding_edge_gpu_sph

You should see iterations counting up, and ``out0000.sph``, ``out0001.sph`` ...
appearing.  When it finishes, check ``energy0.sph``: the total energy in column 5
should be nearly constant and the kinetic energy in column 3 should be tiny.

That is a relaxed star, which is the starting point for everything else.
:doc:`tutorials/relaxing_a_polytrope` walks through the same run explaining what
each number means and how to tell a good model from a bad one.
