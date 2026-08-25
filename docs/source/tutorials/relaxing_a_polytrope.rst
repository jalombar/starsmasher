Relaxing a polytrope
====================

The quickest way to confirm that a StarSmasher build works is to relax a
polytrope.  It needs no stellar-evolution profile, no equation-of-state table
and no start files: the code generates the star itself from a mass and a radius.
It is also the first half of any collision, since a collision starts from two
relaxed stars.

This tutorial has been run exactly as written.  The numbers quoted are from that
run.

Setting up
----------

Copy the example into a working directory of its own::

    cp -r parallel_bleeding_edge relax_polytrope
    cd relax_polytrope
    cp ../example_input/relaxation_preMS/sph.in* .

``sph.init`` selects the calculation.  For a polytrope it contains::

    &INITT
    INAME='1es'
    &END

``sph.input`` describes the star.  The settings that matter here are

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Setting
     - Value used
     - Why
   * - ``n``
     - 10000
     - Number of particles.  Small enough to finish in a minute or two.
   * - ``starmass``
     - 0.2
     - Mass in solar masses.
   * - ``starradius``
     - 0.2
     - Radius in solar radii.
   * - ``nrelax``
     - 1
     - Relax a single star, rather than run a dynamical calculation.
   * - ``treloff``
     - 0
     - Keep the relaxation on for the whole run.
   * - ``tf``
     - 30
     - Stop at t = 30 in code units.
   * - ``dtout``
     - 1
     - Write a snapshot every unit of time.

Build and run
-------------

::

    cd src
    make
    cd ..
    mpirun -np 4 ./relax_polytrope_gpu_sph

The run above took a few minutes on four ranks, completing 5316 iterations and
writing 31 ``out*.sph`` snapshots.

Checking that it worked
-----------------------

``energy0.sph`` has one row per output, with columns for time, gravitational
potential energy :math:`W`, kinetic energy :math:`T`, internal energy :math:`U`,
and the total.  The first and last rows of the run were::

    t=0.000000   W=-0.1714286  T=0.000000     U=0.8571429E-01  Etot=-0.8571428E-01
    t=29.99793   W=-0.1714150  T=0.2975154E-07 U=0.8570280E-01  Etot=-0.8571213E-01

Three things to look for, all satisfied here.

**The star does not move.**  The kinetic energy stays at about
:math:`3\times10^{-8}`, which is zero to the precision that matters.  A
relaxation that is working ends with a star sitting still.

**Energy is conserved.**  The total drifts from :math:`-0.08571428` to
:math:`-0.08571213`, a relative change of :math:`2.5\times10^{-5}` over 5316
steps.

**The structure is right.**  For a polytrope of index :math:`n`, theory gives

.. math::

   W = -\frac{3}{5-n}\,\frac{GM^2}{R}.

With :math:`n = 1.5` and :math:`M = R = 0.2` in code units where :math:`G = 1`,
that is :math:`-0.171429`, which is what the run reports.  The virial theorem
then requires :math:`U = -W/2 = 0.085714`, which is also what the run reports.

If any of those three fails, something is wrong with the build or the settings,
and it is worth resolving before attempting anything larger.

Where to go next
----------------

The final snapshot is a relaxed star, ready to be used as one side of a
collision.  For stars read from a stellar-evolution profile, where the particle
layout needs more care, see :doc:`relaxing_a_star`.
