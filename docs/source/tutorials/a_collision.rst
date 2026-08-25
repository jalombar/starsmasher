A collision
===========

A collision starts from two relaxed stars.  This tutorial collides two
polytropes small enough that the whole thing runs in a few minutes, so that the
mechanics are clear before you commit to anything expensive.

Everything below was run as written.

Relax the two stars
-------------------

Follow :doc:`relaxing_a_polytrope` twice, in separate directories.  The stars
used here were

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 40

   * -
     - :math:`M`
     - :math:`R`
     - Check
   * - star 1
     - 0.4
     - 0.4
     - :math:`W = -0.342857`, theory :math:`-0.342857`
   * - star 2
     - 0.2
     - 0.2
     - :math:`W = -0.171429`, theory :math:`-0.171429`

both with ``n=8000``, ``nnopt=23``, relaxed to ``tf=20`` with ``treloff=0``.
Confirm each is settled before going on: the kinetic energy in column 3 of
``energy0.sph`` ended at :math:`6\times10^{-7}` and :math:`1\times10^{-7}`
respectively, against internal energies of order :math:`0.1`.

Set up the encounter
--------------------

Make a new directory, and bring in the last snapshot of each relaxation under
the names the code expects::

    cp ../star1/out0004.sph sph.start1u
    cp ../star2/out0004.sph sph.start2u

``sph.init`` selects the orbit routine::

    &INITT
    INAME='hyp'
    &END

and ``sph.input`` describes the orbit::

    &input
     tf=12,
     dtout=1,
     nrelax=0,
     sep0=5,
     rp=0.3,
     e0=1.0,
     ngravprocs=1,
    &end

``nrelax=0`` makes this a dynamical calculation rather than a relaxation.
``sep0`` is the separation the stars start at, ``rp`` the periastron separation
of the initial Keplerian orbit, and ``e0`` its eccentricity, so ``e0=1`` is
parabolic.

The stars have radii 0.4 and 0.2, summing to 0.6, and ``rp=0.3`` is well inside
that: they will hit rather than pass by.

Specifying the orbit
~~~~~~~~~~~~~~~~~~~~

``initialize_hyperbolic.f`` works out what it needs from whichever pair you give
it.  Any of ``(e0, vinf2)``, ``(e0, rp)``, ``(semimajoraxis, rp)``,
``(semimajoraxis, e0)`` or ``(rp, vinf2)`` will do.  Leave the others unset.

``vinf2`` is the square of the velocity at infinity and fixes the character of
the encounter: zero is parabolic, positive hyperbolic, negative elliptical.  It
relates to the semimajor axis by :math:`a = -GM/v_\infty^2` with :math:`M` the
total mass, so an elliptical encounter of two stars totalling
:math:`0.6\,M_\odot` with ``vinf2=-0.07`` has :math:`a = 0.6/0.07 = 8.6`.

Run it
------

::

    mpirun -np 2 ./run_gpu_sph

The run above took a few minutes, completing 3405 iterations and writing 13
snapshots.  The energy trace begins and ends::

    t= 0.000000  W=-0.5304919  T=0.01600068  U=0.2573550  Etot=-0.2571362
    t=11.99862   W=-0.4908685  T=0.02670166  U=0.2070898  Etot=-0.2570770

Total energy is conserved to :math:`2\times10^{-4}` across the encounter, which
is the first thing to check.  The kinetic energy rises as the stars fall
together and the potential becomes less negative as material is thrown outward.

.. note::

   ``energy0.sph`` has seven columns here, not nine.  The orbital angular
   frequency and separation are written only for a corotating-frame relaxation,
   ``nrelax >= 2``, not for a dynamical collision.  See :doc:`../using/output`.

Deciding what to resolve
------------------------

The tutorial above resolves both stars, which is the right choice for learning
and often the wrong one for research.

Before resolving a companion, ask whether its interior can exchange information
with the other star's fluid at all.  Compare the companion's radius with the
*smoothing length* of the primary in the region the companion will travel
through.  If the companion is much smaller than one smoothing length, the
primary's fluid cannot resolve its structure however many particles you give it,
and a point mass captures everything the fluid can feel.

The check that decides whether the calculation is worth running at all is a
different one: compare the accretion radius :math:`2GM/v^2` of the companion
against the same smoothing length.  That is the scale on which the companion
gathers material, and it is the interaction you are presumably trying to model.
If it is comfortably larger than a smoothing length, the encounter is resolved
even with a point-mass companion.  If it is *smaller*, no amount of resolving
the companion will help -- the primary itself is too coarsely resolved, and that
is what has to be fixed.

For a low-mass companion falling into an extended giant, the first comparison
typically says "do not bother resolving it" and the second says "the encounter
is nonetheless well resolved".  Both are worth computing before choosing.

Where to go next
----------------

:doc:`relaxing_a_star` covers building progenitors from stellar-evolution
profiles rather than polytropes, which is what most real work needs.
