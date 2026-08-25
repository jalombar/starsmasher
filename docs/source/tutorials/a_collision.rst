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
different one: compare the companion's accretion radius against the same
smoothing length.  That is the scale on which it gathers material, and it is
the interaction you are presumably trying to model.  Use the Bondi-Hoyle-Lyttleton
form

.. math::

   R_\mathrm{acc} = \frac{2GM}{v_\mathrm{rel}^2 + c_s^2},

and keep the sound speed.  Dropping it gives :math:`2GM/v^2`, which is only
valid for supersonic motion; a giant envelope is often transonic, where
:math:`c_s` is not a correction but half the denominator.

If :math:`R_\mathrm{acc}` exceeds a smoothing length, the encounter is resolved
even with a point-mass companion.  If it is *smaller*, no amount of resolving the
companion will help -- the primary itself is too coarse, and that is what has to
be fixed.

A worked example
~~~~~~~~~~~~~~~~

For a 0.4 :math:`M_\odot` companion inside a relaxed 5.4 :math:`M_\odot` giant,
measured from the model rather than estimated, taking :math:`v_\mathrm{rel}` as
the local circular velocity:

.. list-table::
   :header-rows: 1
   :widths: 12 16 14 14 12 16 16

   * - :math:`r`
     - :math:`M_\mathrm{enc}`
     - :math:`v_\mathrm{orb}`
     - :math:`c_s`
     - :math:`h`
     - :math:`R_\mathrm{acc}`
     - :math:`R_\mathrm{acc}/h`
   * - 50
     - 1.129
     - 65.7
     - 64.4
     - 11.8
     - 18.1
     - 1.53
   * - 100
     - 1.776
     - 58.2
     - 51.4
     - 12.9
     - 25.3
     - 1.96
   * - 150
     - 2.728
     - 58.9
     - 42.5
     - 14.0
     - 28.9
     - 2.06
   * - 200
     - 3.741
     - 59.7
     - 34.6
     - 15.3
     - 32.0
     - 2.09
   * - 250
     - 4.532
     - 58.8
     - 28.9
     - 17.3
     - 35.6
     - 2.06
   * - 300
     - 5.009
     - 56.4
     - 22.3
     - 22.6
     - 41.4
     - 1.83
   * - 340
     - 5.220
     - 54.1
     - 16.6
     - 30.2
     - 47.6
     - 1.57

Radii and smoothing lengths in :math:`R_\odot`, masses in :math:`M_\odot`,
velocities in km/s.

Two things follow.  The companion's radius is 0.369 :math:`R_\odot`, between 29
and 82 times smaller than the local smoothing length everywhere it can go, with a
density contrast running from :math:`2.5\times10^6` at the base of the envelope
to :math:`4\times10^8` at the surface.  Resolving it as fluid achieves nothing.

But the accretion radius exceeds the smoothing length everywhere, by a factor of
1.5 to 2.  So the encounter *is* resolved -- marginally at depth.  A factor of
1.5 is passing, not comfortable, and a coarser primary would fail this test.  A
reader who runs the check and gets 1.5 should not conclude they have room to
spare.

.. warning::

   Compute :math:`v_\mathrm{rel}` for your own encounter rather than copying a
   number.  The table above assumes a roughly circular orbit through gas that is
   not itself rotating much; a plunging or eccentric encounter is different.
   Because :math:`v` and :math:`c_s` are comparable through most of this
   envelope, the result is not sensitive to getting either slightly wrong, but it
   is very sensitive to dropping either.

   The Bondi-Hoyle-Lyttleton estimate is derived for a point mass in uniform gas.
   Inside a stratified envelope the density scale height is a further relevant
   length, and where :math:`R_\mathrm{acc}` approaches it the formula is being
   used outside its derivation.  Treat this as the standard first check, not as a
   sufficient condition.

Where to go next
----------------

:doc:`relaxing_a_star` covers building progenitors from stellar-evolution
profiles rather than polytropes, which is what most real work needs.
