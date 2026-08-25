Relaxing a star
===============

Before two stars can be collided, each must be relaxed: laid out as SPH
particles and then allowed to settle into hydrostatic equilibrium.  A collision
started from a poorly relaxed model will show spurious oscillations that have
nothing to do with the encounter, so it is worth spending time here.

This page covers relaxing a star read from a stellar-evolution profile
(``INAME='erg'``).  For a polytrope, see :doc:`relaxing_a_polytrope`.

What a relaxation run does
--------------------------

Setting ``nrelax=1`` turns on an artificial drag force that damps particle
motion, letting the model settle.  ``trelax`` sets the timescale of that drag
and ``treloff`` the time at which it is switched off.

Choosing ``trelax``
~~~~~~~~~~~~~~~~~~~

Set ``trelax`` to an oscillation period visible in the energy output at early
times.  Run the relaxation once, plot the energies from ``energy0.sph``, measure
the period of the oscillation you see, and use that.  The right value depends on
the star's structure and radius, so there is no universal number: a compact
main-sequence star and an extended giant will want values differing by orders of
magnitude.

Judging whether the model is any good
-------------------------------------

The relaxation writes two things worth comparing:

``parent.sph``
   the profile StarSmasher was asked to reproduce.

``col0000.sph``, ``col0001.sph``, ...
   the SPH model itself, dumped periodically.

Comparing ``parent.sph`` against ``col0000.sph`` shows how well the initial
particle layout matches the target.  Comparing it against the last ``col`` file
shows what the relaxation did with it.

Plot, against radius: pressure, density, particle mass, smoothing length,
neighbour number, and the radial components of the hydrodynamic and
gravitational accelerations.  In a good model the hydrodynamic acceleration
mirrors the gravitational one.  Where the two fail to balance, the model is not
in equilibrium and the relaxation will have to fight it.

.. note::

   A useful symptom to recognise: individual particles shooting away at high
   speed and returning, or settling at radii more than twice the stellar radius.
   The energy plot and the surface density can both look perfectly healthy while
   this is happening, because the particles involved carry very little mass.  It
   still means computational effort is being spent on a region that is not being
   modelled properly.

Judging a model after release
-----------------------------

Two mistakes are easy to make when deciding whether a relaxed model is settled,
and both produce confident wrong answers.

**Do not judge it while the drag is still on.**  Under drag everything looks
quiet, because the drag is what is making it quiet.  Configurations that appear
identical during relaxation can differ substantially once released.  Judge after
``treloff``, and allow several oscillation periods before concluding anything:
an effect that has not appeared yet is indistinguishable from one that does not
exist.

**Do not fit a trend through the first period after release.**  Switching the
drag off is itself a disturbance, and the first period contains the excursion it
produces.  Fitting through it measures the release, not the model.  This is not
a small correction: on a real model, a fit of :math:`\ln v_\mathrm{rms}` started
at the moment of release gave a *positive* slope, apparently a star coming
apart, while the same fit started one period later gave a negative slope, a star
settling down.  Same run, opposite conclusions.

Discard the first period, or plot the trace before fitting anything to it.

Comparing two models
~~~~~~~~~~~~~~~~~~~~

Residual velocity in km/s is not comparable between stars of different
compactness.  A compact star has a much higher sound speed, so the same
dimensionless disturbance shows up as a much larger absolute velocity.  Compare
:math:`v_\mathrm{rms}/c_s` instead.  Two models differing by a factor of
seventeen in km/s can differ by only a factor of two once scaled, which is the
number that means anything.

Particle masses and ``equalmass``
---------------------------------

The single most useful control over model quality is ``equalmass``.  Particle
mass is proportional to :math:`\rho^{1-\mathrm{equalmass}}`, so:

* ``equalmass=0`` (the default) gives constant number density.  Every part of
  the star gets a similar *number* of particles, so the tenuous outer layers are
  represented by enormous numbers of very light particles.
* ``equalmass=1`` gives particles of equal mass, concentrating them where the
  mass is and largely ignoring the outermost layers.

For a star spanning many orders of magnitude in density -- a giant, or a massive
main-sequence star -- the default is often a poor choice.  It can produce outer
particles of order :math:`10^{-10}\,M_\odot` alongside inner particles heavier
than :math:`10^{-2}\,M_\odot`, a spread that no relaxation will fix.

Values between 0.5 and 0.75 are a good starting point for such stars, and
sometimes 1 is better still.  The trade-off is explicit: raising ``equalmass``
abandons the outermost layers in exchange for a better model everywhere else,
and it slightly reduces the timestep, so the simulation runs more slowly.
Whether that is a good bargain depends on whether the outer envelope matters for
the science.  Values such as 0.6 or 0.7 are worth trying if neither extreme
looks right.

.. note::

   ``equalmass`` is the same quantity as :math:`\alpha` in Appendix A of
   `Gaburov et al. (2025), ApJ 980, 109
   <https://ui.adsabs.harvard.edu/abs/2025ApJ...980..109G/abstract>`_, where the
   particle-placement scheme is set out in full.

Neighbour number and ``nnopt``
------------------------------

``nnopt`` controls how many neighbours each particle has.  It must be an
integer.

As you increase the number of particles, increase ``nnopt`` as well -- but more
slowly, as :math:`\sqrt{N}`.  Doubling ``N`` means multiplying ``nnopt`` by
:math:`\sqrt{2}`, so a run using ``nnopt=23`` at ``N=100000`` would use
``nnopt=33`` at ``N=200000``.

Increasing ``N`` while leaving ``nnopt`` at a value inherited from a smaller run
is a common mistake.  The symptom is a model that overestimates pressure and
density in the centre and underestimates them through the bulk of the star.

Do not carry ``gam`` between stars
-----------------------------------

``gam`` is easy to mistake for a switch that only matters when the equation of
state is polytropic.  It is not.  Even with ``neos=1`` or ``neos=2`` it sets the
sound speed used by the artificial viscosity, in ``balAV3.f``:

.. math::

   c_i^2 = \mathrm{gam}\,\frac{P_i}{\rho_i},

written in the source as ``ci2 = gam*por2i*rho(i)``, since ``por2`` holds
:math:`P/\rho^2`.  The sound speed then sets the timestep.  So ``gam`` should reflect the star's
actual adiabatic index, not be inherited from whatever file you copied.

The value to use is a pressure-weighted :math:`\Gamma_1` for the star in
question.  A partially ionised giant envelope and a fully convective low-mass
star are different: values near 1.56 and 1.667 respectively are representative.
Copying one star's ``sph.input`` to another and leaving ``gam`` alone can be a
several per cent error in the sound speed, which will not announce itself --
the run proceeds, with a slightly wrong timestep and slightly wrong viscosity.

Cores and compact objects
-------------------------

A giant with a dense core is expensive to model particle by particle.  Two
settings help:

``mco``
   the mass of a core point particle.  Setting this forces a core of the mass
   you want, rather than accepting whatever the particle layout produces.

``hco``
   the softening length of that core particle.  Tuning it is a small but real
   improvement to the model.

Raising ``equalmass`` or ``N`` will generally make an automatically chosen core
point *less* massive, so if a specific core mass is needed, set ``mco``
explicitly.

A worked starting point
-----------------------

For a giant read from a MESA profile, a reasonable first ``sph.input`` is::

    &input
     n=100000,
     nnopt=32,
     nrelax=1,
     trelax=<one oscillation period, measured from a first run>,
     treloff=<a few trelax>,
     tf=<longer than treloff>,
     dtout=<so that you get a few dozen out files>,
     equalmass=0.5,
     neos=1,
     profilefile='profile.data',
    &end

Then iterate: relax, compare ``parent.sph`` against the ``col`` files, and adjust
``equalmass``, ``nnopt`` and ``N`` until the pressure and density profiles agree
and the two accelerations balance.
