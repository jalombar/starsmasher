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

Set ``trelax=0`` and let the code work it out::

    trelax=0,

The drag damps the star's fundamental radial oscillation fastest when its
timescale is near that oscillation's period, and the period scales with the
dynamical time :math:`t_\mathrm{dyn} = \sqrt{R^3/M}`.  So from the star's own
radius and mass the code sets

.. math::

   t_\mathrm{relax} = 4\,t_\mathrm{dyn}, \qquad
   t_\mathrm{reloff} = 10\,t_\mathrm{relax}.

The coefficient 4 is measured, not assumed: on a 5.4 :math:`M_\odot` AGB giant
the period comes out at 3.91 :math:`t_\mathrm{dyn}` and on a 0.4
:math:`M_\odot` M dwarf 4.43, two stars a factor of a thousand apart in radius.
A ten per cent error in ``trelax`` moves the relaxed model's residual motion by
only two or three per cent, which is why one constant serves both.

The ten is an e-folding count: the oscillation envelope decays with a time
constant equal to ``trelax`` itself, so ten of them is ten e-foldings before the
drag is switched off.

``log0.sph`` records what it chose::

    relax: trelax=0, so the drag schedule is derived:
    relax:   radius  =   ...
    relax:   mass    =   ...
    relax:   t_dyn   =   ...
    relax:   trelax  =   ...
    relax:   treloff =   ...

.. important::

   ``trelax=0`` sets ``treloff`` as well, overwriting whatever you put there.
   You cannot derive one and hand-set the other.  It is both or neither.

When to override it
~~~~~~~~~~~~~~~~~~~

Setting ``trelax`` to a positive value keeps its old meaning exactly, and then
``treloff`` is yours to set too.  Two reasons to do so:

**A very soft envelope.**  Where :math:`\Gamma_1` approaches 4/3 the period
lengthens sharply, and :math:`4\,t_\mathrm{dyn}` will underestimate it.  The
symptom is a model still ringing when the drag comes off.

**You want to measure the period yourself.**  Relax briefly with
``trelax=1.d30``, which disables the drag, and take the dominant frequency of
the kinetic energy.  It oscillates at *twice* the mode frequency, so the period
is twice what you read off.  Then set ``trelax`` to that period and ``treloff``
to ten times it.

.. note::

   On a resumed run the schedule is derived again from the star as it now is,
   which has already relaxed and so is slightly smaller than the one it started
   from.  ``log0.sph`` says so when this happens.

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

For a star spanning many orders of magnitude in density, such as a giant or a
massive main-sequence star, the default is often a poor choice.  It can produce outer
particles of order :math:`10^{-10}\,M_\odot` alongside inner particles heavier
than :math:`10^{-2}\,M_\odot`, a spread that no relaxation will fix.

Choosing a value
~~~~~~~~~~~~~~~~

There is no formula.  The choice follows from the density contrast of the star
and from which part of it your science depends on, and it is settled by trying
a value and looking at the result.

**Main-sequence stars and polytropes: start at 0**, the default.  The density
range is modest enough that constant number density resolves the outer layers as
well as the core, which is what you want if the model is headed for a
mass-transfer binary, where the envelope is the part that matters.

**Giants: start at 0.5.**  A giant's core-to-envelope contrast is large enough
that resolution has to be bought somewhere, and buying it in the core is usually
right.  Values of about 0.4 and above are the useful range.

The reason to avoid going low on a giant is a specific failure rather than a
gradual loss of quality.  An under-resolved core cannot hold itself in
equilibrium, and a shock then travels outwards from it: the relaxation does not
settle and the model is unusable.

.. note::

   Some outward motion early on is normal, and is not the same thing.  What
   distinguishes a failure is whether the gas comes back: a model that moves
   out and settles is fine, one that keeps going has effectively exploded.
   Judge this after ``treloff``, not during the drag.  See
   `Judging a model after release`_.

Raising ``equalmass`` also slightly reduces the timestep, so the run is slower.
Values such as 0.6 or 0.7 are worth trying when neither 0.5 nor 1 looks right.

.. note::

   ``equalmass`` is the same quantity as :math:`\alpha` in Appendix A of
   `Gaburov et al. (2025), ApJ 980, 109
   <https://ui.adsabs.harvard.edu/abs/2025ApJ...980..109G/abstract>`_, where the
   particle-placement scheme is set out in full.

Neighbour number and ``nnopt``
------------------------------

``nnopt`` controls how many neighbours each particle has.  It must be an
integer.

As you increase the number of particles, increase ``nnopt`` as well, but more
slowly, as :math:`\sqrt{N}`.  Doubling ``N`` means multiplying ``nnopt`` by
:math:`\sqrt{2}`, so a run using ``nnopt=23`` at ``N=100000`` would use
``nnopt=33`` at ``N=200000``.

Increasing ``N`` while leaving ``nnopt`` at a value inherited from a smaller run
is a common mistake.  The symptom is a model that overestimates pressure and
density in the centre and underestimates them through the bulk of the star.

The rule does not extend downwards indefinitely.  Every SPH quantity is a sum
over neighbours, and a sum over few neighbours carries shot noise that no amount
of relaxing will remove, so there is a floor below which the scaling stops
meaning anything.  Around ``nnopt = 22 + gflag`` is a reasonable one to treat as
the minimum.

.. note::

   That floor is a rule of thumb rather than a measured threshold: it has not
   been tested, and it is offered as a sensible place to stop scaling rather
   than a value with evidence behind it.  If you have reason to go lower,
   compare the resulting model against ``parent.sph`` before trusting it.

Choose ``nnopt`` before you relax anything
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``nnopt`` is fixed at relaxation and travels with the model.  A collision reads
it from ``startfile1`` and overrides whatever ``sph.input`` says, and it stops
outright if the two stars disagree.

The trap is a collision between very different stars, where it is natural to
relax a giant and a dwarf at values suited to each, and then find they cannot
be collided at all.  Decide the value you intend to collide at first, and relax
every body at that value.  See :doc:`../using/input`.

Cores and compact objects
-------------------------

A giant with a dense core is expensive to model particle by particle.  Two
settings help:

``mco``
   the mass of the core point particle.  Setting it forces a core of the mass
   you want.  Leaving it at zero has the code choose one.

``hco``
   the gravitational softening length of that core particle.  Tuning it is a
   small but real improvement to the model.

.. important::

   An automatically chosen ``mco`` is **not** a core mass taken from the
   stellar-evolution model.  The code lays down SPH particles, finds their total
   falls short of the mass the star is supposed to have, and gives the core
   point exactly the difference::

       am(1) = amass - amtot

   It is a bookkeeping value that makes the total come out right, given whatever
   mass the lattice failed to place near the centre.  Read as a physical core
   mass it will mislead you, and it has no reason to agree with the core mass
   your MESA model reports.

   This is also why raising ``equalmass`` or ``N`` makes the automatic core
   *less* massive: both put more particle mass near the centre, so less is left
   over.  If a specific core mass is what you need, set ``mco`` yourself.

Choosing ``mco``, ``equalmass`` and ``N`` together is trial and error.  There is
no settled procedure, and the three interact.

A worked starting point
-----------------------

For a giant read from a MESA profile, a reasonable first ``sph.input`` is::

    &input
     n=100000,
     nnopt=32,
     nrelax=1,
     trelax=0,
     tf=<longer than the treloff the code reports>,
     dtout=<so that you get a few dozen out files>,
     equalmass=0.5,
     neos=1,
     profilefile='profile.data',
    &end

``trelax=0`` sets ``treloff`` too, so it is left out.  Read the value the code
chose out of ``log0.sph`` and make ``tf`` comfortably longer than it, or the run
ends before the star has been released long enough to judge.

Then iterate: relax, compare ``parent.sph`` against the ``col`` files, and adjust
``equalmass``, ``nnopt`` and ``N`` until the pressure and density profiles agree
and the two accelerations balance.
