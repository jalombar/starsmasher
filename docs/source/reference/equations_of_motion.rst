Equations of motion
===================

StarSmasher does not use the textbook SPH equations.  The formulation it
integrates was derived in Appendix A of `Gaburov, Lombardi & Portegies Zwart
(2010), MNRAS 402, 105
<https://ui.adsabs.harvard.edu/abs/2010MNRAS.402..105G/abstract>`_
(`arXiv:0904.0997 <https://arxiv.org/abs/0904.0997>`_).  The difference matters
for exactly the problems this code is built for: collisions and mergers, where
particles of very different mass end up mixed together in the same fluid.

The standard constraint
-----------------------

Standard SPH ties the smoothing length to the local density,

.. math::

   h_i = f(\rho_i, C_i),

where :math:`C_i` must carry dimensions of mass for :math:`h_i` to come out as
a length.  Written the usual way, :math:`h_i = \eta\,(m_i/\rho_i)^{1/3}`, that
mass is the particle's own.  In effect each particle keeps a fixed *mass*
inside its kernel.  That is unobjectionable when every particle has the same
mass, and it is arguably the wrong constraint as soon as they do not.

In a stellar collision, particles from two parent stars mix.  If the stars were
resolved differently, or if ``equalmass`` was used so that particle masses
already span orders of magnitude within a single star, then a particle
drifting into a region where the mean particle mass differs will find itself
with far too few neighbours, or far too many.  Neither is recoverable: too few
and the density estimate is noise, too many and the resolution has quietly been
thrown away.

The constraint StarSmasher uses instead
---------------------------------------

Rather than fixing the mass inside the kernel, StarSmasher fixes an estimate of
the *number of neighbours*.  Each particle carries

.. math::

   N_i = \sum_j G\!\left(|\mathbf{r}_i - \mathbf{r}_j|,\, h_i\right),

and :math:`h_i` is solved for so that :math:`N_i` equals a target value.  That
target is the ``nnopt`` of :doc:`sph_input`.

:math:`G` is best understood by imagining a cruder version of itself.  Suppose
it were a step function, equal to 1 for a particle inside the kernel and 0
outside.  Then :math:`N_i` would be nothing more than a count of particle
:math:`i`'s neighbours, and the constraint would read "give every particle the
same number of neighbours".

:math:`G` is that count made smooth.  It falls off gradually instead of
switching, so a particle drifting towards the kernel edge contributes a
diminishing fraction rather than vanishing the instant it crosses.  That is
what makes the constraint differentiable, and so usable in a variational
derivation, at the cost of :math:`N_i` no longer being an integer count of
anything.

Which :math:`G` is used is set by ``gflag``.  ``gflag=0`` selects the function
of equation (A2) and figure A1 of the paper.  ``gflag=1``, the default, is the
same except that :math:`G=1` inside the inner half of the kernel, which works
better when massive compact objects or core particles are present.  The
namelist default ``nnopt=22+gflag`` follows ``gflag``, so changing one moves
the other.

Smoothing lengths are updated at every timestep, so this is a root find per
particle per step.  Standard SPH solves for :math:`h_i` too, so the expense is
not particular to this formulation.

How many neighbours you actually get
------------------------------------

``nnopt`` is the target of the constraint, not the number of neighbours.  The
number of particles actually inside the kernel is larger, by a factor of roughly
1.4 to 1.7.

The factor is largest at small ``nnopt`` and settles to roughly 1.4 to 1.45
above about ``nnopt=60``.  Measurements across a polytrope, a low-mass MESA
model and a red giant span 1.38 to 1.66, and Gaburov et al. quote ``nnopt=22``
giving 35 to 40 neighbours.

The density estimate is unchanged, and remains the ordinary kernel sum

.. math::

   \rho_i = \sum_j m_j W\!\left(|\mathbf{r}_i - \mathbf{r}_j|,\, h_i\right).

What this costs
---------------

Because :math:`h_i` now depends on particle positions through the neighbour
constraint rather than through the density, the equation of motion acquires
correction terms of a different form from the familiar :math:`\nabla h` terms.
Alongside the usual gradient of :math:`W`, each pair contributes a gradient of
:math:`G`, weighted by factors (written :math:`\omega` and :math:`\chi` in the
paper) that measure how the density and neighbour-number sums respond to a
change in smoothing length.  The result is still derived from a Lagrangian, so
energy and momentum conservation are retained.

The price is smaller than that description suggests.  Standard SPH already
sweeps each particle's neighbours to build :math:`\rho_i`, and already carries
a :math:`\nabla h` correction of its own.  This formulation adds two more
running totals to that same sweep, one for :math:`G` and one for its derivative
with respect to :math:`h`, together with a second table lookup for each pair
because those two are evaluated on a slightly different length than the density
is.  There is no additional pass over pairs and no second neighbour search.

Why it is worth it
------------------

Particles of very different mass can mix freely while each retains a sensible
number of neighbours.  That is precisely what a collision or merger demands, and
it is what makes ``equalmass`` a usable knob: a model can concentrate particles
where the mass is, spanning a wide range of particle masses, without any region
becoming under- or over-resolved as the star is disrupted.

.. note::

   This is also why ``nnopt`` behaves the way :doc:`../tutorials/relaxing_a_star`
   describes.  It is not a loose suggestion about neighbour counts but the
   target of an equation solved for every particle at every step, which is why
   leaving it unchanged while increasing ``N`` degrades the model.

Blackollider
------------

The secondary code in ``Blackollider/`` relates :math:`h` and :math:`\rho`
differently.  Instead of driving the neighbour-number estimate to a single
global ``nnopt``, it solves a per-particle relation of the form

.. math::

   h_i = d_i + \left(\frac{1}{a_i} + b_i\, N_i^{1/3}\right)^{-1},

where :math:`a_i`, :math:`b_i` and :math:`d_i` are properties of the individual
particle.  Because :math:`d_i` is added on, it acts as a floor: particle
:math:`i` can never have a smoothing length below it, whatever the neighbour
count does.

That is what the code is for.  In ``parallel_bleeding_edge`` it can be hard to
stop particles piling up around a massive compact object, and the very small
smoothing lengths that follow bring numerical trouble with them.  A per-particle
lower limit on :math:`h` is a convenient way to prevent that, and this relation
provides one.

For everything else, ``parallel_bleeding_edge`` and the neighbour-number
constraint above are the ones to use.  They are what these pages document and
what the tutorials are written against.
