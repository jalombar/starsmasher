Equations of motion
===================

StarSmasher does not use the textbook SPH equations.  The formulation it
integrates was derived in Appendix A of `Gaburov, Lombardi & Portegies Zwart
(2010), MNRAS 402, 105
<https://ui.adsabs.harvard.edu/abs/2010MNRAS.402..105G/abstract>`_
(`arXiv:0904.0997 <https://arxiv.org/abs/0904.0997>`_), and the difference
matters for exactly the problems this code is built for.

The problem with the usual constraint
-------------------------------------

Standard SPH ties the smoothing length to the local density,

.. math::

   h_i = f(\rho_i, C_i),

where the constant :math:`C_i` necessarily carries dimensions of mass.  In
effect each particle keeps a fixed *mass* inside its kernel.  That is
unobjectionable when every particle has the same mass, and it is the wrong
constraint as soon as they do not.

In a stellar collision, particles from two parent stars mix.  If the stars were
resolved differently -- or if ``equalmass`` was used, so that particle masses
already span orders of magnitude within a single star -- then a particle
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
target is the ``nnopt`` of :doc:`sph_input`: the code solves
``bonetfunc = N_i - nnopt = 0`` for :math:`h_i` at every step.

Because :math:`G` is a weighting function rather than a sharp cutoff, the
constraint is continuous, which is what makes it differentiable and so usable in
a variational derivation.  The nominal value ``nnopt=22`` yields roughly 35 to
40 particles actually inside the kernel.

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
particle.  This gives direct control over the smoothing length of particular
particles -- useful when compact objects are treated as point masses, which is
what Blackollider is for.  For everything else, ``parallel_bleeding_edge`` and
the neighbour-number constraint above are the ones to use.
