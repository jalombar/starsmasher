StarSmasher
===========

StarSmasher is a smoothed particle hydrodynamics (SPH) code for modelling
collisions and close interactions between stars and planets.  It is
parallelised with MPI and computes self-gravity either on NVIDIA GPUs or, where
no such card is available, on the CPU.

The code descends from StarCrash, written by Joshua Faber, James Lombardi and
Frederic Rasio, and is maintained at
`github.com/jalombar/starsmasher <https://github.com/jalombar/starsmasher>`_.

Where to start
--------------

:doc:`quickstart` goes from a clone to a relaxed star in a handful of commands,
and is the fastest way to find out whether the code builds on your machine.

Work then falls into two stages, and the tutorials follow that order.  You relax
a star into an SPH model, once per star, and then you collide the results.
:doc:`tutorials/relaxing_a_polytrope` does the first with a polytrope,
:doc:`tutorials/relaxing_a_star` with a stellar-evolution profile, and
:doc:`tutorials/a_collision` does the second.

:doc:`reference/sph_input` lists every setting with its default and says which
initialization scripts read it.  If you are looking for what a particular
number in an output file means, that is :doc:`using/output`.

If the physics is what you are after rather than the mechanics,
:doc:`reference/equations_of_motion` explains what StarSmasher does differently
from textbook SPH and why it matters for collisions.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   quickstart
   installation

.. toctree::
   :maxdepth: 2
   :caption: Using StarSmasher

   using/index
   tutorials/index

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reference/index

.. toctree::
   :maxdepth: 2
   :caption: Project

   test_suite
   about
   contributing
   faq
