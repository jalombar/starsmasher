StarSmasher
===========

StarSmasher is a smoothed particle hydrodynamics (SPH) code for modelling
collisions and close interactions between stars and planets.  It is
parallelised with MPI and computes self-gravity either on NVIDIA GPUs or, where
no such card is available, on the CPU.

The code descends from StarCrash, written by Joshua Faber, James Lombardi and
Frederic Rasio, and is maintained at
`github.com/jalombar/starsmasher <https://github.com/jalombar/starsmasher>`_.

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
