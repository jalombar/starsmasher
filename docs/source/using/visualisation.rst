Looking at the output
=====================

The snapshots are binary, so you need something that understands the format.
The usual choice is `SPLASH <https://github.com/danieljprice/splash>`_.

Installing SPLASH
-----------------

Two guides in the repository cover this, and there is no point repeating them:

* `Data_visualization/Readme.md
  <https://github.com/jalombar/starsmasher/blob/master/Data_visualization/Readme.md>`_
  is a step-by-step walkthrough written for StarSmasher users, including
  building giza, which SPLASH needs and does not ship.
* Section 3 of `documentation/installation.md
  <https://github.com/jalombar/starsmasher/blob/master/documentation/installation.md>`_
  covers a system-wide install and the ``libgiza.so`` problems that follow one.

Reading StarSmasher snapshots
-----------------------------

SPLASH has a StarSmasher reader built in::

   splash -f starsmasher out*.sph

If your SPLASH predates it, ``Data_visualization/read_data_jamiesph.f90`` is the
reader, with a ``Makefile`` beside it.

The examples in ``example_input`` ship ``splash.defaults`` and ``splash.units``.
Copy them into your run directory before starting SPLASH and it will open with
sensible axes and units rather than raw code units.

What to look at
---------------

For a **relaxation**, the useful comparison is not a picture but a profile.
Plot pressure, density, particle mass, smoothing length and neighbour number
against radius, for ``parent.sph`` and for ``col0000.sph`` together: that shows
how well the initial particle layout matched the profile it was asked to
reproduce.  Repeat with a later ``col`` file to see what the relaxation did.  The
radial components of the hydrodynamic and gravitational accelerations are worth
adding: in a settled model they mirror each other, and where they do not, the
star is not in equilibrium.  See :doc:`../tutorials/relaxing_a_star`.

For a **collision**, the snapshots are the point.  Column density in the orbital
plane over successive ``out`` files shows the encounter, and ``energy0.sph``
tells you whether to trust it.  See :doc:`output`.

.. note::

   A picture can look convincing while the model is poor.  Check the profiles
   and the energy trace as well as the image.

Other tools
-----------

``splot_routines`` in the repository contains analysis scripts, including
routines for masses, orbits and spins.
