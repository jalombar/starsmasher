Input files
===========

Every run reads two files from the directory it is launched in.  Depending on
what you are setting up, it may read others as well.

Always needed
-------------

``sph.init``
~~~~~~~~~~~~

A one-variable namelist naming the initialization script, and nothing else::

    &INITT
    INAME='erg'
    &END

See :doc:`running` for what the codes mean, and :doc:`../reference/sph_init`
for the full list with the routine each one calls.

``sph.input``
~~~~~~~~~~~~~

Everything else: how long to run for, how often to write, how many particles,
which equation of state, the Courant numbers, the orbit.  Every setting has a
default, so the file only has to contain what you want to change::

    &input
    tf=1000.0,
    dtout=10.0,
    n=8000,
    &end

A setting you do not name keeps its default, and a setting the chosen
initialization script does not read is ignored without complaint.
:doc:`../reference/sph_input` lists all of them with their defaults.

.. _code-units:

Code units
----------

Numbers in ``sph.input`` are in code units unless a setting says otherwise.
StarSmasher works in units where

.. math::

   G = M_\mathrm{unit} = R_\mathrm{unit} = 1.

Three settings fix what that means in cgs, and all three can be reset in
``sph.input`` like any other::

    &input
    runit=6.957d10,
    munit=1.9884098706980504d33,
    gravconst=6.67430d-08,
    &end

``runit`` is the number of centimetres in the unit of length, ``munit`` the
number of grams in the unit of mass, and ``gravconst`` Newton's constant in
cgs.  The values above are the defaults: the MESA solar radius, the MESA solar
mass, and the 2018 CODATA value of :math:`G`.  Out of the box, then, a mass of
1 is a solar mass and a length of 1 is a solar radius.

Every other unit is built from those three and nothing else, so changing any
one of them rescales all of them.  The defaults give the values in the last
column.

.. list-table::
   :header-rows: 1
   :widths: 26 36 38

   * - Quantity
     - Unit
     - Default value in cgs
   * - mass
     - :math:`M_\mathrm{unit}`
     - approximately :math:`1.99\times10^{33}` g
   * - length
     - :math:`R_\mathrm{unit}`
     - approximately :math:`6.96\times10^{10}` cm
   * - time
     - :math:`\sqrt{R_\mathrm{unit}^3 / (G M_\mathrm{unit})}`
     - approximately :math:`1.59\times10^{3}` s, or 26.5 minutes
   * - velocity
     - :math:`\sqrt{G M_\mathrm{unit} / R_\mathrm{unit}}`
     - approximately :math:`4.37\times10^{7}` cm s\ :sup:`-1`, or 437 km s\ :sup:`-1`
   * - acceleration
     - :math:`G M_\mathrm{unit} / R_\mathrm{unit}^2`
     - approximately :math:`2.74\times10^{4}` cm s\ :sup:`-2`
   * - density
     - :math:`M_\mathrm{unit} / R_\mathrm{unit}^3`
     - approximately :math:`5.91` g cm\ :sup:`-3`
   * - pressure
     - :math:`G M_\mathrm{unit}^2 / R_\mathrm{unit}^4`
     - approximately :math:`1.13\times10^{16}` dyn cm\ :sup:`-2`
   * - energy
     - :math:`G M_\mathrm{unit}^2 / R_\mathrm{unit}`
     - approximately :math:`3.79\times10^{48}` erg
   * - specific energy
     - :math:`G M_\mathrm{unit} / R_\mathrm{unit}`
     - approximately :math:`1.91\times10^{15}` erg g\ :sup:`-1`
   * - angular momentum
     - :math:`\sqrt{G M_\mathrm{unit}^3 R_\mathrm{unit}}`
     - approximately :math:`6.04\times10^{51}` g cm\ :sup:`2` s\ :sup:`-1`

Put less formally: with the defaults, a density of 1 is one solar mass spread
through one cubic solar radius, and a time of 1 is how long a body in a
circular orbit grazing the solar surface takes to sweep out one radian, a
little under half an hour.  That last one is exact rather than approximate,
which is why a whole such orbit takes :math:`2\pi` code units of time, or
about 2.78 hours.

Temperature is the one quantity that is never scaled.  It is in kelvin
everywhere, in the input and in the output alike.

A polytrope is the exception to all of this, and only when ``neos=0``.  The
polytropic equation of state is :math:`P = A\rho^\gamma`, which brings in no
physical constant, so nothing in the calculation refers to grams or centimetres
at all.  A polytrope of ``starmass=1`` and ``starradius=1`` is a star of one
mass unit and one radius unit, not one solar mass and one solar radius, and you
may read those units as whatever you like.  Setting ``neos=1`` or ``neos=2``
brings physical constants back in, and the model becomes a star of a definite
size again.

Depending on what you are setting up
------------------------------------

These are named by settings in ``sph.input`` rather than by fixed file names, so
each one can be called whatever you like.  All are collected under
:ref:`Input files <sph-input-input-files>` in the ``sph.input`` reference.

A stellar-evolution profile
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``profilefile`` names a plain-text profile of the star you want to reproduce,
a MESA model for instance.  ``erg`` reads it, and without it there is nothing
to build a star from.  ``stellarevolutioncodetype`` says which code wrote it, since
the column layouts differ.

The profile is also written back out as ``parent.sph`` for comparison, which is
how you check that the SPH star matches the model it came from.  See
:doc:`output`.

One or two relaxed stars
~~~~~~~~~~~~~~~~~~~~~~~~

``startfile1`` names the first body of the encounter and ``startfile2`` the
second.

These are not a special format: each is an ordinary ``out*.sph`` snapshot, and
in practice it is the last one a relaxation wrote, that being the relaxed star.
Making a start file means copying that snapshot out of the relaxation directory
and into a fresh one for the collision, under the name the setting expects::

    $ mkdir collision
    $ cp relax_star1/out0042.sph collision/sph.start1u
    $ cp relax_star2/out0038.sph collision/sph.start2u

This is why a collision is normally the *second* thing you run: the first run
makes the star, the second collides it.

.. warning::

   Run the collision in a directory of its own rather than reusing the
   relaxation's.  A relaxation leaves a ``restartrad.sph`` behind, and the code
   picks that up automatically, so a collision started in the same directory
   quietly resumes the relaxation instead.

Seven initialization scripts read ``startfile1``:

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - Code
     - Uses it as
   * - ``hyp``
     - the first of two bodies on a Keplerian orbit
   * - ``bps``
     - a member of the binary, with ``startfile2`` as the single star
   * - ``bph``
     - a member of the binary meeting a black hole
   * - ``tri``
     - one body of the triple
   * - ``2cr``
     - one star of the corotating binary
   * - ``bhe``
     - the star approaching the supermassive black hole
   * - ``res``
     - the model being rescaled

``startfile2`` is the same thing for the second body, read by all of those
except ``res``, which rescales a single model.  ``startfile3`` is the third body
of a triple, read only by ``tri``.

.. note::

   ``hyp`` treats a missing ``startfile2`` as meaning "no second star": the
   second object becomes a single point mass of mass ``mbh``, with softened
   gravity.  Giving ``bbh_m2`` a positive mass asks for a compact-object binary
   there instead.


What a start file decides for you
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A start file carries the model's structure with it, so two settings in
``sph.input`` are read from the file rather than from what you wrote:

``n``
   particle number.  You choose it when you relax a star.  After that it belongs
   to the model, and a collision simply uses however many particles its start
   files contain.  Two copies of a 1974-particle star make a 3948-particle
   collision, whatever ``sph.input`` asks for.

``nnopt``
   taken from ``startfile1``.  ``log0.sph`` records the substitution::

       NOTE: Currently nnopt=          77
             The NNOPT in startfile1=          22
             Changing NNOPT to be          22

Both are decisions made during the relaxation, and neither can be revisited
later.  Changing the resolution of a model means relaxing it again.

.. warning::

   Every body in an encounter must have been relaxed at the same ``nnopt``.
   The run stops if they were not::

       ERROR: nnopt from star 1=          22
              nnopt from star 2=          30

   This matters when the two stars are very different: it is tempting to relax
   a giant and a dwarf at values suited to each, and they will then refuse to
   be collided.  Choose the ``nnopt`` you intend to collide at before relaxing
   either of them.

Setup description files
~~~~~~~~~~~~~~~~~~~~~~~

Some scripts need a short text file describing the arrangement, separate from
the bodies themselves.  Each is named by a setting, so several calculations can
share a directory:

.. list-table::
   :header-rows: 1
   :widths: 20 16 64

   * - Setting
     - Default
     - Read by
   * - ``binaryfile``
     - ``input.bs``
     - ``bps`` and ``2cr``: the binary's masses and separation
   * - ``triplefile``
     - ``input.3s``
     - ``bhe`` and ``tri``: each body's mass, offset and velocity
   * - ``bpbhfile``
     - ``sph.bpbh``
     - ``bph``: orientation angles and the black hole mass
   * - ``imagefile``
     - ``sph.image``
     - ``txt``: an ASCII picture turned into a particle layout

The defaults are the names these scripts used to have built in, so a directory
that worked before still works untouched.

Physics tables
~~~~~~~~~~~~~~

``eosfile`` and ``opacityfile`` name tabulated data, read only when the settings
that need them are switched on, such as a tabulated equation of state or
cooling.
