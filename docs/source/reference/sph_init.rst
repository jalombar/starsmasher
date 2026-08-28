sph.init
========

``sph.init`` is a one-variable namelist that chooses what kind of calculation
to set up.  It is read before ``sph.input``::

    &INITT
    INAME='1es'
    &END

``INAME`` is a three-letter code.  Each one selects a different routine in
``parallel_bleeding_edge/src``:

.. list-table::
   :header-rows: 1
   :widths: 12 30 58

   * - Code
     - Routine
     - Sets up
   * - ``1es``
     - ``initialize_polyes.f``
     - A single polytrope, generated from ``starmass`` and ``starradius``.
   * - ``1mc``
     - ``initialize_polymces.f``
     - A single polytrope with a compact core.
   * - ``2cr``
     - ``initialize_corotating.f``
     - A corotating binary.
   * - ``bps``
     - ``initialize_bps.f``
     - A binary plus a single star.
   * - ``bph``
     - ``initialize_bpbh.f``
     - A binary plus a black hole.
   * - ``hyp``
     - ``initialize_hyperbolic.f``
     - Two relaxed stars on a Keplerian orbit: the usual choice for a
       collision or fly-by.
   * - ``hbs``
     - ``initialize_hyperbolic_binary_single.f``
     - A binary encountering a single star.
   * - ``erg``
     - ``initialize_parent.f``
     - A star read from a stellar-evolution profile (MESA, TWIN or
       Freitag & Benz).
   * - ``meq``
     - ``initialize_multiequalmass.f``
     - A star whose particle-mass scheme varies with position.
   * - ``res``
     - ``initialize_rescale.f``
     - A rescaling of an existing model.
   * - ``tri``
     - ``initialize_triple.f``
     - A triple system.
   * - ``bhe``
     - ``initialize_smbh.f``
     - An encounter with a supermassive black hole.
   * - ``grs``
     - ``initialize_grsph.f``
     - A general-relativistic SPH setup.
   * - ``rin``
     - ``init.f`` (``readin``)
     - A restart from an existing dump.
   * - ``txt``
     - ``initialize_asciiimage.f``
     - Particles laid out from an ASCII image.

.. note::

   Which ``INAME`` you choose decides which settings in :doc:`sph_input`
   actually do anything.  ``vinf2`` and ``e0``, for instance, are read only by
   the orbit routines, and ``equalmass`` only by the routines that build a star
   particle by particle.
