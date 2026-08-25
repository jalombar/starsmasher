Output files
============

A run writes several kinds of file.  All are numbered by *stage*, not by rank: a
fresh run writes ``log0.sph`` and ``energy0.sph``, the first resume writes
``log1.sph`` and ``energy1.sph``, and so on.  A resumed run does not append to
the earlier file.

``energy0.sph``
---------------

One row per output, written by ``enout`` in ``output.f``.  The first seven
columns are always present:

.. list-table::
   :header-rows: 1
   :widths: 12 20 68

   * - Column
     - Symbol
     - Quantity
   * - 1
     - :math:`t`
     - Time, in code units.
   * - 2
     - :math:`W`
     - Gravitational potential energy.
   * - 3
     - :math:`T`
     - Kinetic energy.
   * - 4
     - :math:`U`
     - Internal energy.
   * - 5
     - :math:`E`
     - Total energy, the sum of the three above.
   * - 6
     - :math:`S`
     - Total entropy.
   * - 7
     - :math:`J`
     - Total angular momentum.

Further columns appear conditionally, and it is worth knowing which before
writing a parser.  Three separate settings decide the layout:

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Condition
     - Columns
     - Extra fields, in order
   * - the common case
     - 7
     - none
   * - ``ncooling`` non-zero
     - 8
     - radiated energy
   * - ``reat`` positive
     - 10
     - the potential, kinetic and internal energy of eaten material
   * - both of the above
     - 11
     - radiated energy, then the three eaten-material energies
   * - ``nrelax`` >= 2
     - 9
     - orbital angular frequency, then the separation of the two centres of mass

The last of these is the one that catches people.  The angular frequency and
separation are written for a *corotating-frame relaxation*, which is what
``nrelax >= 2`` selects -- not for a dynamical collision.  A collision run has
``nrelax=0`` and gets the ordinary seven columns even though there are plainly
two stars in it.

Count the columns rather than assuming.

Parsing an energy trace
~~~~~~~~~~~~~~~~~~~~~~~

Two habits will save you a wrong answer.

**Count the columns; do not assume an index.**  The first seven are fixed, but
anything beyond them depends on ``ncooling``, ``reat`` and ``nrelax``, as above.
Code that reaches for column 8 will read the radiated energy in one run and the
orbital angular frequency in another, with nothing to announce the difference.

**Look for later stages.**  A resumed run writes ``energy1.sph``, then
``energy2.sph``, rather than appending to ``energy0.sph``.  Reading only
``energy0.sph`` after a restart silently gives you the first stage and stops,
which looks exactly like a run that ended early.

Reading a relaxation from this file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Column 5 should be nearly constant: it is the cleanest check that a run is
behaving.  Column 3 should fall to something negligible compared with column 4,
which is what a relaxed star looks like.

To measure an oscillation period, use column 2 rather than column 3.  The
gravitational potential energy oscillates once per cycle about a slowly varying
mean, and its zero crossings give a clean period.  The kinetic energy oscillates
at twice the mode frequency with a rapidly decaying envelope, and reading a
period from it is unreliable -- the answer depends on how the data were
detrended.

.. note::

   Judge a relaxation *after* the drag has been switched off at ``treloff``, not
   during it.  Under drag everything looks quiet, and two configurations that
   appear identical while relaxing can differ substantially once released.  Allow
   several oscillation periods after release before drawing a conclusion: an
   effect that has not appeared yet is indistinguishable from one that does not
   exist.

``out*.sph``
------------

Binary snapshots of the full particle state, written every ``dtout``.  These are
what you visualise, and what you use as the starting point for a collision: the
last snapshot of a relaxation is a relaxed star.

``col*.sph``
------------

Written during relaxations.  Comparing ``col0000.sph`` against ``parent.sph``
shows how well the initial particle layout matched the profile it was asked to
reproduce; comparing a later one shows what the relaxation did with it.  See
:doc:`../tutorials/relaxing_a_star`.

``parent.sph``
--------------

The profile the code was asked to reproduce, written when a star is built from a
stellar-evolution model.  It is the reference against which the SPH model should
be judged.

``restartrad.sph``
------------------

Written every few iterations, overwriting the previous one.  If present when a
run starts, it is picked up automatically.  Any ``out*.sph`` can be renamed to
``restartrad.sph`` to restart from that point instead.

``log*.sph``
------------

The running commentary: what was read, what was decided, the neighbour
statistics, the timestep, and the energies at each output.  It is the first
place to look when a run does something unexpected.
