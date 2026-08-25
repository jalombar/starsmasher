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

Three further columns appear conditionally, which is worth knowing before
writing a parser:

* with radiative cooling (``ncooling`` non-zero), the radiated energy follows;
* for a binary, the orbital angular frequency and the separation of the two
  centres of mass follow.

So a single-star run without cooling has seven columns and a cooling binary has
ten.  Count them rather than assuming.

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
