Output files
============

Every run writes these:

``log*.sph``
   How the run was launched, then a summary of the system at each timestep.  The
   first place to look for clues when something goes wrong.

``energy*.sph``
   Energies and related quantities, one row per output, with the time in the
   first column.

``out*.sph``
   The run settings and the full particle state, written every ``dtout``.
   Unformatted, so not readable in a text editor.  These are what you visualise,
   and what a collision starts from.

``restartrad.sph``
   The same format as ``out*.sph``, rewritten every 1000 iterations as a
   checkpoint.  If it is present in the directory when a run starts, the code
   continues from it rather than beginning afresh.

A relaxation writes two more, as a convenience for checking the model:

``parent.sph``
   The profile the code was asked to reproduce, taken from the desired stellar
   model such as a MESA one.  The reference the SPH star should be judged
   against.

``col*.sph``
   Plain text, one line per particle, holding pressure, density, mass, smoothing
   length, neighbour number, and the gravitational and hydrodynamic
   accelerations, among others.  Several of these can be compared directly
   against ``parent.sph``.

All are numbered by *stage*, not by rank: a fresh run writes ``log0.sph`` and
``energy0.sph``, the first resume writes ``log1.sph`` and ``energy1.sph``, and so
on.  A resumed run does not append to the earlier file.

Written every time
------------------

``log*.sph``
~~~~~~~~~~~~

Opens with a record of how the run was set up: the settings it read, the choices
it made from them, and anything it decided for you.  It is the only place some of
those decisions are written down, so it is worth reading once at the start of a
run rather than only when something has gone wrong.

After that it settles into a per-timestep commentary.  When a run misbehaves,
this is the file that says what the code thought it was doing at the time.

.. dropdown:: The timestepping lines, and what to do when energy is not conserved
   :icon: clock

   Every step writes a pair of lines::

      dts= <dt1> <dt2> <dt3> <dt4> <dt5> <dt6> <dt>
      indx <i1>  <i2>  <i3>  <i4>  <i5>  <i6>

   ``dts=`` gives the smallest timestep each criterion would have allowed, and
   then ``dt``, the step actually taken.  ``indx`` gives the particle
   responsible for each of those minima, which is what turns "the timestep
   collapsed" into "the timestep collapsed because of *this* particle".

   Four criteria apply to SPH particles and are combined as
   :math:`dt_{sph} = 1/(1/dt_1 + 1/dt_2 + 1/dt_3 + 1/dt_4)`:

   .. list-table::
      :header-rows: 1
      :widths: 8 34 20 38

      * - 
        - Limits the step by
        - Default
        - 
      * - ``dt1``
        - :math:`cn_1 h / v_{signal}`
        - ``cn1=0.3``
        - the Courant condition: how far a signal crosses a kernel
      * - ``dt2``
        - :math:`cn_2 (h / |a - a_{sm}|)^{1/2}`
        - ``cn2=1d30``
        - acceleration across a smoothing length
      * - ``dt3``
        - :math:`cn_3 u / |du/dt|`
        - ``cn3=0.1``
        - how fast a particle's internal energy is changing
      * - ``dt4``
        - :math:`cn_4 v_{signal} / |a - a_{sm}|`
        - ``cn4=1d30``
        - acceleration against signal speed

   Two more apply to a particle that is a compact object, combined as
   :math:`dt_{co} = 1/(1/dt_5 + 1/dt_6)`, each minimised over every other
   particle *j*:

   .. list-table::
      :header-rows: 1
      :widths: 8 34 20 38

      * - 
        - Limits the step by
        - Default
        - 
      * - ``dt5``
        - :math:`cn_5 r_{ij} / v_{ij}`
        - ``cn5=0.02``
        - separation over relative velocity
      * - ``dt6``
        - :math:`cn_6 (r_{ij} / a_{ij})^{1/2}`
        - ``cn6=0.02``
        - separation over relative acceleration
      * - ``cn7``
        - :math:`r_{ij} = (x_{ij}^2 + y_{ij}^2 + z_{ij}^2 + cn_7 h_i^2)^{1/2}`
        - ``cn7=4``
        - softening inside the separation used above

   The step finally taken is the smallest of :math:`dt_{sph}` and
   :math:`dt_{co}` over all particles.

   Two things follow that are not obvious from the numbers alone.  ``cn2`` and
   ``cn4`` default to ``1d30``, which makes ``dt2`` and ``dt4`` enormous and
   their reciprocals negligible: those two criteria are effectively switched
   off unless you turn them on.  And ``dt`` is smaller than every individual
   column, because reciprocals add, so do not read ``dt`` as the minimum of
   what precedes it.

   **When energy conservation is strongly violated.**  It should vary by about a percent
   or less over a run, and the trace in ``energy*.sph`` is where you see this.  If
   it does not, the usual response is to take smaller steps by tightening the
   Courant numbers, for instance ``cn1=0.2`` and ``cn3=0.1``, and then to
   restart from before the point where the trace went wrong rather than from
   the beginning: copy the ``out*.sph`` written just before that time to
   ``restartrad.sph`` and resubmit.

   The ``indx`` line is worth consulting first.  If one particle is repeatedly
   responsible, the problem is more likely to be that particle, perhaps one with
   a very small smoothing length in a dense region, than the Courant numbers, and
   tightening them will only slow the run down while it fails in the same way.

.. note::

   A run that stops with ``unreasonably small dt. abort run.`` has hit this from
   the other side: some particle drove the timestep to a value that would take
   forever to integrate.  The last ``indx`` line before the abort names it.

``energy*.sph``
~~~~~~~~~~~~~~~

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

Further columns appear conditionally, decided by three separate settings.
Count the columns rather than assuming.

.. dropdown:: Which extra columns appear, and when
   :icon: list-unordered

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
   ``nrelax >= 2`` selects, and not for a dynamical collision.  A collision run has
   ``nrelax=0`` and gets the ordinary seven columns even though there are plainly
   two stars in it.

Parsing an energy trace
^^^^^^^^^^^^^^^^^^^^^^^

Two habits will save you a wrong answer.

**Count the columns rather than assuming an index.**  The first seven are fixed, but
anything beyond them depends on ``ncooling``, ``reat`` and ``nrelax``, as above.
Code that reaches for column 8 will read the radiated energy in one run and the
orbital angular frequency in another, with nothing to announce the difference.

**Look for later stages.**  A resumed run writes ``energy1.sph``, then
``energy2.sph``, rather than appending to ``energy0.sph``.  Reading only
``energy0.sph`` after a restart silently gives you the first stage and stops,
which looks exactly like a run that ended early.

Reading a relaxation from this file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Column 5 should be nearly constant: it is the cleanest check that a run is
behaving.  Column 3 should fall to something negligible compared with column 4,
which is what a relaxed star looks like.

To measure an oscillation period, use column 2 rather than column 3.  The
gravitational potential energy oscillates once per cycle about a slowly varying
mean, and its zero crossings give a clean period.  The kinetic energy oscillates
at twice the mode frequency with a rapidly decaying envelope, and reading a
period from it is unreliable, because the answer depends on how the data were
detrended.

.. note::

   Judge a relaxation *after* the drag has been switched off at ``treloff``, not
   during it.  Under drag everything looks quiet, and two configurations that
   appear identical while relaxing can differ substantially once released.  Allow
   several oscillation periods after release before drawing a conclusion: an
   effect that has not appeared yet is indistinguishable from one that does not
   exist.

``out*.sph``
~~~~~~~~~~~~

The run settings followed by the full state of every particle, written every
``dtout``.  These are the snapshots you visualise, and the starting point for a
collision: the last snapshot of a relaxation is a relaxed star.

They are Fortran ``unformatted`` files, which is to say binary records rather
than text, so they are read with the tools in :doc:`visualisation` rather than a
text editor.

``restartrad.sph``
~~~~~~~~~~~~~~~~~~

The same format as ``out*.sph``, written every 1000 iterations and overwriting
the previous one, so a run that dies loses at most that much work.

If it is present in the directory when a run starts, the code continues from it
instead of starting afresh.  That makes it a deliberate control as well as a
safety net: copying any ``out*.sph`` onto ``restartrad.sph`` restarts the run
from that snapshot, which is how you back up past a point where something went
wrong.

.. warning::

   Because the file is picked up automatically, a leftover ``restartrad.sph``
   in a directory you meant to reuse for a fresh run will silently continue the
   old one.  If a new run begins at a time you did not expect, this is why.

Written for star relaxation runs
--------------------------------

``col*.sph`` and ``parent.sph``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two are written in the same units, and their first five columns hold the
same quantities in the same order, so a relaxation can be checked by plotting
one over the other.  ``parent.sph`` stops after column 6.

.. list-table::
   :header-rows: 1
   :widths: 12 44 44

   * - Column
     - ``col*.sph``
     - ``parent.sph``
   * - 1
     - radius
     - radius
   * - 2
     - pressure
     - pressure
   * - 3
     - density
     - density
   * - 4
     - temperature
     - temperature
   * - 5
     - mean molecular weight
     - mean molecular weight
   * - 6
     - particle mass
     - internal energy
   * - 7
     - smoothing length
     -
   * - 8
     - neighbour number
     -
   * - 9
     - gravitational acceleration, radial component
     -
   * - 10
     - hydrodynamic acceleration, radial component
     -
   * - 11-13
     - x, y, z
     -
   * - 14
     - gravitational potential
     -
   * - 15
     - internal energy variable
     -
   * - 16
     - velocity squared
     -

``col0000.sph`` against ``parent.sph`` shows how well the initial particle
layout matched the profile it was asked to reproduce.  A later one shows what
the relaxation did with it.  See :doc:`../tutorials/relaxing_a_star`.

.. note::

   Column 15 is the code's internal energy *variable*, whose meaning depends on
   ``nintvar``, and it is not the same quantity as column 6 of ``parent.sph``.
   Compare columns 1 to 5, which are unambiguous.
