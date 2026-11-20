Output files
============

Every run writes these:

``log*.sph``
   How the run was launched, then a summary of the system at each timestep.  The
   first place to look for clues when something goes wrong.

``energy*.sph``
   Energies and related quantities, one row per output (with block timesteps,
   one per full synchronization, which is usually many per output), with the
   time in the first column.

``out*.sph``
   The run settings and the full particle state, written every ``dtout``.
   Unformatted, so not readable in a text editor.  These are what you visualise,
   and what a collision starts from.

``restartrad.sph``
   The same format as ``out*.sph``, rewritten every 1000 iterations as a
   checkpoint (with block timesteps, at the first moment all particles are in
   step after 1000 substeps).  If it is present in the directory when a run
   starts, the code continues from it rather than beginning afresh.

A run with a negative ``tf`` or with automatic jumping on writes one more, and a run
that jumps ahead on its orbit writes two others.  The first two are numbered
like ``energy*.sph`` and ``log*.sph``, so a restart starts new ones and keeps
the old:

``ecc*.sph``
   One row per output, tracking the two most massive components of the system:
   their masses and positions, their separation, the mass that has left them,
   and the semimajor axis and eccentricity of the orbit between them.

``jumpahead*.sph``
   The same columns as ``energy*.sph``, three lines per orbital jump, from
   either side of it.

``m1m2rp.sph``
   Written by ``jumpahead``.  Two lines: the particle counts and then the
   masses of the two components, each followed by ``rp``.

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

   Every step (every full synchronization with block timesteps, see below)
   writes a pair of lines::

      dts= <dt1> <dt2> <dt3> <dt4> <dt5> <dt6> <dt>
      indx <i1>  <i2>  <i3>  <i4>  <i5>  <i6>

   ``dts=`` gives the smallest timestep each criterion would have allowed, and
   then ``dt``, the step actually taken.  ``indx`` gives the particle
   responsible for each of those minima, which is what turns "the timestep
   collapsed" into "the timestep collapsed because of *this* particle".

   With block timesteps (``nblock=1``) the timesteps are evaluated at every
   substep, which can mean tens of thousands of times between full
   synchronizations.  Rather than a pair of lines for each, the pair is written
   once per full synchronization, just before the ``block sync`` line.  Each
   column then holds the smallest value that criterion reached over all the
   substeps since the previous synchronization, ``indx`` names the particle
   that set it, and ``dt`` is the smallest step any particle wanted in that
   interval.  The ``block sync`` line and bin histogram that follow say how
   many substeps there were and how the particles were spread across steps.

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
   off unless you turn them on.  When both are ``1d29`` or more the code does
   not compute the smoothed acceleration :math:`a_{sm}` at all (it is a sum over
   the neighbours of every particle, a large part of the step of a dynamical
   run), and the ``dt2`` and ``dt4`` columns are then made with :math:`a` itself,
   that is :math:`a_{sm}=0`: they are still enormous, and the step, the dumps and
   ``energy*.sph`` are not changed.  With either of them below ``1d29`` the
   smoothed acceleration is computed as before.  And ``dt`` is smaller than every individual
   column, because reciprocals add, so do not read ``dt`` as the minimum of
   what precedes it.

   **The Courant numbers and the energy.**  In a run with shared steps (``nblock=0``)
   ``cn1`` and ``cn3`` set the step, because ``cn2`` and ``cn4`` are off by
   default: ``dt1`` is ``cn1`` times the time a signal takes to cross a
   smoothing length, and ``dt3`` is ``cn3`` times the time in which a
   particle's internal energy would change by its own value at its present
   rate.  The defaults are ``cn1=0.3`` and ``cn3=0.1``.  Larger numbers give
   longer and fewer steps and a larger change of the total energy; what follows
   is how much, as measured.

   The first measure is the collision of two Suns of 99954 particles each, on a
   parabolic orbit with ``sep0=15`` and ``rp=1`` (pericentre at t=20.7), with
   shared steps and the gravity tree on the GPU (``ngravtree=2``,
   ``treeacc=0.45``, ``nleafgrav=24``), run to ``tf=30``.  Here
   :math:`\Delta E/E = (E(t)-E(0))/|E(0)|`, with :math:`E` the fifth column
   of ``energy*.sph``:

   .. list-table::
      :header-rows: 1
      :widths: 20 20 30 30

      * - ``cn1`` / ``cn3``
        - steps
        - :math:`\Delta E/E` at ``tf=30``
        - first t with :math:`|\Delta E/E| > 10^{-3}`
      * - 0.3 / 0.1 (default)
        - 6837
        - 5.4e-5
        - never
      * - 0.45 / 0.15
        - 4513
        - 2.7e-4
        - never
      * - 0.45 / 0.2
        - 3871
        - 3.9e-4
        - never
      * - 0.5 / 0.1667
        - 4123
        - 5.7e-4
        - never
      * - 0.55 / 0.1833
        - 3773
        - 1.3e-3
        - 18.6
      * - 0.6 / 0.2
        - 3790
        - 3.4e-3
        - 0.2
      * - 0.6 / 0.06
        - 8596
        - 1.7e-3
        - 0.6

   With 0.45 and 0.2 the steps of this run take 32.4 s instead of 54.7 s on one
   H100 with 63 threads (five runs of each, alternated), and 70.2 s instead of
   121.6 s on one L40S with 8 threads (three runs of each).

   This page counts :math:`|\Delta E/E| \le 10^{-3}` as good and
   :math:`|\Delta E/E| \le 5 \times 10^{-3}` as the limit.  Published StarSmasher
   runs report the total energy conserved to 1.5e-5 at the end of a run
   (Ballone et al. 2023, arXiv:2204.03493), to 0.004% in most runs and always
   within 0.5% (Kremer et al. 2022, arXiv:2201.12368), and to about 0.2% over
   the interval shown (Gaburov, Lombardi and Portegies Zwart 2010,
   arXiv:0904.0997).  ``cn1`` should stay at 0.5 or below: with ``cn1=0.55``
   the standard collision leaves the good range, with ``cn1=0.6`` the energy
   jumps in the first steps, before the stars touch (1.5e-3 at t=0.2 for
   ``cn3=0.2``), and with 0.75 and 0.25 the run stops at t=21.8 with
   ``unreasonably small dt``.

   Other collisions of two stars, the same code and shared steps, to the ``tf``
   in the first column (the stars of 19938 and 50160 particles are relaxed
   Suns made with the same recipe as the stars of 99954; the hyperbolic orbit
   has ``rp=1`` and ``vinf2=1`` in code units):

   .. list-table::
      :header-rows: 1
      :widths: 40 20 20 20

      * - collision (``tf``)
        - 0.3 / 0.1
        - 0.45 / 0.2
        - 0.6 / 0.06
      * - standard, 99954 particles each (30)
        - 6837 steps, 5.4e-5
        - 3871 steps, 3.9e-4
        - 8596 steps, 1.7e-3
      * - 50160 particles each (30)
        - 5661 steps, 4.4e-5
        - 3174 steps, 7.3e-4
        - 7388 steps, 1.5e-2
      * - 19938 particles each (30)
        - 3863 steps, 2.5e-5
        - 2160 steps, 1.3e-3
        - 4914 steps, 1.3e-2
      * - hyperbolic, ``vinf2=1`` (40)
        - 8261 steps, 5.1e-5
        - 4760 steps, 1.6e-3
        - 11137 steps, 1.6e-3
      * - grazing, ``rp=1.6`` (60)
        - 9017 steps, 2.1e-5
        - 5414 steps, 6.7e-4
        - 9621 steps, 3.3e-3
      * - nearly head on, ``rp=0.25`` (45)
        - 12430 steps, 1.7e-3
        - stopped at t=19.7056, ``unreasonably small dt``
        - 15368 steps, 9.9e-3

   The pair 0.45 and 0.2 stays below :math:`10^{-3}` in three of the six
   collisions and between :math:`10^{-3}` and :math:`5 \times 10^{-3}` in two
   (the hyperbolic one and the one of 19938 particles); in the nearly head-on
   collision it stops the run: the internal energy diverges at the largest
   compression (t=19.70), with the SPH on the GPU and on the CPUs alike, while
   the runs started from the same checkpoint with 0.3 and 0.1, 0.45 and 0.1 or
   0.3 and 0.2 go through it.  In that collision the defaults give 1.7e-3, the first time above :math:`10^{-3}` is t=27.9, and none of the
   pairs tried keeps it below :math:`10^{-3}`: with ``cn1`` and ``cn3`` equal to
   0.3 and 0.2 it ends at 4.9e-3, with 0.45 and 0.1 at 4.3e-3, with 0.4 and 0.2
   at 7.3e-3, with 0.45 and 0.15 at 8.7e-3.  In the hyperbolic collision and in
   the one of 19938 particles 0.45 and 0.15 end at 1.2e-3 and 1.1e-3, and 0.4
   and 0.2 at 6.8e-4 and 5.3e-4.

   In every run with 0.45 and 0.2 that finished, :math:`\Delta E/E` was still
   growing at ``tf``, by 7 to 25% over the last tenth of the time, so the
   numbers belong to that ``tf``.

   With the direct sum for the gravity (``ngravtree=0``) the standard collision
   takes 6851 steps and ends at 4.6e-5 with 0.3 and 0.1, 3890 steps and 4.0e-4
   with 0.45 and 0.2, and 8660 steps and 2.5e-3 with 0.6 and 0.06.  With
   ``ngpusph=3`` instead of 7 it takes the same steps and gives the same dump.

   The relaxed polytrope of the test ``energy`` (10000 particles, ``tf=10``)
   drifts by 2.8e-5 with 0.3 and 0.1, by 2.5e-4 with 0.4 and 0.2 and by 1.4e-3
   with 0.45 and 0.2.

   The defaults therefore stay ``cn1=0.3`` and ``cn3=0.1``, the only pair of
   the table below :math:`10^{-3}` in every case but the nearly head-on
   collision.  For a parabolic or grazing collision of well resolved stars
   ``cn1=0.45d0`` and ``cn3=0.2d0`` in ``sph.input`` give the run in 0.56 to
   0.60 of the steps within :math:`10^{-3}`; it is not a pair for a nearly
   head-on or a hyperbolic encounter, nor for stars of a few tens of thousands
   of particles.  These numbers are for collisions of two stars with shared
   steps: for another problem read the energy column of ``energy*.sph`` and
   take smaller numbers if it moves by more than you accept.  The files in
   ``example_input`` write their own values (``cn1=.6d0`` or ``.4d0`` and
   ``cn3=0.06d0``) and do not use the defaults.

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

One row per output, written by ``enout`` in ``output.f``.  With block timesteps
(``nblock=1``) a row is written at every full synchronization instead, every
moment at which all the particles' steps end together.  That is whenever the
longest step any particle is taking ends, so the rows come as often as the
particles allow: in a partial disruption of a red giant, with ``dtout=10``,
every 0.08 to 0.6 time units, rather than every 10.  The spacing changes as the
longest step does, so the rows are not evenly spaced in time.  When some
particles take the longest step allowed, ``dtmaxblk`` (``dtout`` by default),
the rows come once per ``dtmaxblk``; set ``dtmaxblk`` smaller for more.  The summary in
``log*.sph``, the ``dts=`` lines and ``restartrad.sph`` still come only at the
end of each ``dtmaxblk``.  The first seven
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
``dtout``.  The number in the name counts outputs, with four digits up to
``out9999.sph`` and as many as needed (up to seven) after that.  These are the snapshots you visualise, and the starting point for a
collision: the last snapshot of a relaxation is a relaxed star.

They are Fortran ``unformatted`` files, which is to say binary records rather
than text, so they are read with the tools in :doc:`visualisation` rather than a
text editor.

With ``nhdf5=1`` in ``sph.input`` (and a code built with ``make HDF5=1``) each
snapshot is also written as ``out*.h5``, in HDF5, which Python reads directly;
see :doc:`hdf5`.

``restartrad.sph``
~~~~~~~~~~~~~~~~~~

The same format as ``out*.sph``, written every 1000 iterations and overwriting
the previous one, so a run that dies loses at most that much work.  With block
timesteps (``nblock=1``) it is written instead at the first moment all particles
are in step (every ``dtmaxblk``) at least 1000 substeps after the previous one,
so that it holds exactly what an ``out*.sph`` written then would.

If it is present in the directory when a run starts, the code continues from it
instead of starting afresh.  That makes it a deliberate control as well as a
safety net: copying any ``out*.sph`` onto ``restartrad.sph`` restarts the run
from that snapshot, which is how you back up past a point where something went
wrong.

.. warning::

   Because the file is picked up automatically, a leftover ``restartrad.sph``
   in a directory you meant to reuse for a fresh run will silently continue the
   old one.  If a new run begins at a time you did not expect, this is why.

Written when ``tf`` is negative
-------------------------------

``ecc*.sph``
~~~~~~~~~~~~

A negative ``tf`` lets the code revise its own stopping time, and to do that it
analyses the system at every output.  ``ecc*.sph`` is the record of that
analysis: plain text, one row per output, twenty-six columns.  A run with
automatic jumping on writes it too.

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - Columns
     - Meaning
   * - 1
     - time
   * - 2-3
     - masses of the two most massive components, heaviest first
   * - 4-9
     - their centre-of-mass positions, three each
   * - 10
     - the separation between them
   * - 11-12
     - all the mass that has left both stars, as a fraction of the total and in
       code units
   * - 13-14
     - semimajor axis and eccentricity of the orbit between the two
   * - 15
     - the part of that mass which is unbound from the system as well
   * - 16-21
     - their centre-of-mass velocities, three each
   * - 22-23
     - the radius enclosing ``rjumppercent`` percent of each one's mass, binned
       as ``rstar95`` and ``rstar99`` in ``massAndMore.out``, and labelled by the
       percentage (``rstar95_1``, ``rstar95_2`` for the default 95).  Measured
       only when ``rjumpfactor`` sets the jump separation, and 0 otherwise
   * - 24
     - the larger of the two tidal radii :math:`R(M_{\rm other}/M)^{1/3}`,
       0 when not measured
   * - 25
     - ``rjump``, the jump separation in use (``1d30`` when jumping is off)
   * - 26
     - the period of their orbit when it is bound, and 0 otherwise

The first of the two components is the more massive one, which in a star and
black hole encounter means the black hole *together with whatever is bound to
it*, not the point particle alone.  The second is the star.

Columns 12 and 15 are easy to confuse and are not the same quantity.  Column 12
is everything that has come off the stars, whether or not it is still bound to
the system; column 15 counts only the part that is leaving for good.  The
difference between them is material that is unbound from either star but still
bound to the pair, a common envelope in the case this bookkeeping was written
for.  Column 12 is usually the one wanted.

The components are found by ``compbest3``, which iterates a binding test until
membership stops changing.  Columns 13 and 14 are what tell you whether an
encounter has left a bound remnant, and are what you watch when deciding
whether to :doc:`jump ahead <jumping_ahead>`.

``jumpahead*.sph``
~~~~~~~~~~~~~~~~~~

Written only when a jump can happen, and appended to at each one.  The columns
are exactly those of ``energy*.sph``, and the lines come in groups of three:

#. the system as it was before anything was touched,
#. after the orbit has been advanced, but before any mass is discarded,
#. after the mass has been discarded.

Comparing the first two lines isolates what the Kepler solve did; comparing the
second and third isolates what the mass removal did.  Nothing else writes to
this file, so it is three lines per jump and nothing else.

Like ``ecc*.sph``, it carries the number of the ``energy*.sph`` and
``log*.sph`` written alongside it, so each restart starts a new one and the
record of earlier passages is kept.

.. note::

   Fortran buffers this file, so a run in progress may show an ``ecc*.sph`` that
   is empty or many rows behind, and a run that was killed rather than allowed
   to finish may leave it empty altogether.  ``jumpahead*.sph`` below behaves the
   same way.  Read the ``out*.sph`` snapshots instead if you need the current
   state of a live run.

``m1m2rp.sph``
~~~~~~~~~~~~~~

Written by ``jumpahead``, and overwritten by each jump.  Two lines, three
numbers each: the first gives the number of particles in each component and
``rp``, the second the mass of each component and ``rp`` again.  See
:doc:`jumping_ahead`.

Written for star relaxation runs
--------------------------------

``col*.sph`` and ``parent.sph``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two are written in the same units, and their first five columns hold the
same quantities in the same order, so a relaxation can be checked by plotting
one over the other.  ``parent.sph`` stops after column 6.

Everything in both files is in code units, with two exceptions: temperature is
in kelvin and mean molecular weight is in grams.  What a code unit is worth in
cgs, and how to change it, is set out under :ref:`code units <code-units>`.

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
     - specific internal energy
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
