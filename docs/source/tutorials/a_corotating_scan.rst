A corotating scan
=================

A scan walks a synchronized binary inward through a series of separations,
holding it in a corotating frame the whole way.  If the walk is slow enough that
the stars stay in equilibrium, every snapshot is an equilibrium configuration at
its own separation, and the run traces a quasi-equilibrium sequence rather than
an orbit.

That is useful when you want the structure of a binary as a function of
separation.  Where a star first fills its Roche lobe, and how its shape responds
as the companion closes in, both become readable without waiting for an orbit to
decay, and without letting the answer depend on how the system got there.

It is also an efficient way to build initial conditions.  One scan produces a
whole family of relaxed, synchronized binaries at known separations, and any of
them can start a dynamical run.  See `Using a snapshot as initial conditions`_.

A scan is not itself a dynamical calculation.  Nothing orbits, because the code
places the two centres of mass where the schedule says they should be, every
step.

What a scan run does
--------------------

Setting ``nrelax=2`` puts the calculation in a frame corotating with the binary
and adds the centrifugal acceleration.  The angular frequency is not something
you supply.  Each step ``getomega2`` solves for the :math:`\omega` that makes the
net force on the two centres of mass vanish::

    omega2 = -sum m (x ax + y ay) / sum m (x^2 + y^2)

so the frame always rotates at the rate the current configuration requires.  For
an extended star this is not the point-mass Keplerian value, and it should not
be, because the star's own moment of inertia enters the denominator.  On a giant
with a close companion the difference reaches tens of per cent.

Then ``cmadj`` repositions both centres of mass onto the scheduled separation
along the *x* axis and zeroes their centre-of-mass velocities.  The binary is
driven, not orbiting.

.. note::

   ``nrelax=3`` is the same frame with the Coriolis term added.  Use it when the
   gas is meant to move in the rotating frame.  For a scan, where the point is to
   keep everything as close to static as possible, ``nrelax=2`` is what you want.

Setting one up
--------------

Start from a relaxed star, as produced by :doc:`relaxing_a_star`.  Bring its last
snapshot in under the name the code expects::

    cp ../relaxation/out1400.sph sph.start1u

``sph.init`` selects the corotating setup::

    &INITT
    INAME='2cr'
    &END

The second body can be either another relaxed star, copied to ``sph.start2u``,
or a point mass.  **The absence of** ``sph.start2u`` **is what selects the point
mass**: when the code cannot find that file it builds a single particle of mass
``mbh``, in units of ``munit``, with smoothing length ``hco``.  There is no flag
to set.

A point mass is the right choice when the companion is small compared with the
smoothing lengths of the gas it will move through, since a resolved model of it
would buy nothing.  Check that before assuming it, by comparing the companion's
radius against the primary's ``hp`` in the region the companion will reach.

.. note::

   ``hco`` is fixed for the whole run.  Particles with :math:`u = 0` skip the
   smoothing-length update in ``advance.f90``, so unlike every other particle the
   point mass keeps the length it was given.  Choose it deliberately, because it
   is the softening of the companion's gravity for the rest of the calculation,
   and it is inherited by any dynamical run started from these snapshots.  Left
   unset, it defaults to the smoothing length of the primary's last particle,
   which is a surface particle and typically far larger than you want.

The schedule
------------

Two numbers set the geometry, ``sep0`` and ``sepfinal``, and one sets the timing,
``tscanon``.  Before ``tscanon`` the separation is held at ``sep0``, which gives
the star time to settle into the tidally distorted, synchronously rotating shape
the corotating frame demands.  It was relaxed in isolation, so it does not arrive
in that shape.  Several oscillation periods is a reasonable hold.

The scan itself is **exponential in separation, not linear**:

.. math::

   a(t) = a_0 \exp\!\left[-k\,(t - t_\mathrm{scanon})\right], \qquad
   k = \frac{\ln(a_0/a_\mathrm{final})}{t_\mathrm{scanoff} - t_\mathrm{scanon}}

so the separation changes by a constant *fraction* per unit time.  That is
usually what you want, since the tidal field scales with :math:`R/a`.  A linear
form is present in ``advance.f90`` but commented out.  Do not be misled by it
when reading the source.

.. important::

   ``tscanoff`` is not a setting.  The code uses

   .. math::

      t_\mathrm{scanoff} = \min(t_f,\ t_\mathrm{reloff})

   so the scan ends when the run ends or when the drag switches off, whichever
   comes first.  Since ``treloff`` is also the moment the binary is released into
   a dynamical calculation, **a single run cannot both finish a scan and then sit
   at the final separation**.  Choose one:

   * To trace the sequence and stop, set ``treloff`` beyond ``tf``.  The scan
     then spans ``tscanon`` to ``tf``, and the run never goes dynamical.
   * To trace the sequence and then release it, set ``treloff`` to the end of
     the scan and ``tf`` beyond it.

   For a settling phase at the final separation instead, see
   `Resuming freezes the separation`_.

.. note::

   For a corotating run that is not meant to scan at all, set ``sepfinal`` equal
   to ``sep0``.  Leaving ``sepfinal`` at its default stops the run during setup,
   which is the intended guard, and disabling the scan with a negative
   ``tscanon`` alone is not enough: the separation would still jump to
   ``sepfinal`` once ``t`` passed ``tscanoff``.

Choosing the rate
~~~~~~~~~~~~~~~~~

The scan has to be slow compared with the star's fundamental oscillation period,
not merely compared with the run.  A convenient way to state the rate is the
fractional change per period, which the exponential form makes constant at
:math:`1 - e^{-kP}` for a period :math:`P`.

A few per cent per oscillation period is a reasonable starting point, and the
cost of being wrong is one-sided.  Too fast leaves the star chasing an
equilibrium it never reaches, and the sequence is not quasi-static in any useful
sense.  Too slow only costs computer time.  If in doubt, run a short fast scan
first to confirm the mechanics, then commit.

The symptom of scanning too fast is the star ringing at the oscillation period
with an amplitude that grows as the separation closes, rather than settling into
a smoothly changing shape.  Compare snapshots at the same separation from two
scans at different rates.  If they differ, the faster one was too fast.

Choosing the drag
-----------------

The drag damps the oscillations that the changing tidal field excites, and that
is what keeps the sequence quasi-static.  A scan with the drag disabled will
ring.  A timescale near the star's oscillation period is a sensible default, as
it is for a single-star relaxation.

Setting ``trelax=0`` asks for that timescale to come from the model.  For a
corotating binary it is settled in ``initialize_corotating.f`` rather than in
``relax.f``, using star 1 alone, and it tries two things in order.

**The value the relaxation used**, if ``sph.start1u`` still carries one.  The
drag timescale travels in the snapshot header, so a star dumped while its
relaxation was still under drag brings its own answer with it.

**Otherwise star 1's own radius**, meaning the distance from that star's centre
to its farthest particle, together with that star's mass, giving
:math:`t_\mathrm{dyn} = \sqrt{R^3/M}` and
:math:`t_\mathrm{relax} = 4\,t_\mathrm{dyn}`, the same rule a single-star
relaxation uses.  The radius is measured before either star is moved out to
``sep0``, so it is the stellar radius and carries no trace of the separation.
Running the corotating example at ``sep0=8`` and at ``sep0=800`` returns the same
radius and the same ``trelax`` to every digit.

.. note::

   The second branch is the usual one.  A relaxation followed past ``treloff``
   switches its own drag off and records ``trelax=1.d30`` in every snapshot from
   then on, so the last dump of a properly released model carries a sentinel
   rather than a timescale.  Since the last dump is what you normally hand to
   ``startfile1``, expect the derivation rather than the inherited value, and
   read ``log0.sph`` to see which one was used.

On a resumed run neither applies, because the setup routine does not run again.
There the value settled the first time is read back from the dump header.

.. important::

   ``trelax=0`` does **not** set ``treloff`` here, though it does for a
   single-star relaxation.  That is deliberate.  For a corotating binary
   ``treloff`` also fixes ``tscanoff``, so deriving it would let the star
   silently redefine how long the scan lasts.  Set ``treloff`` yourself.

Reading the output
------------------

``energy0.sph``
   The last column is the current separation, so the file doubles as a record of
   the schedule actually followed.  Check it against what you asked for before
   trusting anything else.

``log0.sph``
   Records ``omega2`` each step, and ``binary separation ... at time=`` while the
   scan is running.  Watching ``omega2`` settle during the hold is the cleanest
   way to tell whether the star has finished adjusting to the corotating frame.

``out*.sph``
   The snapshots.  With the schedule above, a snapshot's separation follows from
   its time, so any snapshot can be placed on the sequence without opening it.

A scan stops itself when a particle strays past about 2.5 times either centre of
mass, reporting that it ``has overflowed an outer lagrangian point``.  One
particle is enough to end the run.  A giant's outermost particles sit far outside
the radius that holds most of its mass, so this fires well before the star as a
whole reaches its lobe.  The scan described below ended this way at a separation
of 646, with the radius enclosing 99 per cent of the mass at 0.98 of the Roche
lobe, under one per cent of the mass outside the lobe, and the single outermost
particle at 3.9 times it.  Read it as the limit on how far a scan can be pushed,
not as a measurement of contact.

Using a snapshot as initial conditions
--------------------------------------

Every ``out*.sph`` the scan writes is a relaxed, synchronized binary at a known
separation, and each one can start its own dynamical calculation.  A single scan
therefore replaces a whole set of separate binary relaxations, one per
separation, which is usually the cheapest way to get a series of encounters that
differ only in how close the pair began.

Work in a fresh directory, so that the scan's own output is not overwritten::

    mkdir ../plunge_from_0525
    cp out0525.sph ../plunge_from_0525/restartrad.sph

and in the new ``sph.input`` ask for a dynamical calculation::

    nrelax=0,
    treloff=0,
    trelax=1.d30,

.. important::

   ``treloff`` must be zero or negative, and this is the part that is easy to get
   wrong, because the scan's own ``sph.input`` carries a large one.  Copied over
   unchanged it breaks the release without stopping the run.  With ``nrelax=0``
   the boost into the inertial frame never fires, and the code integrates
   corotating-frame velocities as though they were inertial.  With ``nrelax=3``
   the calculation never leaves the driven state, so ``cmadj`` keeps forcing the
   separation and the outer Lagrangian abort stays armed.

   Note also that ``trelax`` is the e-folding time of the drag, so switching the
   drag off means making it huge rather than small.  A small ``trelax`` damps the
   star to a standstill in a few steps.

No ``sph.init`` is needed, because a run that finds ``restartrad.sph`` reads it
instead of building anything.  Leaving it out is worth doing deliberately, so
that a missing ``restartrad.sph`` stops the run rather than quietly starting a
fresh scan.

Two lines in ``log0.sph`` look alarming and are not.  ``did you mean to have
nrelax>=2?`` is printed whenever the dump was written by a corotating run and
the new one is not, which is exactly the case here, and the ``stop`` beneath it
is commented out.  ``sep0 & sepfinal set to`` some tiny number follows from the
same cause: with ``nrelax`` below 2 every particle is counted into star 1, so
``getcoms`` finds no second centre of mass and the separation it computes is
meaningless.  Neither value is ever used, because ``cmadj`` only runs when
``nrelax`` is 2 or more.  Both lines are absent when releasing with
``nrelax=3``.  The snapshot carries ``omega2`` in its header, so
the code has everything it needs to turn the corotating configuration into an
orbiting one.  It resets the clock and the snapshot counter, restores the
artificial viscosity coefficients, and gives every particle the orbital velocity
its position implies.  ``log0.sph`` reports this as ``assuming a dynamical
calculation is starting``, followed by ``main: going dynamic``.

.. important::

   The released orbit inherits the ``omega`` the corotating run had reached, not
   the one the separation implies.  Take a snapshot before the configuration has
   settled and the pair starts out rotating too slowly to hold itself apart, so
   the orbit is eccentric and the separation falls from the first step.  Started
   from a snapshot only 0.6 drag times into a relaxation, a test pair at a
   separation of 8 came out turning at 6.50 degrees per unit time against the
   7.30 a circular orbit wants, and closed to 7.31 within six time units.

   Check that ``omega2`` in ``log0.sph`` has flattened before picking a snapshot
   to build on.  During a scan it drifts as the separation changes, so what you
   are looking for is a smooth trend rather than the transient of the first few
   drag times.

.. note::

   Confirm the separation of the snapshot you picked before building anything on
   it, either from the last column of ``energy0.sph`` at that time or from the
   schedule.  Naming the new directory after the snapshot, as above, saves
   working it out twice.

.. important::

   ``tf`` has a floor you did not set.  On a restart ``init.f`` runs

   .. math::

      \mathrm{if}\ t \ge |t_f|:\quad t_f \leftarrow \mathrm{sign}
      \left(\max(|t_{f,\mathrm{old}}|,\ |t_f|),\ t_f\right)

   using the **dump's** ``t``, before ``t`` is reset to zero, and
   :math:`t_{f,\mathrm{old}}` is the ``tf`` the scan itself ran with.  So any
   ``tf`` below the dump's own time is silently replaced by the scan's.  A
   snapshot taken at :math:`t = 2.5\times10^{5}` out of a scan that ran to
   :math:`2.75\times10^{5}` cannot be given a shorter run than that, whatever
   ``tf`` says.  Check the ``Reset final time to TF=`` line in ``log0.sph``.

Which frame to release into is a choice, and the two are not the same physical
state.  ``nrelax=0`` overwrites every velocity with the rigid-rotation value its
position implies and zeroes the vertical component, so whatever residual motion
the corotating run still had is discarded.  ``nrelax=3`` keeps the velocities as
they stand and integrates in the rotating frame, but freezes ``omega2`` at the
dump's value, because ``getomega2`` is guarded by ``gonedynamic``.  That frame is
corotating only while the separation stays near where the snapshot was taken.

The residual is small in a settled, detached snapshot and large in one taken
during mass transfer, so the choice matters most exactly where the two runs would
be most interesting to compare.

Resuming freezes the separation
-------------------------------

.. important::

   A scanning run cannot be resumed and carry on scanning.  On restart,
   ``init.f`` sees ``nrelaxold >= 2`` and sets

   .. math::

      a_0 = a_\mathrm{final} = |x_{\mathrm{cm},2} - x_{\mathrm{cm},1}|

   freezing the separation wherever the run had reached.  It says so in
   ``log0.sph``, but the run continues rather than stopping, so an interrupted
   scan becomes a fixed-separation relaxation without ever failing.

   Plan a scan to complete in one uninterrupted block, and check the separation
   column of ``energy0.sph`` against the schedule before believing a resumed run.

This is deliberate, and it is the intended way to add a settling phase.  Let the
scan finish, then restart it, and it relaxes at the separation it ended on.  To
release *that* into a dynamical calculation, set ``treloff`` on the resumed run.

A worked example
----------------

A 5.29 :math:`M_\odot` AGB giant, relaxed from a MESA profile with 99697
particles, scanned against a 0.4 :math:`M_\odot` point mass from 1400 to 600
:math:`R_\odot`.  The giant's fundamental oscillation period is
:math:`P = 1.1\times10^{4}` in code units, which sets the whole schedule::

    &input
     tf=275000,
     dtout=500,
     tscanon=55000,
     sep0=1400,
     sepfinal=600,
     nrelax=2,
     trelax=10000,
     treloff=1.d30,
     startfile1='sph.start1u',
     mbh=0.4,
     hco=6.9,
     nnopt=75,
     nav=3, alpha=1, beta=2, ngr=3,
     neos=1, gam=1.560d0,
     runit=6.957d10,
     munit=1.9884098706980504d33,
    &end

That is a hold of 5 periods and a scan of 20, giving
:math:`k = 3.85\times10^{-6}` and a separation falling by 4.1 per cent per
oscillation period.  ``treloff`` sits above ``tf``, so the scan runs to the end
of the run and the binary is never released.

Everything else in the file is carried over unchanged from the relaxation that
produced ``sph.start1u``.  That matters more than it looks.  The star was relaxed
into equilibrium under a particular equation of state and neighbour number, and
changing either would leave it settling out of a shape it no longer fits.

``log0.sph`` confirms the setup on startup::

    init: new run, iname=2cr
    corotating: sep0=   1400.0000000000000
    particle           1 is a corepoint of mass  0.91193786057154913
    mass1=    5.2938290052005037
    n1=       99697
    black_hole_mass  0.40000000000000002
    mass2=   0.40000000000000002

Read the two masses and ``n1`` every time.  They are the cheapest confirmation
that the star you meant to scan is the one that was loaded, and that the point
mass came out at the mass you asked for.

Where to go next
----------------

A scan gives the sequence.  It does not tell you whether the system would
actually pass through it.  A binary whose primary carries more spin angular
momentum than a third of the orbital angular momentum has no stable synchronous
state, so a real system would run away inward rather than move quasi-statically
along the sequence.  The sequence is still the right object to compute, since it
is what the dynamical run should be compared against, but it is worth knowing
which part of it corresponds to a configuration a binary could occupy.

To continue into the plunge, take a snapshot as described above, and see
:doc:`a_collision` for what a dynamical two-body run writes.
