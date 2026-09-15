Jumping ahead on a wide orbit
=============================

A grazing encounter between a star and a black hole often leaves the star bound,
on an orbit so eccentric and so wide that the next pericentre passage is tens of
thousands of dynamical times away.  Integrating the intervening orbit is
pointless.  For all but the last little while of it the two bodies are a
two-body problem with a star quietly ringing inside it, and SPH is being asked
to do nothing that Kepler could not do exactly.

StarSmasher can skip that part.  The routine is ``jumpahead``, in
``parallel_bleeding_edge/src/skipahead.f``; the file and the routine have never
agreed on a name, and "skip ahead" and "jump ahead" both refer to it.  It
measures the orbit the two components are on, solves the two-body problem
analytically, and puts the system back down at a smaller separation on the
*infalling* branch of the same orbit, with every particle's position and
velocity relative to its own component untouched.  The hydrodynamics resumes
from there.

Where this has been used
------------------------

The technique goes back to the binary-disruption simulations of `Antonini,
Lombardi & Merritt (2011), ApJ 731, 128
<https://ui.adsabs.harvard.edu/abs/2011ApJ...731..128A/abstract>`_
(`arXiv:1008.5369 <https://arxiv.org/abs/1008.5369>`_), whose Section 3.3,
"Timescale considerations and orbital advancement", sets out the argument: the
thermal timescale of a bound star is :math:`10^5` to :math:`10^7` yr, far longer
than an orbital period, so the star's structure barely changes over an orbit and
nothing is lost by advancing it analytically.  They wait at least eight days
after periapsis before measuring the orbital elements, and check the
approximation against runs that integrate the orbit in full.

`Godet et al. (2014), ApJ 793, 105
<https://ui.adsabs.harvard.edu/abs/2014ApJ...793..105G/abstract>`_
(`arXiv:1408.1819 <https://arxiv.org/abs/1408.1819>`_) use the same treatment for
the repeated partial stripping of a donor by an intermediate-mass black hole in
HLX-1.  Their Section 6.1 states it compactly: once the donor "has retreated
sufficiently far from the black hole to become stabilized (typically about 100
dynamical timescales after periapsis), we employ the analytic Kepler two-body
result to advance the orbit to the same separation but now with the donor
infalling toward the BH", and "during this advancement of the orbit, we excise
from the simulation any particles that have been stripped from the star".
Sections 5 and 6 of the same paper discuss how sensitive the orbital evolution
is to the phase of the star's oscillation at pericentre, which is worth reading
before deciding what to do about artificial viscosity.

`Kıroğlu et al. (2023), ApJ 948, 89
<https://ui.adsabs.harvard.edu/abs/2023ApJ...948...89K/abstract>`_ call it
orbital regularization in their Section 2.2, and say most clearly why it is
needed: for :math:`M_{\rm BH} > 100\,M_\odot` at :math:`r_p = r_{\rm T}` the
remnant comes away with :math:`e > 0.999` and :math:`a > 10^4\,R_\odot`, an
orbital period of order ten years, or :math:`\sim 10^5` dynamical timescales.
They also record the choices that go with it: jump once the star has receded far
enough that the orbital elements are well determined and the star is close to
hydrostatic equilibrium; turn artificial viscosity off for the leg that follows,
so that oscillations in the remnant can be followed cleanly; and accept that all
debris bound to the black hole is treated as accreted, which is harmless when
:math:`M_{\rm BH} \gg M_*`.

Turning it on
-------------

``tjumpahead``
   The time at which the jump happens.  The default, ``1d30``, never fires.  Any
   other value is a deliberate request and is honoured.  The code stores it
   negated, and that negative sign is what tells ``changetf``, which revises the
   run's own schedule at every output, to leave the jump time alone.  ``main.f``
   compares against its absolute value, so you write it positive and never see
   the sign.

``throwaway``
   Whether the debris is discarded.  The default is ``.true.``, which is what
   the papers do.  See `Discarding the debris`_.

``tf``
   Not required for a jump, but set it negative anyway.  A negative ``tf`` lets
   the code revise its own stopping time, and makes it analyse the system at
   every output and write ``ecc.sph``, which is how you follow the orbit.

That is the whole interface.  Put ``tjumpahead`` in ``sph.input`` and run; the
jump fires on the first iteration past it.  It can go in from the start, or be
added later and the run resumed from ``restartrad.sph``, which is useful when
you would rather look at the first passage before committing to a jump time.

.. note::

   ``changetf`` can in principle schedule a jump on its own: when it sees a
   bound pair receding with an apocentre past 1000 code units it logs ``FUTURE
   CANDIDATE FOR JUMPING AHEAD``.  The line that would set a jump time there is
   commented out, so nothing follows from it, and the choice stays with you.

Choosing when to jump
---------------------

Jump too early and the orbital elements are still changing and the star is not
back in equilibrium; jump too late and you have paid for the integration you
were trying to avoid.

A separation that works
~~~~~~~~~~~~~~~~~~~~~~~

The prescription behind the published intermediate-mass black hole runs is a
separation,

.. math::

   r_{\rm jump} = 1.7 \left(\frac{M_{\rm BH}}{M_*}\right)^{1/3}
                  \max(r_p,\, r_{\rm T}),

which for :math:`r_p \le r_{\rm T}` is the same as
:math:`1.7\,(M_{\rm BH}/M_*)^{2/3} R_*`.  The coefficient comes from watching
what the debris does.  Measuring the star-black hole separation in
:math:`r_p = r_{\rm T}` runs with a 1 :math:`M_\odot` star, and writing
:math:`\sigma` for :math:`(M_{\rm BH}/M_*)^{1/3}\max(r_p,r_{\rm T})`:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Separation
     - What has happened by then
   * - :math:`1.4\,\sigma`
     - the first bound debris has fallen back to the black hole
   * - :math:`1.7\,\sigma`
     - the debris stream is just starting to self-intersect, and the disc radius
       has stopped growing
   * - :math:`2.0\,\sigma`
     - that material has been round once more

1.7 is the last moment before shocks from self-intersection begin to matter,
which is what makes it the right place to jump in a run with the artificial
viscosity turned off: there is nothing yet for the viscosity to do.  It is the
same criterion Kıroğlu et al. describe from the other end, noting that for
:math:`r_p \gtrsim r_{\rm T}` that separation is reached just as the stream
starts to self-intersect.

For :math:`M_{\rm BH}/M_* = 5, 10, 100, 200, 500` and 1000 the coefficient
:math:`1.7(M_{\rm BH}/M_*)^{1/3}` comes to 2.9, 3.7, 7.9, 9.9, 13.5 and 17.  At
the low end that is only a few tidal radii, close enough that the star may not
have finished being disrupted, so look at a snapshot before trusting the number.

Turning a separation into a time
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``tjumpahead`` is a time, so the prescription has to be converted.  The first
passage is near enough parabolic for Barker's equation, the parabolic
counterpart of Kepler's equation (`Pathan 2008, Math. Gaz. 92, 39
<https://www.cambridge.org/core/journals/mathematical-gazette/article/abs/eulers-and-barkers-equations-a-geometric-derivation-of-the-time-of-flight-along-parabolic-trajectories/3AEDFC36C19C75A01F91984247A603E4>`_,
or Section 4.5 of Bate, Mueller & White, *Fundamentals of Astrodynamics*).
Written in terms of :math:`x = r/r_p` it gives the time since pericentre as

.. math::

   t = \frac{\sqrt{2}}{6\pi}\,(x+2)\sqrt{x-1}\; t_{\rm orb}
     = 0.07503\,(x+2)\sqrt{x-1}\; t_{\rm orb},
   \qquad
   t_{\rm orb} = 2\pi\sqrt{\frac{r_p^3}{GM_{\rm BH}}}.

Putting :math:`x_{\rm jump} = 1.7(M_{\rm BH}/M_*)^{1/3}` into it gives
:math:`t_{\rm jump} = 1.9\,t_{\rm orb}` for :math:`M_{\rm BH} = 100\,M_\odot`
and :math:`5.7\,t_{\rm orb}` for :math:`1000\,M_\odot`, measured from
pericentre.  The same expression run with the starting separation gives the time
from the start of the run to pericentre, so the two together give a
``tjumpahead`` before the run has been started.

A worked example
----------------

The run below takes under an hour on one GPU and uses only files that come with
the repository.  It sends the relaxed 8 :math:`M_\odot` star from
``example_input/collision`` past a 100 :math:`M_\odot` black hole, which is what
``hyp`` produces when ``sph.start2u`` is absent: the second body becomes a single
point mass of mass ``mbh``.

Working out the jump time
~~~~~~~~~~~~~~~~~~~~~~~~~

The star is :math:`R_* = 3.17\,R_\odot`, so
:math:`r_{\rm T} = R_*(M_{\rm BH}/M_*)^{1/3} = 7.4`.  Taking :math:`r_p = 12`,
a grazing pass at about :math:`1.6\,r_{\rm T}`, the prescription gives
:math:`r_{\rm jump} = 1.7 \times 2.32 \times 12 = 47`.  With
:math:`t_{\rm orb} = 2\pi\sqrt{12^3/100} = 26` code units, Barker's equation puts
:math:`x = 47/12` at 20 code units after pericentre, and the starting separation
:math:`x = 60/12` at 27 before it.  So the jump wants to happen around
:math:`t = 47`; the orbit here is bound rather than exactly parabolic, which
brings pericentre in a little earlier, and 45 is a round number close enough.

Setting up
~~~~~~~~~~

Put the executable and the star in an empty directory::

    $ mkdir jump && cd jump
    $ cp ../parallel_bleeding_edge/parallel_bleeding_edge_gpu_sph .
    $ cp ../example_input/collision/sph.start1u .

``sph.init``::

     &INITT
     INAME='hyp' ! "hyperbolic" collision (also works for parabolic and eccentric encounters)
     &END

``sph.input``.  Note the negative ``tf``, and that ``sph.start2u`` is
deliberately absent (so that a black hole is used instead)::

     &input
     tf=-9999, ! a negative final time can help to make things more automatic
     dtout=10, ! time in code units between output files
     sep0=60, ! initial separation
     rp=12.0d0, ! periapsis distance
     e0=0.995d0, ! initial eccentricity
     mbh=100.0d0, ! black hole mass in solar masses
     tjumpahead=45.0, ! time in code units to do an orbital jump
     &end

``throwaway`` needs no line: its default is already ``.true.``.  Run it::

    $ mpirun -np 4 ./parallel_bleeding_edge_gpu_sph

The first passage
~~~~~~~~~~~~~~~~~

Reading the ``out*.sph`` snapshots back gives the orbit tightening as it passes
pericentre, which is the tidal energy going into the star:

.. code-block:: text

    out0000.sph  t=  0.06  r= 59.93  rdot=-1.69  a= 2515.6  e=0.99523
    out0001.sph  t= 10.02  r= 42.12  rdot=-1.90  a= 2528.9  e=0.99526
    out0002.sph  t= 20.02  r= 21.82  rdot=-2.11  a= 4166.8  e=0.99712
    out0003.sph  t= 30.01  r= 15.67  rdot=+1.79  a= 1269.8  e=0.99054
    out0004.sph  t= 40.00  r= 36.09  rdot=+1.98  a= 1199.7  e=0.99002

Pericentre is at :math:`t\approx24`, and by :math:`t=40` the semimajor axis has
dropped from 2516 to 1200.  ``ecc.sph`` carries the same story one row per
output.  Once the star is receding, ``log0.sph`` adds this at every output;
the numbers below are the ones from :math:`t=40`:

.. code-block:: text

    bound orbit with orbital period=   44277.381515463341
    the stars will take a long time to orbit, we might give up: ecc=  0.99313867634119102

The jump
~~~~~~~~

Everything that follows is from ``log0.sph``.  First the announcement, then
``compbest3`` sorting the particles, which takes three passes to converge:

.. code-block:: text

    jumpping ahead at time t=   45.000532609753130
    ...
              nit        nchng           m1           m2           m3
                0        19702  7.98693      100.012      0.00000
                1           88  7.98637      100.013      0.00000
                2            0  7.98637      100.013      0.00000

The star is left with :math:`7.98637\,M_\odot`, and the black hole's component
has gained :math:`0.0130\,M_\odot` of debris bound to it.  Almost none of the
stripped mass is unbound: the ``mejecta=`` reported at the previous output is
:math:`3.9\times10^{-4}`.

Then the orbit, measured and solved:

.. code-block:: text

    eccentricity: r12=   45.7046568     am1=   7.98637468     am2=   100.012968
    reduced mass mu=   7.3957953400447867
    total orbital energy= -0.2275
    total angular momentum=    376.1      -0.4639E-03   0.7433E-03    376.1
    components of eccentricity vector=  -0.9932      -0.3172E-02  -0.1219E-05
    apastron separation rmax=   3498.6677841672463
    sep0=   22.852328392088165
    1st simple check: -0.22751752517881471      -0.22751752517881466
    semilatusrectum=   23.940548138731852   mu=   7.3957953400447867
      sep0=   22.852328392088165   ecc=  0.99315723880756224
    cos=   4.7947739105261220E-002  sin= -0.99884984572992441
    rdot(from semilatusrectum): -2.1069863876044512  rdot(from e): -2.1069863876044521

The separation has been halved, 45.70 to 22.85.  :math:`\sin\theta` is negative
and so is :math:`\dot r`, which is what puts the star on the infalling branch,
and the two independent ways of computing :math:`\dot r` agree to sixteen
digits.  The orbital energy agrees with the value reconstructed from :math:`e`,
:math:`\mu` and :math:`L` to the same precision.

Then the three energy analyses, before the jump, after it, and after the debris
has been discarded.  These are also the three lines written to
``jumpahead.sph``, in the columns of ``energy*.sph``:

.. code-block:: text

                  t           epot           ekin           eint           etot         ajtot
       44.99822      -43.92045      17.29951       13.53356      -13.08738      376.6425
       44.99822      -61.39812      34.77565       13.53356      -13.08891      376.6425
       44.99822      -61.29525      34.73058       13.53313      -13.03154      376.2588

Line 1 to line 2 is the jump.  The internal energy is identical to every digit
printed, because nothing was done to the star, and so is the total angular
momentum.  The potential and kinetic energies both change, and should: the star
has been moved from :math:`r=45.7` to :math:`r=22.9`, so it is deeper in the
black hole's potential and moving faster.  The total energy moves in the fifth
digit.

Line 2 to line 3 is the mass removal:

.. code-block:: text

    now we throw away some particles:
    mass thrownaway=   1.3416840884748770E-002
    new ntot=       17729

``etot`` falls from :math:`-13.0889` to :math:`-13.0315` and ``ajtot`` from
376.64 to 376.26, which is what the deleted particles carried off with them.

``m1m2rp.sph`` records the result:

.. code-block:: text

       17728.000000000000        1.0000000000000000        12.000000000000000
       7.9863746823827269        100.01296764968698        12.000000000000000

17728 SPH particles in the star, one point mass, and ``rp=12``; then the
component masses.  The second component weighs :math:`100.0130`, the black hole
plus the debris bound to it, and that is the mass the Kepler solve used.  The
point particle itself is still exactly :math:`100`, so the leg that follows is
integrated with the original ``mbh``.

What it bought
~~~~~~~~~~~~~~

With :math:`a = 1755` and :math:`e = 0.99316` the orbital period is 44464 code
units.  Kepler's equation gives the time of flight from :math:`r=45.7` outbound,
around apocentre at 3499, back to :math:`r=22.9` infalling: 44439 code units, or
99.94% of a full period.  The whole run up to the jump spans 45.

The snapshots after the jump show the star going straight back in and round
again:

.. code-block:: text

    out0005.sph  t= 50.00  ntot=17729  r= 13.44  rdot=-1.30  a=1748.2  e=0.99313
    out0006.sph  t= 60.00  ntot=17729  r= 24.75  rdot=+2.10  a= 881.6  e=0.98641

The second pericentre passage is over by :math:`t\approx53`, and the semimajor
axis has halved again.  Without the jump it would have arrived somewhere near
:math:`t = 45000`.

How it works
------------

The two-body solve
~~~~~~~~~~~~~~~~~~

From the two component masses and their centre-of-mass positions and velocities
the routine forms the reduced mass :math:`\mu`, the orbital energy
:math:`E_{\rm orb}`, the angular momentum :math:`\mathbf{L}`, and the
eccentricity vector

.. math::

   \mathbf{e} = \frac{\mathbf{v}\times\mathbf{L}}{k} - \hat{\mathbf{r}},
   \qquad k = G m_1 m_2,

the Laplace-Runge-Lenz vector scaled to have magnitude :math:`e`.  It is
conserved in the Kepler problem and points at pericentre, so it fixes the
orientation of the orbit as well as its shape.  As a check, :math:`e^2` computed
from the energy and angular momentum is compared against
:math:`|\mathbf{e}|^2`, and the run stops if they disagree by more than
:math:`10^{-14}`.

The new separation is then chosen, and the true anomaly that goes with it
follows from the orbit equation,

.. math::

   \cos\theta = \frac{1}{e}\left(\frac{\alpha}{r} - 1\right),
   \qquad \alpha = \frac{L^2}{\mu k},

with :math:`\theta` taken *negative*.  :math:`\theta = 0` is pericentre, so a
negative :math:`\theta` is the pre-pericentre branch: the star is placed where it
is falling in, not where it is climbing out.  The radial velocity :math:`\dot r`
is taken negative for the same reason.

.. note::

   The papers advance the orbit to the *same* separation the star had reached,
   and so did the version of the code the published runs were made with.  The
   code as it stands halves it: ``sep0=0.5d0*r12``, with the unhalved
   ``sep0=r12`` commented out on the line below.  Halving skips more of the
   orbit, at the cost of resuming the hydrodynamics closer in, and it breaks the
   property the published runs relied on, that the separation after a jump
   equals the separation before it, so each passage takes the same time to reach
   the next pericentre.  That is the line to change to get it back.

Each particle is then translated and boosted by the difference between its
component's new and old centre-of-mass state.  Every particle keeps its position
and velocity *relative to its own component*, so the star's internal structure,
its spin, and whatever oscillation the encounter left it ringing with all carry
across untouched.  That is the approximation the whole method rests on, and it
is what Antonini et al. justify with the thermal timescale argument.

.. warning::

   **The simulation clock does not advance.**  ``jumpahead`` changes ``t`` by
   about one timestep and no more.  The orbital time that was skipped, which is
   most of a period, is simply not counted.  Every time in ``log*.sph``,
   ``energy*.sph``, ``ecc.sph`` and the ``out*.sph`` headers after a jump is a
   time with the wide part of the orbit cut out of it.  If you need real elapsed
   time, add the Kepler time of flight yourself, from the :math:`a` and
   :math:`e` printed in the log.

.. dropdown:: Synchronising the leapfrog first

   The integrator is a leapfrog, so velocities and internal energies are half a
   step ahead of the positions.  ``jumpahead`` rolls them back by
   :math:`\mathrm{d}t/2` so that everything refers to one instant, recomputes
   densities, smoothing lengths and gravity, and writes the energy summary to
   the log between ``***analyze system right before jump ahead:***`` and
   ``***done analyzing system right before jump ahead***``.  The matching pair
   after the jump is what you compare it against.  At the end, ``lfstart``
   restarts the leapfrog.

.. dropdown:: Splitting the system into two components

   The Kepler problem needs two bodies, so the particles have to be divided.
   ``compbest3`` does it, the same component finder that fills ``ecc.sph``: it
   iterates a binding test until membership stops changing, and puts material
   bound to neither body into a fourth, unbound component.  If a third component
   turns out to be more massive than the lighter of the first two, it is renamed
   into its place, so the jump always describes the two most important bodies.

   This happens whichever way ``throwaway`` is set.  It has to: the two-body
   solve needs the centre of mass of each *body*, and sorting particles by
   whether they are point masses is not good enough once anything has been
   stripped.  A tidal tail drags the centre of mass of "every SPH particle" well
   away from the centre of mass of the star, and it would also mistake a red
   giant's core particle for the compact object.

.. dropdown:: Rotating the new orbit into place

   The positions and velocities are worked out in the orbital plane, so they
   have to be rotated into the simulation frame.  Three Euler angles do it:
   :math:`\theta_1` from :math:`\hat{\mathbf{L}}\cdot\hat{\mathbf{z}}`,
   :math:`\phi_1` from the line of nodes, and :math:`\psi_1` chosen so that
   applying the rotation to :math:`(-e,0,0)` reproduces the eccentricity vector
   measured before the jump.  The new orbit therefore lies in the same plane as
   the old one and has pericentre in the same direction.  Nothing about the
   orbit changes except where along it the star sits.

.. dropdown:: The five checks

   Before any particle is moved, the reconstructed separation, orbital energy,
   angular momentum, eccentricity and semi-latus rectum must each match the
   measured value to one part in :math:`10^{8}`, and the run stops if any does
   not.  A jump that runs to completion has already proved that it conserved the
   orbit.

.. _Discarding the debris:

Discarding the debris
~~~~~~~~~~~~~~~~~~~~~

With ``throwaway=.true.``, the default, every particle that is neither a point
mass nor a member of the surviving body's component is deleted, the particle
count is compacted, and the mass removed is reported as ``mass thrownaway=``.
That covers both the material thrown to infinity and the material left bound to
the black hole.  If fewer than ``nnopt`` particles survive, the run stops rather
than continue with a star it cannot resolve.  When there is no point mass
anywhere, as in a star-star encounter, there is no accretor either, so both
bodies are kept and only material bound to neither goes.

With ``throwaway=.false.`` nothing is deleted, and the particles in neither body
are each moved with whichever body they are more tightly bound to.  The count is
logged as ``throwaway is off, so N particles belonging to neither body are
carried along``.  That is an approximation, since such a particle is not on
either body's orbit, but a better one than leaving the debris behind while both
bodies are moved out from under it.

The Kepler solve uses the component masses from ``compbest3``, so debris bound
to the black hole is counted in the black hole's mass *for the purpose of the
orbit*.  The point particle's own mass is not increased, so the leg after the
jump is integrated with the original ``mbh``.  That is the approximation Kıroğlu
et al. describe and defend on the grounds that :math:`M_{\rm BH}\gg M_*`; if
your mass ratio is less extreme, this is the place to look first.

How much of the bound debris ought to count as accreted is a real question, and
the answer is "not all of it".  `Ayal, Livio & Piran (2000), ApJ 545, 772
<https://ui.adsabs.harvard.edu/abs/2000ApJ...545..772A/abstract>`_ find that
around 75% of the returned debris becomes unbound again for supermassive black
holes, and measurements on :math:`M_{\rm BH} = 500\,M_\odot` runs here put the
fraction that stays bound at 10% or less.  Since the black hole is far more
massive than anything being argued about, the hydrodynamics is essentially
unaffected by the choice, and accretion rates can be rescaled afterwards by
whatever factor a later reader prefers.

.. warning::

   Discarding the debris is not a neutral bookkeeping choice for the *orbit*.
   Test runs that kept the material bound to the black hole ended up on markedly
   more tightly bound orbits, and in cases where discarding it led to the star
   being ejected after a few passages, keeping it postponed the ejection or
   prevented it.  The question to ask is how the orbital period compares with
   the timescale on which the disc would really be cleared, by accretion and by
   feedback winds, neither of which is in the simulation.  For intermediate-mass
   black holes the orbital period is long enough that discarding is the better
   approximation.  For stellar-mass black holes it is much less clear, and the
   two settings are worth running against each other.

Checking that a jump was sane
-----------------------------

The routine stops the run itself if the reconstructed orbit does not match the
measured one, so the arithmetic needs no checking.  What is worth checking is
everything it does not guarantee:

* The three lines in ``jumpahead.sph``.  ``eint`` and ``ajtot`` should be
  unchanged across the jump to every digit printed.  ``epot`` and ``ekin`` will
  both change, since the separation changed.  ``etot`` should barely move; if it
  moves by a percent, the two components are not well described as point
  masses, which usually means the jump was made too soon after pericentre.
* The component masses in ``m1m2rp.sph``.  The excess of the accretor's
  component over ``mbh`` is the debris treated as accreted, and is worth
  knowing: it is the mass the orbit was solved with but the point particle was
  not given.
* That the star was really relaxed when you jumped.  Nothing in the code checks
  this, and it is the assumption that matters most.
* That enough particles are left.  Each jump with ``throwaway=.true.`` removes
  some, and several passages in a row can leave a star without enough particles
  to find neighbours; an average neighbour count in single figures is not doing
  SPH any more, whatever the run says.  The code stops with ``not enough
  particles to continue`` only when fewer than ``nnopt`` survive, which is later
  than you want to find out.  ``ls -l out*.sph`` is the quick check: with
  ``throwaway`` on the snapshot size drops at every jump, so the file sizes show
  both where the jumps were and how much was lost.

There is also a self-test built into the routine.  ``skipahead.f`` carries two
commented-out lines, ``sintheta=-sintheta`` and ``rdot=-rdot``, flagged in the
source as a pair.  Uncomment both and the Kepler solve puts the system back
exactly where it already was, so the jump becomes the identity while everything
around it still runs: the component split, the rotation, the particle move, the
five checks.  Anything that moves is then a bug.  It is the first thing to reach
for after touching the routine.

The same component finder runs in the post-processing.  ``compbest3.f`` is
shipped with the splot routines in ``splot_routines/`` as well, where option 71
calls it to build ``massAndMore.out``.  Bound masses quoted from a jump and
bound masses measured afterwards therefore come from the same algorithm.

.. seealso::

   :doc:`../reference/sph_input` for ``tjumpahead``, ``throwaway`` and ``tf``,
   and :doc:`output` for the files named here.
