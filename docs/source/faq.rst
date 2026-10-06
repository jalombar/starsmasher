Frequently asked questions
==========================

**The build stops with** ``cannot find -lcudart``.

Your CUDA installation is not where the build is looking.  Run ``make config``
to see what it found, then ``make CUDAPATH=/path/to/cuda``.

**The build stops with** ``Unsupported gpu architecture``.

Your CUDA release no longer supports the architecture being requested.  Recent
releases have dropped older ones.  CUDA 13, for instance, no longer accepts
``sm_61``.  The gravity library normally detects your card from ``nvidia-smi``,
and if that is unavailable or not working you can pass the compute capability
yourself with
``make COMPUTE_CAPABILITY=75`` (where "75" is replaced with ten times the compute capability).

**I have no NVIDIA card.  Can I still use StarSmasher?**

Yes.  ``make cpu`` builds a version that needs no GPU.  It is not as slow as you
might expect, because the GPU only takes over the gravity.  See
:doc:`using/running` for measured timings.

**The run stops immediately with** ``init: error reading input file sph.init``.

A run needs ``sph.init``.  It is three lines.  See :doc:`reference/sph_init`.

**The run stops with** ``Cannot match namelist object name``.

``sph.input`` sets a variable that is not in the namelist.  Check it against
:doc:`reference/sph_input`.

**My relaxed star has particles flying off.**

Often the particle-mass distribution.  For a giant star spanning many orders of
magnitude in density, ``equalmass=0`` doesn't resolve the central regions well, which launches a shock wave outward and the low-mass outer particles never
settle.  See :doc:`tutorials/relaxing_a_star`.

**How many neighbours does** ``nnopt`` **give me?**

More than ``nnopt``: the actual count is larger by roughly 1.4 to 1.7 depending
on the kernel.  ``nnopt`` is the target of the smoothing-length constraint, not a
neighbour count.  See :doc:`reference/equations_of_motion`.

**Which of** ``parallel_bleeding_edge`` **and** ``Blackollider`` **should I use?**

``parallel_bleeding_edge`` unless you are modelling collisions with significantly massive compact
objects treated as point masses.  The two differ mainly in how the smoothing
length is related to the density.

**My restarted run isn't appending to energy0.sph.**

It's not supposed to.  Output is numbered per stage: a
resumed run writes ``log1.sph`` and ``energy1.sph``, for example, rather than appending to the
originals.

**With block timesteps, what does "mean refreshed" in** ``log0.sph`` **mean?**

With ``nblock=1``, every ``dtmaxblk`` the code writes lines like these to
``log0.sph``::

    block sync t=   9.8007812500000000       finest bin=           4  substeps=                   16  wake-ups=                    0  mean active=   2869.5000000000000       mean refreshed=   518.68750000000000
    bin histogram (bins 0..finest):     0     0     0  4193   752

The "mean refreshed" value (about 519 here) can be much larger than the number
of neighbours of any one particle (about 37 on average).  Here is where it
comes from.

*Two groups of particles.*  The bin histogram counts particles by step size,
from bin 0 (a step of ``dtmaxblk``, here 0.025) down to the finest bin (a step
of ``dtmaxblk/2**bin``).  At this moment the particles were in two groups:

* a fast group of 752 particles in bin 4, which need short steps of 0.0015625,
* a slow group of 4193 particles in bin 3, which can take steps twice as long,
  0.003125.

*Substeps.*  The code moves forward in the smallest step in use, 0.0015625.
Each of these is a substep, and 16 of them make up one ``dtmaxblk``.  The fast
group is updated (forces calculated, velocities changed) in every substep.  The
slow group is updated only in every other substep, because its step is two
substeps long.  A particle that is updated in a substep is *active*.  One that
is not is *inactive*: it keeps moving in a straight line at its current
velocity, but nothing is recalculated for it.

So the substeps alternate.  In the odd ones, the slow and fast steps end
together, so all 4945 particles are active and none is inactive.  In the even
ones, only the 752 fast particles are active, and the 4193 slow ones are
inactive.

*Why refreshing is needed.*  To calculate the force on an active particle, the
code uses its neighbours' densities and smoothing lengths.  If a neighbour is
inactive, those values were calculated a substep ago, and the neighbour has
moved since.  Using them would give a slightly wrong force.

*What "refreshed" means.*  In each substep, every inactive particle that is a
neighbour of an active particle has its density and smoothing length
recalculated at its current position.  That is a refresh.  A refreshed particle
stays inactive: no force is calculated for it and its velocity is not changed.
It just gets up-to-date values, so that its active neighbours get correct
forces.  Refreshing is always on.

*Where 519 comes from.*

* In the 8 odd substeps, everyone is active, so nothing needs refreshing: 0.
* In the 8 even substeps, the 752 fast particles are active.  About 1040 slow
  particles are close enough to at least one fast particle to be its
  neighbour, so those 1040 are refreshed.
* Averaged over all 16 substeps: (8 × 0 + 8 × 1040) / 16 ≈ 519.

*Why it is bigger than the neighbour count.*  "About 37 neighbours" is per
particle.  There are 752 fast particles scattered through part of the star,
and their neighbourhoods overlap only partly.  Counting every slow particle
that is a neighbour of at least one fast particle gives about 1040.

The other values on the ``block sync`` line: ``t`` is the time of the
velocities and internal energies, half the finest step after the time of the
positions (9.800 here), as in every dump.  ``finest bin`` is the smallest step
in use, ``substeps`` is how many substeps were taken since the last sync
(0.025/0.0015625 = 16), ``wake-ups`` counts particles whose step was cut short
because a neighbour suddenly needed a much shorter one, and ``mean active`` is
the average number of particles updated per substep.  Here
(8 × 4945 + 8 × 752) / 16 ≈ 2850 is close to the 2869.5 reported (a few
particles changed bins during the interval).  With shared timesteps, all 4945
particles would be updated in every substep.
