Running a simulation
====================

Two files control a run.  ``sph.init`` chooses the kind of calculation, and
``sph.input`` sets everything else.  Both must be present.

The minimum ``sph.init`` is three lines, for example::

    &INITT
    INAME='erg'
    &END

``INAME`` is a three-letter code naming the setup routine.  Two of them cover
most work:

``erg``
   Build a star from a stellar-evolution profile, a MESA model for instance,
   and relax it into an SPH star.  This is where nearly every project starts, and it
   is the subject of :doc:`../tutorials/relaxing_a_star`.

``hyp``
   Put two already-relaxed bodies on a Keplerian orbit and let them meet: a
   collision or fly-by between two stars, or between a star and a compact
   object.  See :doc:`../tutorials/a_collision`.

The usual sequence is therefore two runs: an ``erg`` run per star to make it,
then one ``hyp`` run to collide them.

.. dropdown:: The other thirteen initialization scripts
   :icon: list-unordered

   .. list-table::
      :header-rows: 1
      :widths: 12 88

      * - Code
        - Sets up
      * - ``1es``
        - a single polytrope, from ``starmass`` and ``starradius``
      * - ``1mc``
        - a single polytrope with a compact core
      * - ``2cr``
        - a corotating binary
      * - ``bps``
        - a binary plus a single star
      * - ``bph``
        - a binary plus a black hole
      * - ``hbs``
        - a binary encountering a single star
      * - ``meq``
        - a star whose particle-mass scheme varies with position
      * - ``res``
        - a rescaling of an existing model
      * - ``tri``
        - a triple system
      * - ``bhe``
        - an encounter with a supermassive black hole
      * - ``grs``
        - a general-relativistic setup
      * - ``rin``
        - a restart from an existing dump
      * - ``txt``
        - particles laid out from an ASCII image

   :doc:`../reference/sph_init` gives the same list with the routine each one
   calls.

Each code selects one initialization script, and they live in files named
``initialize_*.f`` in ``parallel_bleeding_edge/src``: ``initialize_parent.f``
for ``erg``, ``initialize_hyperbolic.f`` for ``hyp``, and so on.

Which script you choose also decides which *other* files the run needs: a
stellar-evolution profile for ``erg``, a relaxed star or two for ``hyp``.  See
:doc:`input` for what those are and where they come from.

The important thing about these routines is how little they do.  A setup routine
runs once, before the first step, and its only job is to decide where the
particles start: their positions, masses, velocities and internal energies.  It
then hands that state to the integrator and takes no further part.  Everything
after that is the same code regardless of which ``INAME`` you chose: the
hydrodynamics, the gravity, the timestep control, and what gets written to
disk.

Two consequences are worth keeping in mind.  A setting that controls the
*evolution* means the same thing for every ``INAME``, so advice about ``cn1`` or
``nnopt`` carries across unchanged.  A setting that describes the *starting
state* is read only by the routine that needs it, so ``e0`` does nothing in a
relaxation and ``equalmass`` does nothing in a collision, and neither will warn
you.

See :doc:`../reference/sph_input` for every setting and its default.

Choosing the number of ranks
----------------------------

::

    mpirun -np 8 ./my_run_gpu_sph

Worth testing on your own machine, but sixteen ranks is often a sensible
ceiling.  Beyond that, communication costs eat into the gains, and leaving cores
free lets you run several simulations at once.

``ngravprocs`` sets how many ranks compute gravity.  It cannot exceed the number
of GPUs present: ask for more and it is clamped to what you have.  Leave it at
0 to have it detected automatically.

GPU or CPU
----------

``make`` builds the GPU version and ``make cpu`` the CPU one.  The GPU takes
over only the gravity: neighbour finding and the SPH sums run on the CPU either
way.  That is why the two are indistinguishable on a small model and diverge
sharply on a large one.  Gravity is the part that grows fastest with the number
of particles, and it is the only part the card touches.

What was measured
~~~~~~~~~~~~~~~~~

Relaxing a one-solar-mass polytrope, ``INAME='erg'`` aside, with::

    &input
    n=N, nnopt=22, tf=1.0, dtout=1000.0, nrelax=1, treloff=0,
    starmass=1.0, starradius=1.0,
    &end

on four MPI ranks, on one machine: an RTX 5080 with the CPU build using all
four ranks for gravity.  Both builds were checked to produce the same answer,
with identical ``W``, ``U`` and total energy to every digit written, so the two
columns are the same work, not merely similar work.

.. list-table::
   :header-rows: 1
   :widths: 14 14 18 18 16 20

   * - :math:`N`
     - Iterations
     - GPU, per step
     - CPU, per step
     - CPU/GPU
     - CPU, whole run
   * - 1 000
     - 23
     - 82 ms
     - 75 ms
     - 0.9
     - 1.7 s
   * - 4 000
     - 30
     - 68 ms
     - 69 ms
     - 1.0
     - 2.1 s
   * - 16 000
     - 43
     - 68 ms
     - 171 ms
     - 2.5
     - 7.4 s
   * - 64 000
     - 64
     - 129 ms
     - 1817 ms
     - 14
     - 116 s
   * - 128 000
     - 80
     - 225 ms
     - 7104 ms
     - 32
     - 568 s

The per-step columns are the comparison worth making.  The iteration count is
not a constant of the experiment: every run above stops at the same ``tf``, but
a bigger model takes shorter timesteps and so needs more of them, 23 steps at
:math:`N=1000` against 80 at :math:`N=128\,000`.  A table of total wall times
alone would mix that in with the cost of a step and show neither clearly.

What this means in practice
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below about :math:`N=4000` it makes no difference.  The GPU is fractionally
slower, because there is not enough work to cover the cost of moving it to the
card.  By :math:`N=16\,000` the CPU is taking two and a half times as long, and
by :math:`N=128\,000` it is thirty times.

So a relaxation of a few thousand particles is perfectly comfortable without a
card, and the tutorials in this documentation are sized to run either way.  A
production model of :math:`10^5` particles is not: the same job that occupies a
GPU for a few minutes will occupy the CPU build for hours.

These figures are one machine and one kind of run.  The protocol is written
above so you can repeat it on yours rather than trusting the numbers, and it is
worth doing before committing to a long job.

.. warning::

   If you have built both versions, ``./*_sph`` matches two files and the shell
   hands you the alphabetically first, which is the CPU build.  It runs, just
   more slowly, and nothing says so.  The give-away is the first lines of
   output: the GPU build reports ``SPHgrav found 1 CUDA devices`` and
   ``is running on ... with gpu 0``, and the CPU build reports neither.

Restarting
----------

A ``restartrad.sph`` file is written every few iterations, overwriting the
previous one.  If one is present when a run starts, it is used automatically.
You can also restart from any snapshot by renaming it to ``restartrad.sph``.

.. note::

   Output files are numbered per stage, not per rank.  A fresh run writes
   ``log0.sph`` and ``energy0.sph``, the first resume writes ``log1.sph`` and
   ``energy1.sph``, and so on.  A resumed run does not append to
   ``energy0.sph``, which matters if you are parsing an energy trace across a
   restart.
