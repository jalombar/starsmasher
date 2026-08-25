Running a simulation
====================

Two files control a run.  ``sph.init`` chooses the kind of calculation, and
``sph.input`` sets everything else.  Both must be present: without ``sph.init``
the code stops with ``init: error reading input file sph.init`` and nothing more
helpful.

The minimum ``sph.init`` is three lines::

    &INITT
    INAME='1es'
    &END

See :doc:`../reference/sph_init` for the codes, and
:doc:`../reference/sph_input` for every setting and its default.

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
over only the gravity; neighbour finding and the SPH sums run on the CPU either
way, so the difference is smaller than you might expect.  Measured on one
machine, at four ranks with ``nnopt`` scaled as :math:`\sqrt{N}`:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - :math:`N`
     - GPU
     - CPU
     - CPU/GPU
   * - 1000
     - 3.6 s
     - 3.1 s
     - 0.86
   * - 16000
     - 35.9 s
     - 38.7 s
     - 1.08
   * - 64000
     - 200 s
     - 242 s
     - 1.21
   * - 128000
     - 543 s
     - 723 s
     - 1.33

Below about :math:`N=4000` the CPU build was the faster of the two.
Extrapolating the trend, it would take around :math:`N=6\times10^5` before it was
twice as slow.  Where that crossover falls is hardware dependent, and more cores
push it higher rather than lower, because one GPU serves the whole job however
many ranks you run.

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
   ``log0.sph`` and ``energy0.sph``; the first resume writes ``log1.sph`` and
   ``energy1.sph``, and so on.  A resumed run does not append to
   ``energy0.sph``, which matters if you are parsing an energy trace across a
   restart.
