Measuring the best ``ngravtree`` and ``ngpusph``
================================================

``ngravtree = -1`` and ``ngpusph = -1`` choose the gravity solver and the part of the SPH done on the GPUs from the number of particles, the threads, the ranks and the GPUs
(:doc:`gravity`, :doc:`gpu_sph`).  They never choose from a time, so that a run and its restart choose the same.  The rules come from measurements on three machines, and on
other hardware another pair of values can be faster.  ``parallel_bleeding_edge/tools/sph_autotune.py`` tries the pairs on your input and on your machine, for a few tens of
steps each, and tells you which two lines to write in ``sph.input``.

It is a plain Python 3 program (3.6 or later, nothing to install).

Running it
----------

In a folder that holds a run ready to start: ``sph.input``, and ``sph.init`` for a new run or ``restartrad.sph`` to restart from a dump (the code opens ``sph.init`` only when there is
no ``restartrad.sph``), the input files, and the executable::

   python3 /path/to/parallel_bleeding_edge/tools/sph_autotune.py --exe ./my_run_gpu_sph --threads 16

It runs the executable the way you give in ``--launch`` (by default ``mpirun -np 1 --bind-to none``: write your own launch command there, for example ``mpirun -np 2 --map-by
slot:PE=8 --bind-to core``), once for each candidate.  The options:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Option
     - Meaning
   * - ``--dir``
     - the folder of the run (default: the current folder)
   * - ``--exe``
     - the executable (default: the only ``*_gpu_sph``, ``*_cpu_sph`` or ``run_sph`` of the folder)
   * - ``--launch``
     - the command that starts it, as in your own runs; the executable is added at the end
   * - ``--threads``
     - OpenMP threads for each rank: sets ``OMP_NUM_THREADS`` (default: the environment as it is)
   * - ``--steps``, ``--warmup``
     - steps timed for each candidate (30) after the first steps that are not timed (5)
   * - ``--reps``
     - repetitions of each candidate (1); with more, they alternate in opposite orders and the spread between them decides what is equivalent
   * - ``--gravity``, ``--sph``
     - the values to try (``0,2`` and ``0,7``; for a build without GPU code ``0,1`` and ``0``); ``--gravity keep`` leaves the gravity of ``sph.input`` and tries the SPH only
   * - ``--work``, ``--keep``
     - where the candidates run (default: ``autotune_<date>_<time>`` in the folder of the run), and keep their folders
   * - ``--timeout``
     - seconds allowed to each run of a candidate (600)

What it does
------------

Nothing in your folder is changed.  Each candidate runs in a sub-folder of the work folder: small files are copied, large ones linked (the files that the code writes in its
folder, such as ``restartrad.sph``, are copied however large they are, and the ``out*.sph``, ``energy*.sph`` and ``log*.sph`` of earlier runs are not taken), and ``sph.input``
is written again with ``ngravtree`` and ``ngpusph`` set and the rest as it was, comments included.  The tool reads the standard output of the run, notes when each step ends,
stops the run (``mpirun`` and the ranks) after the warm-up and the steps, and removes the sub-folder.  One candidate runs at a time.  It also starts the code with
``ngravtree = -1`` and ``ngpusph = -1``, and with the gravity of your ``sph.input`` and ``ngpusph = -1`` (the default case), only as far as the line in ``log0.sph`` that says
what the code chose, so that the table says how far each of those choices is from the fastest pair.

``log0.sph`` of every candidate is read, to see what ran.  A candidate that the code refuses (``ngravtree = 2`` in a build without CUDA), or that it changes (``ngpusph = 7`` in a
periodic box, which the code does on the CPUs, and the log says so), is shown with the reason and left out of the ranking.

A candidate whose timed steps were not steady (their mean more than 10% above their median) is measured a second time: another program that takes the machine for a moment makes
some of the steps slow, and the candidate would lose for no reason of its own.  If the second run is steady it takes the place of the first, and the table says so.  If it is not
steady either, that is how the run goes (a step that writes a dump, neighbour lists that are built again every few steps, block time steps), and both runs count.

The example that follows is from a workstation with a Threadripper 7980X and a RTX 4070 Ti SUPER: a relaxed star of 99954 particles restarted from a dump, with ``--threads 63`` and the launch command shown in the second line (the paths are shortened)::

   sph_autotune: folder /data/star100k
     executable /data/build/test_gpu_sph
     launch "mpirun --oversubscribe -np 1 --bind-to none -x OMP_NUM_THREADS -x OMP_PLACES -x OMP_PROC_BIND <executable>", OMP_NUM_THREADS=63
     30 steps after 5 of warm-up, 1 repetition(s), work folder /data/star100k/autotune_20261113_123712

   [auto] ngravtree=-1 ngpusph=-1: only to read what the code chooses ...
         chose ngravtree=2 ngpusph=7
   [default] ngravtree=0 ngpusph=-1: only to read what the code chooses ...
         chose ngpusph=0
   GPU build; candidates: (0,0), (0,7), (2,0), (2,7)

   [1/4] ngravtree=0 ngpusph=0 ...
         16.02 ms a step (median 16.15, 30 steps, 2 s)
   [2/4] ngravtree=0 ngpusph=7 ...
         22.58 ms a step (median 22.72, 30 steps, 2 s)
   [3/4] ngravtree=2 ngpusph=0 ...
         13.88 ms a step (median 13.85, 30 steps, 2 s)
   [4/4] ngravtree=2 ngpusph=7 ...
         14.95 ms a step (median 14.87, 30 steps, 2 s)

   Result of sph_autotune for /data/star100k, 99954 particles, 1 rank(s) with 63 thread(s) (from log0.sph)

   candidate                 steps   mean ms   median      min      max  vs best    status  note
   ---------------------------------------------------------------------------------------------
   ngravtree=0 ngpusph=0        30     16.02    16.15    14.23    20.37   +15.4%  measured  <- your sph.input, chosen by ngpusph=-1 with your gravity
   ngravtree=0 ngpusph=7        30     22.58    22.72    18.54    28.54   +62.6%  measured
   ngravtree=2 ngpusph=0        30     13.88    13.85    12.23    15.53  fastest  measured
   ngravtree=2 ngpusph=7        30     14.95    14.87    12.71    21.44    +7.7%  measured  <- chosen by the automatic values

   The automatic values (ngravtree=-1, ngpusph=-1): the code chooses ngravtree=2 ngpusph=7: 14.95 ms against 13.88 for the fastest, 7.7% slower.
   The default case (ngravtree=0 as in your sph.input, ngpusph=-1): the code chooses ngravtree=0 ngpusph=0: 16.02 ms against 13.88 for the fastest, 15.4% slower.

   Fastest: ngravtree=2, ngpusph=0 (13.88 ms a step).

   The fastest uses the tree, and your sph.input has the direct sum.  The tree is an approximation (treeacc), so the forces differ from the direct sum by the error of
   the approximation, not only in the last digits.  The fastest with the gravity of your sph.input:
   To keep the same gravity nothing needs to change: your sph.input as it is now (ngravtree=0, ngpusph=-1, which the code runs as ngravtree=0, ngpusph=0) gives 16.02 ms a step, the fastest.
   The file has -1, so the code chooses this pair itself; to fix it by hand (the choice can change on another machine or with another number of threads) write:
       ngravtree = 0
       ngpusph = 0
   To take the fastest and accept the other gravity solver write:
       ngravtree = 2
       ngpusph = 0

   The difference between the fastest overall and the next candidate is 7.7%: with one repetition a candidate can change by several percent from one run to the next; run again with --reps 3 to be sure.

   The values written by hand stay as they are, so a run is reproducible.  sph_autotune changed nothing in your folder or in the code: write the lines yourself.
   The automatic values (-1) choose again from the particles, threads, ranks and GPUs, and may choose differently on another machine.
   This is the cost of the first 35 steps of your run; measure again if the cost changes during the run (a collision), or if the machine, the threads or the particles change.

   Your folder is unchanged (name, size and modification time of 2 files, before and after).
   Each candidate had 14.1 MB copied and 0 large files linked; the folders of the candidates were removed.  Time of sph_autotune: 11 s.

   Written: /data/star100k/autotune_20261113_123712/autotune.txt and /data/star100k/autotune_20261113_123712/autotune.csv

``autotune.txt`` and ``autotune.csv`` in the work folder hold the same table.

How to read the advice
----------------------

* The pair that your ``sph.input`` gives now is always measured and ranked with the others, even when it is not one of the candidates (``ngpusph = 3``, ``ngravtree = 1`` in a GPU build).  With
  ``-1`` in the file it is the pair that the code chooses for those values, read from ``log0.sph``.  The table marks it ``<- your sph.input``; if the code refuses or changes it, the table says why.
* The ranking is by the **mean** step time of the timed steps (the median is shown too).  Two candidates within 3% of each other, or within the spread of the repetitions if
  ``--reps`` is more than 1, are equivalent.  Among the candidates equivalent to the fastest the advice is, in this order: the pair that your ``sph.input`` gives now, then the pair that
  ``ngravtree = -1`` and ``ngpusph = -1`` choose, then the fastest.  If you already have an equivalent pair the tool says that there is nothing to change: you keep the same digits as before.  When the
  file has ``-1``, it also gives the two lines that fix by hand what the code chooses now.
* The tree is an approximation (``treeacc``): if the fastest candidate uses the tree and your ``sph.input`` has the direct sum, or the other way round, the forces change by the
  error of the approximation, not only in the last digits.  The tool then gives two answers: the fastest with the gravity of your ``sph.input``, and the fastest of all.
  Choosing between them is yours.
* ``ngpusph`` 0 and 7 do not give the same last digits of the forces (:doc:`gpu_sph`).  Both are right, but the run that follows is not the same, bit for bit, as the run you had.
* The values you write stay as they are, so the run is reproducible.  The tool changes nothing, in your folder or in the code: you write the two lines.

Limits
------

* It times the first steps of the run it is given.  When the cost changes with time (a collision: the star is a different problem at the pericentre), the best pair can too:
  start the tool from a restart file taken where the run spends its time.
* Measure again if the machine, the number of threads, the number of ranks or the number of particles changes.  The setup of the run (reading the dump, or making the star)
  is repeated for each candidate, so use a restart file for a run that would otherwise start from a stellar profile.
* With block time steps (``nblock = 1``) a step of the tool is a sub-step, and the cost of a sub-step depends on how many particles are active: take ``--steps`` of a few hundred, so that
  the steps cover whole cycles.  Steps that differ this much are not steady for the tool, so a candidate is then measured twice and both runs count.
* Dumps written during the timed steps (a small ``dtout``) add the same time to every candidate, so the differences look smaller than they are.
* With one repetition a candidate can change by several percent from one run to the next (on the workstation of the example, with 63 threads, the same candidate gave from 15.0 to 16.5 ms a step in
  runs made one after the other).  When the advice changes your values for a gain below 10%, or the nearest candidate that is not equivalent is less than 10% behind the one advised, the tool says
  that the difference is not sure and proposes ``--reps 3``.
* The tool does not know what else the machine is doing.  It measures a candidate again only when some of its steps were much slower than the others.  Another program that runs during all the
  timed steps of a candidate slows them all alike, and the tool does not notice: in a test, one second of another load on every core made a candidate look twice as slow as it is, with no
  warning.  A load that raises the mean by less than 10% is not noticed either.  Run the tool when the machine does nothing else.  If the machine cannot be quiet use ``--reps 3`` and look at
  the spread in the table: a large spread says that a repetition was disturbed.
* It stops the run it started when the steps are done or at ``--timeout`` (and on Ctrl-C, which also removes the sub-folders), but ``kill -9`` on the tool itself leaves the run it was timing.
