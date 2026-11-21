Choosing the gravity solver
===========================

By default the gravity of the particles is the direct sum over all pairs.  It
runs on the GPUs of an executable built with ``make`` and on the CPUs of one
built with ``make cpu``.  ``ngravtree`` in ``sph.input`` selects another solver.

Gravity with a tree on the CPUs
-------------------------------

``ngravtree = 1`` computes the gravity with a tree on the CPU cores.  Its cost
grows as :math:`N \log N` instead of :math:`N^2`, and it uses the threads of
every rank (see :doc:`../installation` for the OpenMP build), so with many
particles and few or no GPUs it is much faster than the direct sum.  It needs
``nselfgravity = 1``, which the code sets for you.

The tree is a k-d tree whose leaves hold at most ``nleafgrav`` particles.  It is
built and walked as the tree of Phantom is (Price et al. 2018, PASA 35, e031):
a leaf collects the nodes that are far enough away, replaces each by its
quadrupole in one Taylor series of the force about the centre of the leaf, and
sums the near particles directly with the pair formula of the direct sum.  Two
settings control it, and both are in the reference of ``sph.input``:

* ``treeacc`` is the opening angle :math:`\theta`.  A node is used as a whole if
  :math:`\theta^2 r^2 > (s_A + s_B)^2` and if every pair between the two nodes
  is farther apart than the softening lengths, so that the force between them is
  Newtonian.  Smaller values are more accurate and slower, and ``treeacc = 0``
  gives the direct sum.  The error of the force depends on the model: a star
  that is centrally condensed needs a smaller value than a uniform lattice.
* ``nleafgrav`` is the largest number of particles of a leaf.

The result does not depend on the number of threads (the same bits for any
``OMP_NUM_THREADS``), as in the rest of the code.  The tree, unlike the direct
sum, is not antisymmetric in the forces it gives: the total momentum and angular
momentum of the particles drift slowly, while the total energy is conserved as
well as with the direct sum.  Use the direct sum where they must stay exact.

Things to know when running it:

* It does not use a GPU.  Use ``make cpu`` on a node without one; the
  executable of ``make gpu`` also runs it, but stops at start-up when no GPU is
  visible.
* Every rank builds the whole tree.  With the OpenMP build, one rank with all
  the cores (or one per socket) is the right layout.
* At most ``ngravprocsmax`` (32) ranks compute gravity; the others only do the
  hydrodynamics.
* ``gravtree.f90`` and ``gravpair.f90`` are compiled with ``-mavx2 -mfma
  -fno-trapping-math`` on x86-64, and nothing else is, so the default gravity
  keeps its bits.  ``make TREEFLAGS=`` (empty) builds them for a processor
  without AVX2.

Overlapping the GPUs with the smoothing-length solve
----------------------------------------------------

With the direct sum on the GPUs the card is idle while the CPUs solve the
smoothing lengths.  ``ngravoverlap = 1`` uses that time.  Before the lengths are
solved, the GPU starts the tiles of particles that are far from each block
whatever the lengths turn out to be (almost all the pairs of a large model);
when they are solved the code checks that none has grown past the limit that
was assumed, and the GPU does the remaining tiles and adds the two parts in a
fixed order.  If a length did grow past its limit, which is rare, the gravity of
that step is done again in the usual way and a line is written to the standard
output.

The result is deterministic, but not the same bits as ``ngravoverlap = 0``: the
force sums are grouped differently, which changes them in the last digits.  The
energy is conserved as well as without it.  It is not used with block time steps
(``nblock = 1``) or with the tree on the CPUs.

The default, ``ngravoverlap = -1``, chooses once at the start: 1 for the direct
sum on one GPU and for the tree on the GPUs, and 0 for the direct sum on more than
one GPU (``ngravprocs`` above 1, or ``ngpuperrank`` other than 1), where the
overlap made the step longer for 100k to 500k particles (it made it shorter at 2M,
where the tree is the better solver, which ``ngravtree = -1`` takes), and 0 where
there is no CUDA library (``make cpu``, the Metal build).  The choice is written in the
log of rank 0.  ``ngravoverlap = 0`` and ``ngravoverlap = 1`` set it by hand;
``ngravoverlap = 0`` gives back the bits of the direct sum without the overlap.

Several GPUs for one process
----------------------------

By default each process that does gravity (up to ``ngravprocs`` of them) uses one
GPU.  ``ngpuperrank = N`` lets a process use ``N`` GPUs: the process with local
number *r* on its node uses the GPUs *r*\ N to *r*\ N + N - 1.  Each of them sorts
all the particles in the same order, takes a run of whole blocks of the sorted
list and sums over all the other particles with the same division as a single
GPU, so **the forces have the same bits for any** ``N``, however the blocks are
shared.  The blocks are shared in proportion to the speed each GPU showed in the
last calls, which moves only the borders between the runs; a warning is written
when the GPUs of a process are out of balance (the slowest more than 10% above
the mean) in three calls in a row.  ``ngpuperrank = 0`` divides all the GPUs of
the node among the gravity processes.

This is the way to use all the GPUs of a node from a single process that has all
the cores (``-np 1`` with the OpenMP build): the preparation of the particles
and the sort are done by every GPU, so the gain is smaller than the number of
GPUs for a small number of particles.  It works with ``ngravoverlap = 1``.  The
tree gravity on the CPUs does not use the GPUs, and the block time steps use
the first GPU of the process only.

The tree on the GPUs
--------------------

``ngravtree = 2`` is the same tree gravity as ``ngravtree = 1`` (the same opening
criterion, the same Taylor series of the far field about the centre of mass of each
leaf, the same pair formula for the leaves that are near), with the visit of the
tree done by the GPUs.  The CPUs make the tree and its moments, and the GPUs do the
rest: a warp for each leaf, a thread for each particle.  It needs a build with CUDA
(``make cpu`` and the Metal build stop with a message) and ``nleafgrav`` of at most
32; ``treeacc`` and ``nleafgrav`` mean the same as for ``ngravtree = 1``, and the
default is still ``ngravtree = 0``, the direct sum.

The forces differ from those of ``ngravtree = 1`` with the same ``treeacc`` and
``nleafgrav`` only by the rounding of the float arithmetic of the GPU, which is far
below the error of the tree itself.  The force on a particle depends on the tree and
on nothing else: **the same bits for any number of GPUs** (``ngpuperrank``), **of
processes and of threads**, and with block steps the active particles have the bits
they have when all are active.  (The dumps of two runs with a different number of
processes still differ, with every gravity solver, because the sums of the SPH part
depend on it.)

A leaf that holds a particle whose softening length is larger than a star (the
particles thrown out of a collision) cannot use any node near it: its list is the
whole system.  Such a leaf is divided in chunks of its list that different warps, and
different GPUs, do, and the sums of the chunks are added in a fixed order, so that the
bits do not depend on who did what; the other leaves are done as they would be
otherwise.  With several GPUs for one process the leaves are divided among the GPUs in
proportion to the speed each showed in its first pass.

``ngravoverlap = 1`` with ``ngravtree = 2`` makes the tree and queues the visit on the
GPUs as soon as the smoothing lengths are solved, and collects the forces after the
loop over the pairs, so that the CPUs and the GPUs work together.  The bits are those
of ``ngravoverlap = 0``.

As the tree on the CPUs, the tree on the GPUs does not conserve the momentum and the
angular momentum exactly, and conserves the energy as well as the direct sum does.  Use
one process with the OpenMP build: every process builds the whole tree, so one process
for each GPU repeats the work.

Letting the code choose
-----------------------

``ngravtree = -1`` chooses the solver once, when the number of particles N is known, as the fastest one, and writes the choice in the log of rank 0, for example::

   gravity: ngravtree=-1 (auto): N=499770 >= 9966 (the limit measured with ngpusph=-1 and 63 threads per rank) and a GPU -> ngravtree=2 (tree walked on the GPU)

The rule is:

* no GPU (the CPU build, or no GPU on the node): the tree on the CPUs (1);
* no gravity (``ngr = 0``), no self-gravity or an invalid tree setting: the direct sum, with the reason in the log, rather than a stop;
* a GPU that cannot walk the tree (the Metal build), or ``nleafgrav`` above 32: the direct sum;
* a GPU and ``ngravtreeauto`` written (0 or more): the direct sum below ``ngravtreeauto`` particles and the tree walked on the GPU (2) from it up, whatever the SPH;
* a GPU and ``ngravtreeauto = -1`` (the default): the tree walked on the GPU (2) from the number of particles at which it is faster than the direct sum, and the direct
  sum below it.  The limit depends on where the SPH is done (:doc:`gpu_sph`), because the CPUs and the GPU share the step:

  ================================================== ======================================
  the SPH                                            the tree from N =
  ================================================== ======================================
  ``ngpusph = -1`` as well, or 3 or 7                9966
  ``ngpusph = 1``                                    25000
  ``ngpusph = 0`` (or forced to 0)                   75000 with 48 threads or more per rank, 450000 with fewer
  ================================================== ======================================

  With ``ngpusph = -1`` as well the two are chosen together, as the fastest pair: the tree on the GPU with the SPH on the GPU (``ngpusph = 7``) wins from 9966 particles on, with
  7 to 63 threads per rank, by 5% to 16% at 9966 particles, 12% to 58% at 99954 and 2.9 to 4.7 times at 499860 (the best direct sum against the tree with the SPH on the GPU).  With two GPUs (two ranks, or one rank with ``ngpuperrank = 2``) the tree
  wins from about 15000 particles with the SPH on the GPU (10% at 19938, 2 ranks of 8 threads) and from 450000 with ``ngpusph`` 0 or 1; with three or more GPUs (ranks times ``ngpuperrank``) the tree wins from 9966 particles with one rank and
  from 99954 with several ranks, with every ``ngpusph`` (each rank builds the whole tree, which is why it wins later with several ranks); with ``ngpuperrank = 0`` (not measured) the limit stays 299000.  Below 9966 particles, the smallest star measured, the choice stays what it was in the earlier versions: the direct sum.

The limits come from the step times of the direct sum and of the tree on a GPU, with ``ngpusph`` 0, 1 and 7, for relaxed Suns of 9966 to 499860 particles, 7, 16, 32 and 63 threads
per rank, one rank and one GPU, measured on two kinds of node: a Threadripper 7980X with a RTX 4070 Ti SUPER, and an EPYC 9335 with a RTX PRO 6000; the limits
of the SPH with the direct sum were then moved after a check on a third node, a Xeon 8480C with a H100 (see :doc:`gpu_sph`).  The limits of the tree above do not depend on how
fast the GPU is; the limits of the SPH with the direct sum on one GPU, with 32 threads or more, do, through the class of the GPU (mid or strong, from its number of multiprocessors, ``ngpuclass``: see :doc:`gpu_sph`),
because the nodes do not agree (with 32 threads and 100000 particles the direct sum and the CPUs for the SPH win on the first, the GPU for the SPH on the second and the third).  Each limit is
the number that gives the smallest loss in the worst cell of the machines it was measured on, and the number of multiprocessors is only a coarse hint of the speed of a GPU.  On other hardware another value may be better:
write ``ngravtreeauto`` then (and ``ngpuclass``, or ``ngpusph``, for the SPH).

With three or more GPUs the limits come from relaxed Suns of 9966 to 499860 particles and the collision of two Suns (199908 particles) at four moments, with 7 to 63 threads per rank, on 4 H100 (a node with two
Xeon 8480C: one rank with ``ngpuperrank = 4``, 4 ranks with a GPU each, 2 ranks with two) and on 8 RTX PRO 6000 (a node with an EPYC 9335: one rank with ``ngpuperrank = 8``, 8 ranks with one GPU each, 4 ranks with two, 2 ranks
with four), mean step of the window, one repetition (three in some cells).  With one rank and 63 threads the tree with ``ngpusph = 7`` takes 6.34 ms at 99954 particles on 4 H100 against 11.25 for the direct sum with 7
(and 11.41 with 1, the choice of the earlier versions), and 17.34 against 53.57 at 499860 (8 RTX PRO 6000: 16.91 against 27.06).  With 4 ranks of 15 threads the direct sum is ahead at 30090 particles (3.22 ms against 3.51 for the tree) and the tree from
99954 (8.72 against 9.11), and with 8 ranks of 7 threads the direct sum is ahead up to 199908 (12.86 against 14.01), because every rank builds the whole tree.  The limits are numbers of measured stars, not interpolated: 9966 and 99954.
Cells where the choice of the new limits still takes more than 10% longer than the fastest pair the mode leaves free (4 or 8 GPUs, one repetition unless said): one rank, 16 threads, 4 GPUs, 9966 particles, the tree with 7
takes 4.25 ms against 3.54 for the direct sum with 7 (three repetitions); 8 GPUs, 32 or 63 threads, 9966 particles, the SPH on the CPUs is faster (the tree with ``ngpusph = 0``: 2.46 and 2.71 ms against 3.27 and 3.99 with 7), and at 30090
particles with 63 threads ``ngpusph = 1`` is faster (the tree: 3.90 ms against 4.33 with 7); 8 ranks of 7 threads, 99954 particles, the direct sum takes 5.96 ms against 7.29 for the tree (three repetitions); 2 ranks of 2 GPUs with 31 threads,
30090 particles, the tree takes 3.23 ms against 3.61 for the direct sum, and at 99954 the direct sum takes 10.29 ms with 1 against 11.89 with 7.  With three GPUs, or with ``ngpuperrank = 0``, nothing was measured: three GPUs take the numbers
of the measured ones.

The choice depends on nothing but N, the threads per rank, the ranks, the GPUs and the options, never on a time measured during the run, so a run and its restart choose the same, and the
results are bit for bit those of the value chosen, set by hand.  The part of the SPH done on the GPUs is chosen the same way by ``ngpusph = -1`` from the gravity chosen here and the
threads per rank (see :doc:`gpu_sph`).  ``ngravoverlap = -1`` (the default) then takes the value that suits the choice, 1 for the tree on the GPU.  With several MPI processes the rule is the
same, but the tree is built by each of them, so one process is advised.  The default stays ``ngravtree = 0``, which keeps the results of earlier versions.
