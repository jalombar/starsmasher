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
