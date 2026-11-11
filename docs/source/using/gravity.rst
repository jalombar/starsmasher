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
