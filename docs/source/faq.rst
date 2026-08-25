Frequently asked questions
==========================

**The build stops with** ``cannot find -lcudart``.

Your CUDA installation is not where the build is looking.  Run ``make config``
to see what it found, then ``make CUDAPATH=/path/to/cuda``.

**The build stops with** ``Unsupported gpu architecture``.

Your CUDA release no longer supports the architecture being requested.  Recent
releases have dropped older ones -- CUDA 13, for instance, no longer accepts
``sm_61``.  The gravity library normally detects your card from ``nvidia-smi``;
if that is unavailable, pass the compute capability yourself with
``make COMPUTE_CAPABILITY=75``.

**I have no NVIDIA card.  Can I still use StarSmasher?**

Yes.  ``make cpu`` builds a version that needs no GPU.  It is not as slow as you
might expect, because the GPU only takes over the gravity; see
:doc:`using/running` for measured timings.

**The run stops immediately with** ``init: error reading input file sph.init``.

A run needs ``sph.init`` as well as ``sph.input``.  It is three lines; see
:doc:`reference/sph_init`.

**The run stops with** ``Cannot match namelist object name``.

``sph.input`` sets a variable that is not in the namelist.  Check it against
:doc:`reference/sph_input`.

**It segfaults as soon as it starts, in the CUDA library.**

Usually ``ngravprocs`` exceeding the number of GPUs.  Current versions clamp it
and say so; older ones did not.  Leave ``ngravprocs`` at 0 to have it detected.

**My relaxed star has particles flying off.**

Often the particle-mass distribution.  For a star spanning many orders of
magnitude in density, ``equalmass=0`` gives very light outer particles that never
settle.  See :doc:`tutorials/relaxing_a_star`.

**How many neighbours does** ``nnopt`` **give me?**

More than ``nnopt``: the actual count is larger by roughly 1.4 to 1.7 depending
on the model.  ``nnopt`` is the target of the smoothing-length constraint, not a
neighbour count.  See :doc:`reference/equations_of_motion`.

**Which of** ``parallel_bleeding_edge`` **and** ``Blackollider`` **should I use?**

``parallel_bleeding_edge`` unless you are modelling collisions with compact
objects treated as point masses.  The two differ mainly in how the smoothing
length is related to the density.

**My restarted run behaves differently from an uninterrupted one.**

Check which files you are reading.  Output is numbered per stage, not per rank: a
resumed run writes ``log1.sph`` and ``energy1.sph`` rather than appending to the
originals.
