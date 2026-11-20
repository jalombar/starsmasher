# StarSmasher test suite

Short runs that check StarSmasher against things that must be true, rather than
against stored output. A golden-file test tells you something changed; these
tell you whether the answer is *right*.

Each test states a physical or internal-consistency expectation and fails with
the numbers when it is not met.

## Running

```
python3 tests/run_tests.py                 # everything
python3 tests/run_tests.py --list          # what exists
python3 tests/run_tests.py virial energy   # named tests only
python3 tests/run_tests.py --np 2          # fewer MPI ranks
python3 tests/run_tests.py --hdf5          # build with make HDF5=1, so that hdf5_dump runs
```

By default the suite builds the code once into a scratch directory and reuses
it. Pass `--exe /path/to/..._sph` to test an executable you have already built.

## The tests

| Name | Checks |
|---|---|
| `virial` | A relaxed n=1.5 polytrope satisfies `W = -3/(5-n) GM^2/R` and `U = -W/2`. |
| `energy` | Total energy is conserved over the run to within a tolerance. |
| `quiet` | A relaxed star stays put: kinetic energy stays negligible. |
| `ranks` | The answer does not depend on how many MPI ranks are used. |
| `cpu_gpu` | The CPU and GPU builds agree. |
| `block_cpu_gpu` | With block timesteps, where gravity is computed on the active particles only, the CPU and GPU builds agree. |
| `gpu_hilo` | GPU gravity with and without the float remainders of the positions (`ngravhilo`) agrees for a single star. |
| `gpu_tree` | The tree gravity on the GPU (`ngravtree=2`) agrees with the tree on the CPUs (`ngravtree=1`) for the same opening angle and leaf size: energies, accelerations and potentials of a relaxed star; a second GPU run divides most leaves in chunks (`STARSMASHER_TREE_HEAVY`, below). |
| `gravity_auto` | `ngravtree=-1` chooses the gravity from the number of particles (GPU build: the direct sum below `ngravtreeauto`, the tree walked on the GPU from it up; CPU build: the tree on the CPUs), says so in `log0.sph`, and gives exactly the energies of the value it chose, set by hand. |
| `gpu_sph` | `ngpusph=1` (the smoothing lengths, densities and neighbour lists on the GPU) gives exactly the bits of `ngpusph=0` (energies of every row, bytes of every dump), on one rank and on several; with `ngpusphcheck=1` the CPUs and the GPU agree particle by particle; the CPU build says it has no GPU code and gives the same bits. |
| `gpu_forces` | `ngpusph=3` (the pair sums of the forces on the GPU too) agrees with `ngpusph=1` to the rounding of the sums, with shared and with block steps, on one rank and on several: energies of every row to one unit of the last digit that `energy0.sph` prints (7 digits), and with `ngpusphcheck=2` accelerations to 1e-12 and du/dt to 1e-10 of the largest; the CPU build says that it has no GPU code and gives the same bits. |
| `gpu_ranks` | `ngpusphcheck=2` and `3` with more ranks than GPUs (`CUDA_VISIBLE_DEVICES=0`, 2 and 3 ranks): the ranks without a GPU of their own take part in the comparison of the forces, the run finishes (it used to hang) and agrees with `ngpusph=1` to 1e-9. |
| `gpu_resident` | `ngpusph=7` (the candidates of the smoothing-length search kept on the GPU between steps, Verlet lists, with a small margin so that they are built again in the run) gives exactly the bits of `ngpusph=3` (energies of every row, bytes of every dump), on one rank and on several, also with `ngpusphcheck=3`; the CPU build says that it has no GPU code and gives the same bits. |
| `gpu_mixed` | `ngpusphprec=1` (the terms of the pairs of the forces on the GPU in single precision, sums in double) stays within the tolerances: with `ngpusphcheck=2` the accelerations differ from the CPU's by less than 1e-4 and du/dt by less than 1e-3 of the largest at every step, and the energies of the run are within 1e-4 of those in double precision; the CPU build says that it has no GPU code and gives the same bits. |
| `sph_auto` | `ngpusph=-1` chooses from the configuration (the CPU build 0, a GPU build with the tree on the CPUs 7, with the direct sum by threads and particles), writes the choice and the reason in `log0.sph`, and gives exactly the bits of a run with the value it chose set by hand, on one rank and on several. |
| `vdotsm` | With `cn2` and `cn4` finite (so that `tstep` uses the smoothed accelerations of `vdotsm`, which sums over the neighbour lists on the CPU) `ngpusph=1` gives exactly the bits of `ngpusph=0` and `ngpusph=7` those of `ngpusph=3` (energies of every row, bytes of every dump, number of steps), on one rank and on several; the same run with the default `cn2` and `cn4` comes out different. |
| `gpu_multidev` | The SPH and the forces with several devices for the rank, tried on one GPU with `SPHGPU_NDEV` (logical devices dealt out round the GPUs): with 1, 2 and 4 logical devices the run reaches the end, with `ngpusphcheck=2` accelerations to 1e-12 and du/dt to 1e-10 of the largest, and the energies of every row are those of one device to 1e-9. |
| `guard_neighbors` | Exhausting the neighbour list stops the run with a message rather than corrupting memory. |
| `guard_ngravprocs` | Asking for more gravity processes than there are GPUs is clamped, not fatal. |
| `periodic` | A fully periodic box conserves momentum and energy and has no density seam at its faces. |
| `periodic_block` | Block timesteps in a periodic box agree with shared ones and keep every particle in the box. |
| `guard_periodic` | A periodic box too small for the kernel stops the run with a message. |
| `hdf5_dump` | `nhdf5=1` writes an `out*.h5` that equals `out*.sph` bit for bit, and `nhdf5=2` the same files without `out*.sph`. Skipped unless the build has HDF5 (`make HDF5=1`, or `--hdf5`) and `h5py` is installed. |

The last four exist because each corresponds to a real defect: results that
depended on rank count, a CPU/GPU disagreement, an unbounded write past the end
of an array, and a segfault when `ngravprocs` exceeded the GPU count.

## Running the suite with another solver

`STARSMASHER_TEST_EXTRA` adds entries to the namelist of every run, separated by
commas; the entries a test sets itself win. For example

```
STARSMASHER_TEST_EXTRA="ngravtree=1,treeacc=0.45" python3 tests/run_tests.py
```

runs the whole suite with the gravity tree on the CPUs. The reference values and
the tolerances of the tests do not change, so this checks that the other solver
passes the same physical and consistency tests as the default one. The
exceptions are `ngravtree=2` and `ngravtree=-1`: the GPU build then runs the tree on
the GPU (with `-1`, below `ngravtreeauto`, the direct sum) and the CPU build the tree on
the CPUs, and `cpu_gpu` and `block_cpu_gpu` compare the two builds with the limits 5e-6
and 2e-4 instead of 1e-6 and 1e-5 (the worst differences measured, with the star of 2000
and of 10000 particles, are 1.2e-6 and 5.6e-5).

## Environment variables of the tests

`STARSMASHER_TEST_EXTRA` is described above. `STARSMASHER_TREE_HEAVY` is read by the
tree gravity on the GPUs (`ngravtree=2`): the number of entries of the list of near
leaves above which a leaf is divided in chunks (768 by default). The `gpu_tree` test
sets it to 24 for its second GPU run, so that most leaves of its small star take the
path of the leaves with a very long list, which an ordinary model never takes.
