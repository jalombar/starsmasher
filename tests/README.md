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
| `guard_neighbors` | Exhausting the neighbour list stops the run with a message rather than corrupting memory. |
| `guard_ngravprocs` | Asking for more gravity processes than there are GPUs is clamped, not fatal. |
| `periodic` | A fully periodic box conserves momentum and energy and has no density seam at its faces. |
| `periodic_block` | Block timesteps in a periodic box agree with shared ones and keep every particle in the box. |
| `guard_periodic` | A periodic box too small for the kernel stops the run with a message. |

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
