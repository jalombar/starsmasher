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
| `guard_neighbors` | Exhausting the neighbour list stops the run with a message rather than corrupting memory. |
| `guard_ngravprocs` | Asking for more gravity processes than there are GPUs is clamped, not fatal. |

The last four exist because each corresponds to a real defect: results that
depended on rank count, a CPU/GPU disagreement, an unbounded write past the end
of an array, and a segfault when `ngravprocs` exceeded the GPU count.
