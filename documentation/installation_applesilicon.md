# Installation on apple silicon

This guide describes how to build StarSmasher on an Apple Silicon Mac. On this platform the GPU build uses Metal instead of CUDA, and the main Makefile detects `Darwin/arm64` automatically. You should not need to edit source files or Makefiles before compiling.

These instructions apply to `parallel_bleeding_edge`, which is the actively developed version of the code.

## [0] Prerequisites

Install Apple's command line tools:

```bash
xcode-select --install
```

Install Homebrew if it is not already installed:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Install the required build tools:

```bash
brew install gcc open-mpi
```

Check that the MPI Fortran compiler is available:

```bash
mpif90 --version
```

No CUDA installation is required on Apple Silicon.

## [1] Obtaining starsmasher

Clone the repository:

```bash
git clone https://github.com/jalombar/starsmasher.git
cd starsmasher
```

## [2] Compiling starsmasher

Go to the source directory:

```bash
cd parallel_bleeding_edge/src
```

Check what the Makefile detects:

```bash
make config
```

On an Apple Silicon Mac, the output should include:

```text
OS               Darwin
architecture     arm64
GPU backend      metal
```

Build the GPU version:

```bash
make gpu
```

The executable is written to:

```text
parallel_bleeding_edge/parallel_bleeding_edge_gpu_sph
```

You can also build the CPU version:

```bash
make cpu
```

The CPU executable is written to:

```text
parallel_bleeding_edge/parallel_bleeding_edge_cpu_sph
```

### [2.1] Optional cpu-only builds

The CPU-only build is available for debugging, validation, and comparisons, but it is not the recommended way to run StarSmasher on Apple Silicon. Apple Silicon systems already include an integrated GPU in the processor, and the current M1 Pro benchmarks show the Metal GPU build running about 5x faster than the CPU-only build for the tested StarSmasher simulations. As a rough estimate for gravity-dominated runs, Pro-class Apple Silicon chips are expected to run about 5x faster with the Metal GPU build than with the CPU-only build on the same system, while Max-class chips are expected to run about 10x faster.

For that reason, the GPU build is recommended for normal Apple Silicon runs, and the examples below use `parallel_bleeding_edge_gpu_sph`.

## [3] Running starsmasher

Run StarSmasher from a simulation directory containing the usual input files. For example:

```bash
mpirun -np 2 /path/to/starsmasher/parallel_bleeding_edge/parallel_bleeding_edge_gpu_sph
```

Do not force `DYLD_SHARED_REGION=private` when running the Metal backend. On macOS, that setting can prevent Metal from seeing the GPU device.

## [4] Optional build settings

The Apple Silicon build uses smaller default array sizes than the Linux/NVIDIA build so it can compile and run on typical Mac laptops:

```text
NMAX=105000
NNMAX=96
NTAB=200000
```

You can override these values at compile time if your Mac has enough memory:

```bash
make gpu NMAX=200000 NNMAX=128 NTAB=400000
```

## [5] Installing splash

SPLASH can be used to visualize StarSmasher output files:

```bash
brew install splash
```

Open StarSmasher snapshots with:

```bash
splash -f starsmasher out*.sph
```

If Homebrew upgrades HDF5 and an existing SPLASH binary can no longer find the expected HDF5 library, reinstall SPLASH:

```bash
brew reinstall splash
```

## [6] Troubleshooting

If `mpif90` is not found, make sure Homebrew is on your `PATH`. On Apple Silicon this usually means adding Homebrew's shell environment to your shell startup file:

```bash
eval "$(/opt/homebrew/bin/brew shellenv)"
```

If `make gpu` does not select Metal on an Apple Silicon Mac, check the detected platform:

```bash
uname -s
uname -m
make config
```

You can force the Metal backend explicitly:

```bash
make gpu GPU_BACKEND=metal
```
