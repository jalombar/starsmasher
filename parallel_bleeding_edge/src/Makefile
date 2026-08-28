# Makefile for StarSmasher (parallel_bleeding_edge)
#
# Everything below is detected from your environment. Nothing is hardcoded to a
# particular machine, CUDA version, or MPI installation. To see what was
# detected before building, run:
#
#     make config
#
# Any detected value can be overridden on the command line, e.g.
#     make CUDAPATH=/opt/cuda OLEVEL=-O2
#     make cpu            # build without GPU support

# ---------------------------------------------------------------- MPI --------
# mpif90 is used as both compiler and linker, so it supplies its own include
# paths and MPI libraries. Do not add -lmpi_mpifh or similar by hand: the right
# set differs between OpenMPI and MPICH, and mpif90 already knows it.
MPIF90 ?= mpif90
LD      = $(MPIF90)
LDFLAGS =

MPIF90_PATH := $(shell command -v $(MPIF90) 2>/dev/null)
ifeq ($(MPIF90_PATH),)
$(error Cannot find '$(MPIF90)' in PATH. Install an MPI implementation, or set MPIF90=/full/path/to/mpif90)
endif

# Which compiler is mpif90 wrapping? Determines the fixed-form flags we need.
MPI_VENDOR := $(shell $(MPIF90) --version 2>/dev/null | head -n1 | awk '{print $$1}')

ifeq ($(MPI_VENDOR),GNU)
  # gfortran. -fallow-argument-mismatch is required from gfortran 10 onward
  # (the MPI calls in this code pass different types to the same argument) but
  # does not exist before it, so probe rather than assume.
  ALLOW_MISMATCH := $(shell echo "PROGRAM T; END PROGRAM T" | \
      $(MPIF90) -fallow-argument-mismatch -x f95 -o /dev/null - 2>/dev/null \
      && echo -fallow-argument-mismatch)
  FORM_FLAGS = -ffixed-line-length-132 $(ALLOW_MISMATCH)
else
  # ifort/ifx and most others spell the fixed-form line length limit as -132.
  FORM_FLAGS = -132
endif

FC = $(MPIF90) $(FORM_FLAGS)

# --------------------------------------------------------------- CUDA --------
# Derive the CUDA root from nvcc's location rather than naming a version.
NVCC_PATH := $(shell command -v nvcc 2>/dev/null)
ifneq ($(NVCC_PATH),)
  CUDAPATH ?= $(shell dirname $(shell dirname $(NVCC_PATH)))
endif

# The runtime library lives in lib64 on most installs, but CUDA 12+ may only
# provide targets/<arch>/lib. Take the first directory that actually exists.
ifneq ($(CUDAPATH),)
  CUDA_LIBDIR := $(firstword $(wildcard $(CUDAPATH)/lib64 \
                                        $(CUDAPATH)/targets/*/lib \
                                        $(CUDAPATH)/lib))
endif

# ------------------------------------------------------------- code model ----
# The static arrays in starsmasher.h are large: at the default nmax the block
# sits close to the 2 GB limit that -mcmodel=small imposes on x86-64, and
# raising nmax, nnmax, ntab or kdm pushes it over, giving a "relocation
# truncated" or "additional relocation overflows" error at link time that does
# not say what actually went wrong. Use the medium model where it is supported;
# it is not a valid option on every architecture, so probe for it.
ifeq ($(shell uname -m),x86_64)
  MCMODEL ?= $(shell echo "PROGRAM T; END PROGRAM T" | \
      $(MPIF90) -mcmodel=medium -x f95 -o /dev/null - >/dev/null 2>&1 \
      && echo -mcmodel=medium)
else
  MCMODEL ?=
endif

# ------------------------------------------------------------------ flags ----
OLEVEL   ?= -O4
FFLAGS    = $(OLEVEL) $(MCMODEL)
CFLAGS    = $(OLEVEL)
CXXFLAGS  = $(OLEVEL)

GRAVLIB = SPHgrav_lib2
LIBS    = -lm -lstdc++

# ---------------------------------------------------------------- objects ----
FOBJS = kdtree2.o advance.o balAV3.o                             \
        getderivs.o getTemperature.o                             \
        initialize_hyperbolic.o init.o kernels.o                 \
        main.o output.o pressure.o                               \
        ran1.o spline.o splint.o initialize_polyes.o             \
        initialize_polymces.o initialize_corotating.o            \
        relax.o zeroin.o resplintmu.o grav.o initialize_triple.o \
        temperaturefunction.o initialize_bps.o                   \
        initialize_bpbh.o calccom.o elements.o                   \
        initialize_smbh.o initialize_parent.o                    \
        initialize_multiequalmass.o initialize_rescale.o         \
        initialize_hyperbolic_binary_single.o skipahead.o        \
        compbest3.o changetf.o tstep.o initialize_grsph.o        \
        eatem.o useeostable.o hunt.o usekappatable.o             \
        initialize_asciiimage.o

CPUOBJS = $(FOBJS) cpu_grav.o
GPUOBJS = $(FOBJS) gpu_grav.o

# The executable is named after the directory holding src, so that a build in
# my_run/src produces my_run_gpu_sph, and is placed one level above src.
GPUEXEC = $(shell basename $(shell dirname $(shell pwd)))_gpu_sph
CPUEXEC = $(shell basename $(shell dirname $(shell pwd)))_cpu_sph

%.o: %.f starsmasher.h
	$(FC) -c $(FFLAGS) $<

%.o :: %.f90 starsmasher.h
	$(FC) -c $(FFLAGS) $<

# ----------------------------------------------------------------- targets ---
gpu: check-cuda $(GPUOBJS)
	$(MAKE) -C $(GRAVLIB) CUDAPATH=$(CUDAPATH)
	$(LD) -o $(GPUEXEC) $(LDFLAGS) $(GPUOBJS) $(LIBS) \
		-L$(GRAVLIB) -lSPHgrav -L$(CUDA_LIBDIR) -lcudart
	cp -p $(GPUEXEC) ..
	@echo "***MADE VERSION THAT USES GPUS***  ->  ../$(GPUEXEC)"

cpu: $(CPUOBJS)
	$(LD) -o $(CPUEXEC) $(LDFLAGS) $(CPUOBJS) $(LIBS)
	cp -p $(CPUEXEC) ..
	@echo "***MADE VERSION THAT DOES NOT NEED GPUS***  ->  ../$(CPUEXEC)"

check-cuda:
ifeq ($(NVCC_PATH),)
	@echo "ERROR: nvcc not found in PATH, so the GPU version cannot be built."; \
	echo "  - If CUDA is installed, add its bin directory to PATH, or pass"; \
	echo "    CUDAPATH=/path/to/cuda explicitly."; \
	echo "  - If you have no NVIDIA GPU, build the CPU version instead:  make cpu"; \
	exit 1
else ifeq ($(CUDA_LIBDIR),)
	@echo "ERROR: found nvcc at $(NVCC_PATH), but no CUDA runtime library"; \
	echo "  directory under $(CUDAPATH) (looked for lib64, targets/*/lib, lib)."; \
	echo "  Pass the correct root explicitly, e.g. make CUDAPATH=/opt/cuda"; \
	exit 1
else
	@:
endif

config:
	@echo "MPI"
	@echo "  mpif90          $(MPIF90_PATH)"
	@echo "  wrapping        $(MPI_VENDOR)"
	@echo "  fixed-form flag $(FORM_FLAGS)"
	@echo "CUDA"
	@echo "  nvcc            $(if $(NVCC_PATH),$(NVCC_PATH),NOT FOUND - only 'make cpu' will work)"
	@echo "  CUDAPATH        $(CUDAPATH)"
	@echo "  runtime libdir  $(CUDA_LIBDIR)"
	@echo "Compilation"
	@echo "  FFLAGS          $(FFLAGS)"
	@echo "  code model      $(if $(MCMODEL),$(MCMODEL),default (small))"
	@echo "Output"
	@echo "  GPU executable  ../$(GPUEXEC)"
	@echo "  CPU executable  ../$(CPUEXEC)"

clean:
	/bin/rm -rf *.o *__*.f90 *.mod $(GPUEXEC) $(CPUEXEC) ../$(GPUEXEC) ../$(CPUEXEC)
	$(MAKE) clean -C $(GRAVLIB)

.PHONY: gpu cpu clean config check-cuda
