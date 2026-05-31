# Build StarSmasher from this directory.

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)
APPLE_SILICON := $(filter Darwin-arm64,$(UNAME_S)-$(UNAME_M))

GPU_BACKEND ?= auto
ifeq ($(GPU_BACKEND),auto)
  ifneq ($(APPLE_SILICON),)
    GPU_BACKEND := metal
  else ifneq ($(shell command -v nvcc 2>/dev/null),)
    GPU_BACKEND := cuda
  else
    GPU_BACKEND := none
  endif
endif

MPIF90 ?= mpif90
LD = $(MPIF90)
LDFLAGS =

MPITYPE = $(shell $(MPIF90) --version 2>/dev/null | head -n1 | awk '{print $$1}')
ifeq ($(MPITYPE),GNU)
  FC_FIXED_FLAGS = -ffixed-line-length-132 -fallow-argument-mismatch
else
  FC_FIXED_FLAGS = -132
endif

ifeq ($(UNAME_S),Darwin)
  SDKROOT ?= $(shell xcrun --show-sdk-path 2>/dev/null)
  MACOSX_DEPLOYMENT_TARGET ?= $(shell sw_vers -productVersion | awk -F. '{print $$1"."$$2}')
  BUILD_ENV = SDKROOT=$(SDKROOT) MACOSX_DEPLOYMENT_TARGET=$(MACOSX_DEPLOYMENT_TARGET)
  FC = env $(BUILD_ENV) $(MPIF90) $(FC_FIXED_FLAGS)
  LINK = env $(BUILD_ENV) $(LD)
  PLATFORM_FFLAGS = -frecursive
else
  FC = $(MPIF90) $(FC_FIXED_FLAGS)
  LINK = $(LD)
  PLATFORM_FFLAGS =
endif

OLEVEL ?= -O4
FFLAGS = $(OLEVEL) $(PLATFORM_FFLAGS)
CFLAGS = $(OLEVEL)
CXXFLAGS = $(OLEVEL)
LIBS = -lm

BINDIR ?= ../bin
CPUEXEC ?= starsmasher_cpu
GPUEXEC ?= starsmasher_gpu
LIMITS_FILE = starsmasher_limits.generated.h

ifneq ($(APPLE_SILICON),)
  DEFAULT_NMAX = 105000
  DEFAULT_NNMAX = 96
  DEFAULT_NTAB = 200000
else
  DEFAULT_NMAX = 900000
  DEFAULT_NNMAX = 170
  DEFAULT_NTAB = 10000000
endif

NMAX ?= $(DEFAULT_NMAX)
NNMAX ?= $(DEFAULT_NNMAX)
NTAB ?= $(DEFAULT_NTAB)

ifeq ($(GPU_BACKEND),metal)
  GRAVLIB = SPHgrav_lib_metal
  GPU_GRAV_OBJ = gpu_grav_metal.o
  GPU_LIBS = -L$(GRAVLIB) -lSPHgrav -framework Foundation -framework Metal -lobjc -lc++
else ifeq ($(GPU_BACKEND),cuda)
  GRAVLIB = SPHgrav_lib2
  GPU_GRAV_OBJ = gpu_grav.o
  CUDAPATH ?= $(shell dirname $$(dirname $$(command -v nvcc 2>/dev/null)))
  CUDA_EXTRA_LIBS ?=
  GPU_LIBS = -L$(GRAVLIB) -lSPHgrav -L$(CUDAPATH)/lib64 -lcudart -lstdc++ $(CUDA_EXTRA_LIBS)
else ifeq ($(GPU_BACKEND),none)
  GRAVLIB =
  GPU_GRAV_OBJ =
  GPU_LIBS =
else
  $(error Unknown GPU_BACKEND "$(GPU_BACKEND)"; use auto, metal, cuda, or none)
endif

FOBJS = eos_table_data.o kdtree2.o advance.o balAV3.o              \
        getderivs.o getTemperature.o                              \
        initialize_hyperbolic.o init.o kernels.o                  \
        main.o output.o pressure.o                                \
        ran1.o spline.o splint.o initialize_polyes.o              \
        initialize_polymces.o initialize_corotating.o             \
        relax.o zeroin.o resplintmu.o grav.o initialize_triple.o  \
        temperaturefunction.o initialize_bps.o                    \
        initialize_bpbh.o calccom.o elements.o                    \
        initialize_smbh.o initialize_parent.o                     \
        initialize_multiequalmass.o initialize_rescale.o          \
        initialize_hyperbolic_binary_single.o skipahead.o         \
        compbest3.o changetf.o tstep.o initialize_grsph.o         \
        eatem.o useeostable.o hunt.o usekappatable.o              \
        initialize_asciiimage.o

CPUOBJS = $(FOBJS) cpu_grav.o
GPUOBJS = $(FOBJS) $(GPU_GRAV_OBJ)

.PHONY: all info cpu gpu clean gravlib require_gpu_backend

all: cpu

info:
	@echo "os=$(UNAME_S)"
	@echo "arch=$(UNAME_M)"
	@echo "gpu_backend=$(GPU_BACKEND)"
	@echo "mpif90=$(MPIF90)"
	@echo "nmax=$(NMAX)"
	@echo "nnmax=$(NNMAX)"
	@echo "ntab=$(NTAB)"

$(LIMITS_FILE): Makefile
	@tmp="$@.tmp"; \
	{ \
	  echo "      integer nmax,nnmax"; \
	  echo "      parameter(nmax=$(NMAX),nnmax=$(NNMAX))"; \
	  echo "      integer ntab"; \
	  echo "      parameter(ntab=$(NTAB))"; \
	} > "$$tmp"; \
	if ! cmp -s "$$tmp" "$@"; then mv "$$tmp" "$@"; else rm -f "$$tmp"; fi

%.o: %.f starsmasher.h $(LIMITS_FILE)
	$(FC) -c $(FFLAGS) $<

%.o: %.f90 starsmasher.h $(LIMITS_FILE)
	$(FC) -c $(FFLAGS) $<

useeostable.o initialize_polyes.o: eos_table_data.o

$(BINDIR):
	mkdir -p $@

cpu: $(BINDIR)/$(CPUEXEC)

$(BINDIR)/$(CPUEXEC): $(CPUOBJS) | $(BINDIR)
	$(LINK) -o $@ $(LDFLAGS) $(CPUOBJS) $(LIBS)
	echo ***MADE VERSION THAT DOES NOT NEED GPUS***

require_gpu_backend:
	@if [ "$(GPU_BACKEND)" = "none" ]; then \
	  echo "No GPU backend available. Use GPU_BACKEND=metal on Apple Silicon or GPU_BACKEND=cuda with CUDA."; \
	  exit 1; \
	fi

gravlib: require_gpu_backend
	$(MAKE) -C $(GRAVLIB) CUDAPATH="$(CUDAPATH)"

gpu: require_gpu_backend gravlib $(BINDIR)/$(GPUEXEC)

$(BINDIR)/$(GPUEXEC): $(GPUOBJS) | $(BINDIR)
	$(LINK) -o $@ $(LDFLAGS) $(GPUOBJS) $(LIBS) $(GPU_LIBS)
	echo ***MADE VERSION THAT USES $(GPU_BACKEND) GPUS***

clean:
	/bin/rm -f *.o *.mod *__*.f90 *_*.mod $(LIMITS_FILE)
	/bin/rm -f $(BINDIR)/$(CPUEXEC) $(BINDIR)/$(GPUEXEC)
	@if [ -d SPHgrav_lib2 ]; then $(MAKE) --no-print-directory -C SPHgrav_lib2 clean; fi
	@if [ -d SPHgrav_lib_metal ]; then $(MAKE) --no-print-directory -C SPHgrav_lib_metal clean; fi
