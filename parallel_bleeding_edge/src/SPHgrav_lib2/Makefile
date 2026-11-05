.SUFFIXES: .cu 

CXX  := g++
CC   := gcc
AR   := ar ruv
RANLIB := ranlib

# Path to your CUDA installation: CUDAPATH should contain the directory lib64 that contains libcudart.so as well as the bin directory that contains nvcc
CUDAPATH       := $(shell dirname $(shell dirname $(shell which nvcc)))
CUDAINCLUDE    := -I$(CUDAPATH)/include
NVCC           := $(CUDAPATH)/bin/nvcc

# Automatically determine the compute capability of the GPU, in the undotted form
# nvcc wants (12.0 -> 120).  Ask nvidia-smi first, then the query tool that ships
# with CUDA, which reads the driver directly and still answers when nvidia-smi is
# too old to know --query-gpu=compute_cap.  Keep only all-digit answers, so that a
# warning printed on stdout ("Failed to initialize NVML: ...") cannot end up
# inside the -arch flag.  With more than one GPU, take the lowest: code built for
# a lower capability still runs on a higher card, because the driver recompiles
# the embedded PTX at start-up, but not the other way round.
COMPUTE_CAPABILITY := $(shell { nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | tr -d '.' ; \
                                $(CUDAPATH)/bin/__nvcc_device_query 2>/dev/null ; echo ; } \
                              | grep -x '[0-9][0-9]*' | sort -n | head -1)

# Check if the capability is empty (no GPU detected) and set default value to 90 if so
ifneq ($(COMPUTE_CAPABILITY),)
    NVCCFLAGS := -arch=sm_$(COMPUTE_CAPABILITY)
else
    NVCCFLAGS := -arch=sm_90  # 80 for a100, 90 for h100
endif

NVCCFLAGS += -O4 -g  $(CUDAINCLUDE)  -I./ -Xptxas -v ##########,-abi=no 
NVCCFLAGS += -Xcompiler="-Wall"
CUDA_LIBS = -L$(CUDAPATH)/lib64 -lcudart

CXXFLAGS  +=  -Wall -g -O4
# the host code of the library uses OpenMP when StarSmasher is built with OPENMP=1 (the top Makefile passes OMPCXX)
ifneq ($(OMPCXX),)
override NVCCFLAGS += -Xcompiler=$(OMPCXX)
override CXXFLAGS += $(OMPCXX)
endif
# LDFLAGS   +=  -fopenmp
LDGPUGLAGS := $(LDFLAGS) $(CUDA_LIBS)

OBJS = grav_force_direct.cu_o grav_tree_gpu.cu_o
GRAVLIB = libSPHgrav.a

all: $(GRAVLIB)

# the flags of the last build, so that a change of OMPCXX (the top Makefile passes it) or of the optimisation rebuilds the objects
FLAGSTAMP = cu_flags.stamp
$(FLAGSTAMP): FORCE
	@echo "$(subst ",,$(NVCCFLAGS) $(CXXFLAGS))" > $@.tmp; if cmp -s $@.tmp $@; then /bin/rm -f $@.tmp; else mv $@.tmp $@; fi
FORCE:

$(GRAVLIB): $(OBJS)
	    /bin/rm -f $@
	    $(AR) $@ $(OBJS)

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.cu_o:  %.cu cutil.h cuVec3.h cuVector.h $(FLAGSTAMP)
	 $(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.cu_o $(GRAVLIB) $(FLAGSTAMP)
