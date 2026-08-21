.SUFFIXES: .cu 

CXX  := g++
CC   := gcc
AR   := ar ruv
RANLIB := ranlib

# Path to your CUDA installation: CUDAPATH should contain the directory lib64 that contains libcudart.so as well as the bin directory that contains nvcc
CUDAPATH       := $(shell dirname $(shell dirname $(shell which nvcc)))
CUDAINCLUDE    := -I$(CUDAPATH)/include
NVCC           := $(CUDAPATH)/bin/nvcc

# Automatically determine the compute capability of the GPU
COMPUTE_CAPABILITY := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '.')


# Check if the capability is empty (no GPU detected) and set default value to 90 if so
ifneq ($(COMPUTE_CAPABILITY),)
    NVCCFLAGS := -arch=sm_$(COMPUTE_CAPABILITY)
else
    NVCCFLAGS := -arch=sm_89  # 5.2 for carrv122, 2.0 for gonzales, 6.1 for physics, 8.0 for a100, 9.0 for h100
endif

NVCCFLAGS += -O4 -g  $(CUDAINCLUDE)  -I./ -Xptxas -v ##########,-abi=no 
NVCCFLAGS += -Xcompiler="-Wall"
CUDA_LIBS = -L$(CUDAPATH)/lib64 -lcudart

CXXFLAGS  +=  -Wall -g -O4
# LDFLAGS   +=  -fopenmp
LDGPUGLAGS := $(LDFLAGS) $(CUDA_LIBS)

OBJS = grav_force_direct.cu_o
GRAVLIB = libSPHgrav.a

all: $(GRAVLIB)

$(GRAVLIB): $(OBJS)
	    /bin/rm -f $@
	    $(AR) $@ $(OBJS)

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.cu_o:  %.cu
	 $(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.cu_o $(GRAVLIB)

$(GRAVLIB): cutil.h cuVec3.h cuVector.h
