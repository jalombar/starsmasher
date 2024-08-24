# Take into account the localtion of CUDAPATH and the
# "-abi=no" flag in NVCCFLAGS
# Pau Amaro Seoane, Berlin 24 August 2024

.SUFFIXES: .cu 

CXX  := g++
CC   := gcc
AR   := ar ruv
RANLIB := ranlib

# Path to your CUDA installation
# Make sure you do not have different versions installed.
# If that's the case, use the same CUDAPATH here and in the
# main folder, and choose the most recent version for both.

CUDAPATH       := /usr/local/cuda-12.6
CUDAINCLUDE    := -I$(CUDAPATH)/include
NVCC           := $(CUDAPATH)/bin/nvcc

# Remove the "-abi=no" flag to avoid bug
#
NVCCFLAGS := -arch=sm_75
#NVCCFLAGS += -O4 -g  $(CUDAINCLUDE)  -I./ -Xptxas -v,-abi=no 
NVCCFLAGS += -O4 -g  $(CUDAINCLUDE)  -I./ -Xptxas -v

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

