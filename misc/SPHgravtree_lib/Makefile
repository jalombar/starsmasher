.SUFFIXES: .cu 

CXX  := g++
CC   := gcc
AR   := ar
RANLIB := ranlib

CUDAPATH       := /usr/local/cuda
CUDAINCLUDE    := -I$(CUDAPATH)/include
NVCC           := $(CUDAPATH)/bin/nvcc

#use this for Fermi cards
#NVCCFLAGS := -arch=sm_20
#use this for GTX280/290/295, or any GT200 architecture with DP support
#NVCCFLAGS := -arch=sm_13
NVCCFLAGS += -O4 -g  $(CUDAINCLUDE)  -I./ -Xptxas -v,-abi=no 
# NVCCFLAGS += -maxrregcount=32
NVCCFLAGS += -Xcompiler="-Wall"
CUDA_LIBS = -L$(CUDAPATH)/lib64 -lcudart

CXXFLAGS  +=  -Wall -g -O4
# LDFLAGS   +=  -fopenmp
LDGPUGLAGS := $(LDFLAGS) $(CUDA_LIBS)

OBJS = grav_force_tree.cu_o
GRAVLIB = libSPHgrav.a

all: $(GRAVLIB)

$(GRAVLIB): $(OBJS)
	make -C ./libSequoia
	/bin/rm -f $@
	mkdir libSequoia-tmp
	cd libSequoia-tmp; $(AR) x ../libSequoia/libsequoia.a
	$(AR) ruv $@ $(OBJS) ./libSequoia-tmp/*
	/bin/rm -rf libSequoia-tmp

.cpp.o: 
	$(CXX) $(CXXFLAGS) -c $< -o $@

%.cu_o:  %.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.o *.cu_o $(GRAVLIB)

$(GRAVLIB):  libSequoia/include/sequoiaInterface.h

