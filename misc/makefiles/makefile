FC = $(shell which ifort)
LD = $(shell which mpif90)
MPIPATH = $(shell dirname $(shell dirname $(LD)))
FCFLAGS = -132 -I$(MPIPATH)/include #-mcmodel=large
LDFLAGS = -lifcore -lsvml -lifport -limf

ifeq ($(FC),)
  FC = $(LD)
    ifeq ($(FC),)
      $(error mpif90 not found)
    endif
  FCFLAGS = -ffixed-line-length-132
  LDFLAGS = -lgfortran
#else
endif
ifeq ($(LD),)
   $(error MPI not found)
endif

CUDAPATH = $(shell dirname $(shell dirname $(shell which nvcc)))
#CUDAPATH = /sw/kfs/cuda/5.0/linux_binary

OLEVEL = -O4 #-g 
FFLAGS = $(OLEVEL)
CFLAGS = $(OLEVEL)
CXXFLAGS = $(OLEVEL)
GRAVLIB = SPHgrav_lib           # use this one for the direct force
#GRAVLIB = SPHgravtree_lib       # use this one for the  tree  force
LIBS =  -lm #-L$(GRAVLIB) -lSPHgrav -L$(CUDAPATH)/lib64 -lcudart

FOBJS = kdtree2.o advance.o balAV3.o \
        getderivs.o getTemperature.o \
        initialize_hyperbolic.o init.o kernels.o \
        main.o output.o pressure.o\
        ran1.o spline.o splint.o\
        initialize_polyes.o initialize_polymces.o\
        initialize_corotating.o relax.o zeroin.o resplintmu.o grav.o initialize_triple.o temperaturefunction.o initialize_bps.o initialize_bpbh.o calccom.o elements.o initialize_smbh.o initialize_egrgparent.o initialize_hyperbolic_binary_single.o skipahead.o compbest3.o componentize3.o changetf.o tstep.o initialize_grsph.o eatem.o useeostable.o hunt.o usekappatable.o initialize_asciiimage.o
#  grav.o 
COBJS = rtc.o #peano_sort.o
CPUOBJS = $(FOBJS) $(COBJS) cpu_grav.o
GPUOBJS = $(FOBJS) $(COBJS) gpu_grav.o
%.o: %.f starsmasher.h
	$(FC) $(FCFLAGS) -c $(FFLAGS) $<
%.o :: %.f90 starsmasher.h
	$(FC) $(FCFLAGS) -c $(FFLAGS) $<
#%.o: %.c starsmasher.h
#	$(CC) -c $(CFLAGS) $<

GPUEXEC = $(shell basename $(shell dirname $(shell pwd)))_gpu_sph
CPUEXEC = $(shell basename $(shell dirname $(shell pwd)))_cpu_sph

gpu: $(GPUOBJS)
	make -C $(GRAVLIB)
	$(LD) -o $(GPUEXEC) $(LDFLAGS) $(GPUOBJS) $(LIBS) -L$(GRAVLIB) -lSPHgrav -L$(CUDAPATH)/lib64 -lcudart
	mv $(GPUEXEC) ..
	echo ***MADE VERSION THAT USES GPUS***

cpu: $(CPUOBJS)
	$(LD) -o $(CPUEXEC) $(LDFLAGS) $(CPUOBJS) $(LIBS) 
	mv $(CPUEXEC) ..
	echo ***MADE VERSION THAT DOES NOT NEED GPUS***

clean:  
	/bin/rm -rf *o *__*.f90 *_*.mod ../*_sph
	make clean -C $(GRAVLIB)
