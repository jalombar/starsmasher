LD = $(shell which mpif90)
LDFLAGS = 

MPIPATH = $(shell dirname $(shell dirname $(LD)))
MPITYPE = $(shell mpif90 --version | head -n1 | awk '{print $$1}')
ifeq ($(MPITYPE),GNU)
   # mpif90 was built using gfortran
   FC = mpif90 -ffixed-line-length-132
   # use the following two lines instead if mpif90 was not built with ifort but you want to use ifort
   #FC = ifort -132 -I$(MPIPATH)/include
   #LDFLAGS = -lifcore -lsvml -lifport -limf
else
   # mpif90 was (presumably) built using ifort
   FC = mpif90 -132
endif

CUDAPATH = $(shell dirname $(shell dirname $(shell which nvcc)))

OLEVEL = -O4 #-g 
FFLAGS = $(OLEVEL)
CFLAGS = $(OLEVEL)
CXXFLAGS = $(OLEVEL)
GRAVLIB = SPHgrav_lib2           # use this one for the direct force
#GRAVLIB = SPHgravtree_lib       # use this one for the  tree  force
LIBS =  -lm -lstdc++ #-L$(GRAVLIB) -lSPHgrav -L$(CUDAPATH)/lib64 -lcudart

FOBJS = kdtree2.o advance.o balAV3.o \
        getderivs.o getTemperature.o \
        initialize_hyperbolic.o init.o kernels.o \
        main.o output.o pressure.o\
        ran1.o spline.o splint.o\
        initialize_polyes.o initialize_polymces.o\
        initialize_corotating.o relax.o zeroin.o resplintmu.o grav.o initialize_triple.o temperaturefunction.o initialize_bps.o initialize_bpbh.o calccom.o elements.o initialize_smbh.o initialize_parent.o initialize_hyperbolic_binary_single.o skipahead.o compbest3.o changetf.o tstep.o initialize_grsph.o eatem.o useeostable.o hunt.o usekappatable.o initialize_asciiimage.o
#  grav.o 
CPUOBJS = $(FOBJS) cpu_grav.o
GPUOBJS = $(FOBJS) gpu_grav.o
%.o: %.f starsmasher.h
	$(FC) -c $(FFLAGS) $<
%.o :: %.f90 starsmasher.h
	$(FC) -c $(FFLAGS) $<

GPUEXEC = $(shell basename $(shell dirname $(shell pwd)))_gpu_sph
CPUEXEC = $(shell basename $(shell dirname $(shell pwd)))_cpu_sph

gpu: $(GPUOBJS)
	make -f Makefile -C $(GRAVLIB)
	$(LD) -o $(GPUEXEC) $(LDFLAGS) $(GPUOBJS) $(LIBS) -L$(GRAVLIB) -lSPHgrav -L$(CUDAPATH)/lib64 -lcudart -L/cm/shared/apps/openmpi/open64/64/1.10.1/lib64/ -lmpi_mpifh
	mv $(GPUEXEC) ..
	echo ***MADE VERSION THAT USES GPUS***

cpu: $(CPUOBJS)
	$(LD) -o $(CPUEXEC) $(LDFLAGS) $(CPUOBJS) $(LIBS) 
	mv $(CPUEXEC) ..
	echo ***MADE VERSION THAT DOES NOT NEED GPUS***

clean:  
	/bin/rm -rf *o *__*.f90 *_*.mod ../*_sph
	make clean -C $(GRAVLIB)
