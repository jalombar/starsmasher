# Note: This Makefile will output the executable to a "bin"
#       directory one level up. Make sure it exists.
#       Pau Amaro Seoane, 24 August 2024

# Choose your compiler: gfortran or ifort       
# gfortran
FC      = mpif90 -ffixed-line-length-132 -fallow-argument-mismatch
LDFLAGS = 

# ifort 
#FC = mpif90 -132

LD = $(shell which mpif90)

MPIPATH  = $(shell dirname $(shell dirname $(LD)))
CUDAPATH = /usr/local/cuda-12.6/

OLEVEL = -O4 #-g 
FFLAGS = $(OLEVEL)
CFLAGS = $(OLEVEL)
CXXFLAGS = $(OLEVEL)
GRAVLIB = SPHgrav_lib2 # use this one for the direct force
LIBS =  -lm -lstdc++   #-L$(GRAVLIB) -lSPHgrav -L$(CUDAPATH)/lib64 -lcudart

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
#  grav.o 
CPUOBJS = $(FOBJS) cpu_grav.o
GPUOBJS = $(FOBJS) gpu_grav.o
%.o: %.f starsmasher.h
	$(FC) -c $(FFLAGS) $<
%.o :: %.f90 starsmasher.h
	$(FC) -c $(FFLAGS) $<
#%.o: %.c starsmasher.h
#	$(CC) -c $(CFLAGS) $<

GPUEXEC = test_gpu_sph_REMOVE
CPUEXEC = test_cpu_sph_REMOVE

gpu: $(GPUOBJS)
	make -f Makefile -C $(GRAVLIB)
	$(LD) -o $(GPUEXEC) $(LDFLAGS) $(GPUOBJS) $(LIBS)           \
	   	-L$(GRAVLIB) -lSPHgrav -L$(CUDAPATH)/lib64 -lcudart \
	       	-L/cm/shared/apps/openmpi/open64/64/1.10.1/lib64/   \
		-lmpi_mpifh
	mv $(GPUEXEC) ../bin/
	echo ***MADE VERSION THAT USES GPUS***

cpu: $(CPUOBJS)
	$(LD) -o $(CPUEXEC) $(LDFLAGS) $(CPUOBJS) $(LIBS) 
	mv $(CPUEXEC) ../bin/
	echo ***MADE VERSION THAT DOES NOT NEED GPUS***

clean:  
	/bin/rm -rf *o *__*.f90 *_*.mod ../*_sph
	make clean -C $(GRAVLIB)
