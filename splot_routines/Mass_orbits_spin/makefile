# Makefile for splot
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
FFLAGS= -O3 #-w90 -w95 -w
#FFLAGS= -O2 #-W1 -ip -ftz #-axP
LIBS      = #-L../. -lGPUsph_gs -L/usr/local/cuda/lib64 -lcudart #-lcuda #-L/opt/MDGRAPE3/lib -lmdgrape3

OBJECTS = pplot.o orig.o compbest3.o calccom.o elements.o robusttrajectories.o

EXEC	= $(shell basename $(shell pwd))

$(EXEC) : $(OBJECTS)
	$(LD) -o $@ $(LDFLAGS) $(LIBS) $(OBJECTS) $(LIBS) $(FFLAGS)
	cp $(EXEC) ~/bin

.SUFFIXES:
.SUFFIXES: .o .c .f

%.o: %.f spha.h
	$(FC) -c $(FFLAGS) $<

clean:
	/bin/rm -rf *o *__*.f90 *_*.mod
