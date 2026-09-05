CXX := clang++
AR := ar ruv
RANLIB := ranlib

OLEVEL ?= -O4
CXXFLAGS := $(OLEVEL) -Wall -std=c++17 -fobjc-arc

OBJS := grav_force_direct_metal.o
GRAVLIB := libSPHgrav.a

all: $(GRAVLIB)

$(GRAVLIB): $(OBJS)
	/bin/rm -f $@
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

%.o: %.mm
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	/bin/rm -f $(OBJS) $(GRAVLIB)
