#!/bin/sh

mpirun -np 8 ../bin/*_sph > output.txt 2> stderr.txt &
