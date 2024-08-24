#!/bin/sh

BIN=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/bin/test_gpu_sph

mpirun -np 8 $BIN > output.txt 2> stderr.txt &
