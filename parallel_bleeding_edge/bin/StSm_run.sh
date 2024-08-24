#!/usr/bin/env zsh

# Location of the executable
BIN=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/bin/test_gpu_sph

# Log the start time and command
echo "Starting MPI run at $(date)" | tee -a mpi_run.log
echo "Running command: mpirun -np 8 $BIN" | tee -a mpi_run.log

# Run the MPI program
mpirun -np 8 $BIN > output_$(date +%Hh%Mm%Ss_%m_%d_%Y).txt 2> stderr_$(date +%Hh%Mm%Ss_%m_%d_%Y).txt &

# Get the PID of the background process
pid=$!

# Wait for the process to finish
wait $pid

# Log the completion of the job
if [ $? -eq 0 ]; then
    echo "MPI run completed successfully at $(date)" | tee -a mpi_run.log
else
    echo "MPI run failed at $(date)" | tee -a mpi_run.log
fi
