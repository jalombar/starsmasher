#!/bin/bash

# Define the location of your tools, ICs, and the directory where you
# have your simulations

TOOLS=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/tools
MESA=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/MESA/MESA_initial_3D_models/
WORK=$HOME/treball/StarSmasher_Collisions

# Initialize the base directory name in the work folder
cd $WORK || { echo "Failed to change directory to $WORK"; exit 1; }

# Echo the location of WORK
echo "Working directory: $WORK"

# Generate a random name with the prefix "Collision_"

while true; do
    suffix=$(openssl rand -hex 4)
    dir_name="Collision_$suffix"

    # Check if the directory already exists
    if [ ! -d "$dir_name" ]; then
        echo "Directory name generated: $dir_name"
        
        # Create the directory
        mkdir "$dir_name" || { echo "Failed to create directory $dir_name"; exit 1; }
        
        break
    fi
done

# Copy all necessary collision files from tools/ into that directory
cp "$TOOLS"/*_collision "$dir_name/" || { echo "Failed to copy collision files"; exit 1; }

# Change to the target directory
cd "$dir_name/" || { echo "Failed to change directory to $dir_name"; exit 1; }

# Rename the files
mv sph.init_collision sph.init || { echo "Failed to rename sph.init_collision to sph.init"; exit 1; }
mv sph.input_collision sph.input || { echo "Failed to rename sph.input_collision to sph.input"; exit 1; }

# Check if the 3D MESA SPH initial profiles exist, else exit.
# If they exist, copy them into the Collision directory

if ls $MESA/sph.start* 1> /dev/null 2>&1; then
    cp $MESA/sph.start* . || { echo "Failed to copy 3D MESA SPH initial profiles"; exit 1; }
    echo "Copied the following files from $MESA to $WORK/$dir_name:"
    ls $MESA/sph.start*
else
    echo "Error: 3D initial profiles do not exist in $MESA. Create them. Exiting now."
    exit 1
fi

# Echo the final working directory and simulation directory
echo "New simulation directory: $WORK/$dir_name"

# Get into the Collision directory and start the simulation
echo "Everything ready. Now call StSm_run.sh"

