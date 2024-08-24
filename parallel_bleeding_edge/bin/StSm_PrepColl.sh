#!/bin/bash

# Define the location of your tools, ICs and the directory where you
# have your simulations

TOOLS=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/tools
MESA=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/MESA/MESA_initial_3D_models/
WORK=$HOME/treball/StarSmasher_Collisions

# Initialize the base directory name in the work folder
cd $WORK

# Generate a random name with the prefix "Collision_"

while true; do
    suffix=$(openssl rand -hex 4)
    dir_name="Collision_$suffix"

    # Check if the directory already exists
    if [ ! -d "$dir_name" ]; then
        echo "Directory name generated: $dir_name"
        
        # Create the directory
        mkdir "$dir_name"
        
        break
    fi
done

# Copy all necessary collision files from tools/ into that directory

cp $TOOLS/*_collision "$dir_name/"
cd "$dir_name/"
mv sph.init_collision  sph.init
mv sph.input_collision sph.input

# Check if the 3D MESA SPH initial profiles exist, else exit.
# If they exist, move them into the Collision directory

if ls $MESA/sph.start* 1> /dev/null 2>&1; then
    mv $MESA/sph.start* "$dir_name/"
else
    echo "Error: 3D initial profiles do not exist in $MESA. Create them. Exiting now."
    exit 1
fi

# Get into the Collision directory and start the simulation

cd "$dir_name/"
echo "Everything ready. Now call StSm_run.sh"
