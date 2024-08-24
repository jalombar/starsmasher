#!/bin/bash

# Define the location of your tools and ICs
TOOLS=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/tools
MESA=/projects/StarSmasher/starsmasher/parallel_bleeding_edge/MESA_ICs

# Initialize the base directory name
dir_name="Collision"
suffix=""

# Check if the directory already exists and increment the suffix if it does
if [ -d "$dir_name" ]; then
    i=2
    while [ -d "${dir_name}_$i" ]; do
        i=$((i + 1))
    done
    suffix="_$i"
    dir_name="${dir_name}${suffix}"
fi

# Step 1: Create the directory with the final name
#         and echo its name.
mkdir "$dir_name"
echo "Created directory: $dir_name"

# Step 2: Copy all necessary collision files from tools/ into that directory
cp $TOOLS/* "$dir_name/"

# Step 3: Move the 3D MESA SPH initial profiles into the Collision directory
mv $MESA/sph.start* "$dir_name/"

# Step 4: Get into the Collision directory and start the simulation
cd "$dir_name/"
echo "Run StSm_run.sh"
