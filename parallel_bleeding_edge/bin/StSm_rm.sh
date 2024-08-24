#!/bin/bash

# Sometimes we do not need to store all
# snapshots. In that case, you can remove
# intervals of them with this script.
# Use it like:
# $ snapshrm_StSm.sh 250 3022

# Set the range for removal
# as arguments $1 and $2
start=$1
end=$2

# Loop through each file and remove it if it falls within the interval
for file in out*.sph; do
    # Extract the number from the filename
    number=$(echo $file | sed -E 's/out([0-9]+)\.sph/\1/')

    # Check if the number is within the specified range
    if [ "$number" -ge "$start" ] && [ "$number" -le "$end" ]; then
        echo "Removing $file"
        rm "$file"
    fi
done
