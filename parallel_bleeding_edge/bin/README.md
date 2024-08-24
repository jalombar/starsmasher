# StarSmasher Collision Simulation Scripts

This repository contains several shell scripts designed to facilitate the management and execution of collision simulations using the StarSmasher software. Below is an overview of each script and its functionality. Note that I am assuming you use the `zsh` shell. Otherwise, please change the kernel line in all scripts.

## Scripts Overview

### 1. `StSm_clean.sh`

**Purpose:**  
This script is used to clean up the working directory by removing all `.sph` files, as well as `output.txt` and `stderr.txt`. These files are typically generated during simulation runs, and this script helps maintain a tidy workspace by deleting them.

**Usage:**  
Simply run the script to remove the files:
```
./StSm_clean.sh
```

### 2. `StSm_plot.sh`

**Purpose:**  
This script checks whether the `splash` visualization tool is installed on your system. If `splash` is not found, it provides installation instructions and exits. If `splash` is installed, the script runs it with the `starsmasher` files passed as arguments.

**Usage:**  
Execute the script and provide any necessary arguments for `splash`:
```
./StSm_plot.sh [additional_arguments]
```

### 3. `StSm_PrepColl.sh`

**Purpose:**  
This script sets up a new collision directory with a unique name and prepares it for running a StarSmasher simulation. It generates a random directory name prefixed with "Collision_", copies necessary files from the tools directory into the new directory, and moves the 3D MESA SPH initial profiles if they exist. If the initial profiles are not found, the script exits with an error message.

**Usage:**  
Run the script to create a new collision directory and set up the environment:
```
./StSm_PrepColl.sh
```

### 4. `StSm_rm.sh`

**Purpose:**  
This script is designed to selectively remove snapshots from a StarSmasher simulation run. You can specify a range of snapshot numbers to remove, which helps manage disk space by deleting unnecessary files.

**Usage:**  
Specify the range of snapshots to remove by providing the start and end numbers as arguments:
```
./StSm_rm.sh [start] [end]
```
For example, to remove snapshots between 250 and 3022:
```
./StSm_rm.sh 250 3022
```

### 5. `StSm_run.sh`

**Purpose:**  
This script runs the StarSmasher simulation using MPI (Message Passing Interface). It logs the start and end times of the simulation, captures output and error logs with unique filenames based on the current timestamp, and checks for successful completion of the run.

**Usage:**  
Run the script to start the StarSmasher simulation:
```
./StSm_run.sh
```

The script will automatically log the details of the run and handle the output and error files.
