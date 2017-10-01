# OBTAINING AND COMPILING

## With git

On your workstation, clone the repo.
```
git clone https://github.com/jalombar/starsmasher.git
```

## Installing on Allegheny
To get the code, first get an account on carrv109.allegheny.edu and ask Jamie to open up connections to carrv109 from your preferred IP address.
Then install with git and go to the install folder (with a default name of 'starsmasher').
```
cd starsmasher
```

This main directory contains at least four subdirectories:

* example_input: contains example input files for various cases
* misc: miscellaneous files such as assorted makefiles and the StarCrash manual
* parallel_bleeding_edge: latest version of code (this is probably what you want)
* splot_directories: post-simulation code used for data extraction and analysis

For runs on a gpu-enabled cluster, the directory you care about is parallel_bleeding_edge.
Let's try to compile the code in there.

```
cd parallel_bleeding_edge/src/SPHgrav_lib
```

## Compiling the GPU-enabled code

### Compile the gravity library (CUDA)
The SPHgrav_lib subdirectory contains code written by Evghenii Gaburov (and somewhat modified by Jamie Lombardi and Sam Knarr) for calculating softened gravitational forces and potentials on NVIDIA GPUs.
If you are running on something other than keeneland, you will need to customize the Makefile for your system: any path that contains "cuda" or "CUDA" probably needs to be changed.
To make the library, simply issue the command

```
make
```

(The final thing that make does is "./makelib.sh" which is what will make the library.)
If it doesn't work, and you want to have a fresh start, first use "make clean".
If all goes well the library libGPUsph_gatherscatter.a will be generated.

### Compile starsmasher
The SPH source code is in the .. directory.
This is a parallel code.
Once the GPU library has been compiled, we will move to the source code directory, copy the appropriate makefile (assuming keeneland makefile for this) to 'makefile', and compile with make.

```
cd ..
cp makefile.keeneland makefile
make
```

You will likely need to customize the makefile for your computer.  
For example, you may need to change the compiler or update paths to point to the appropriate locations.
The executable name will end with _sph and will depend on the name of the directory that src appears within.
The executable will automatically be moved up a level in the directory structure (that is, to the ".."  directory).  
You probably want to remove the executable for now so that you don't later copy it to other places where you don't want it:

```
cd ..
rm *_sph
```

