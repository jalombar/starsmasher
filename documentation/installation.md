# Before installing StarSmasher
Before installing StarSmasher, you need few things to make a clean installation and run the program.  
The 1st thing is that you must have one or more NVIDIA Graphic Cards. Why? 
There are 2 version of the code. The GPU version and the CPU version. The GPU version of the code uses Nvidia graphic card to calculate gravity and let to the CPU the rest (and less) of the work. The CPU version of the code make the CPU handle all the calculations.	
Depending on your machine, in most cases (unless you have a very high number of CPU per node) the GPU version of the code is much faster than the CPU version of the code. To have a comparison, my laptop has an NVIDIA GeForce GTX 950M, and I can handle 150000 Particles with the GPU version of the code. Instead, with the CPU version (and 2 quad-core intel i7), I can only get 15000 particles at the same amount of time respect to the GPU version.	
Then you must need an NVIDIA Graphic Card to run StarSmaser, not others. If you don't have get a computer or a machine with NVIDIA Graphic card.

The 2nd thing is that the CPU version is actually bugged and you won’t be able to create a star. During the relaxation of it from a MESA profile (as you will see in the tutorial “How to create a Star using MESA”), the Star will just explode if you use more than 6000 particles. Unless the bug will be solved, we suggest you to use the GPU code of the program, because it’s much faster and safe.
The 3rd thing that we suggest is to use Linux as an OS to run StarSmasher, especially Ubuntu due to the ease of installation of any libraries. This doesn’t mean that you won’t be able to install StarSmasher in others Linux OS, MacOS or Windows, just know that you can met difficulties and waste a lot of time if you are not skilled and experienced on installing libraries.	
For example, I tried (for 2 months) to install StarSmasher in Windows via Cygwin with a very bad result. It’s been impossible to install the GPU version of the code (only the CPU after lot of difficulties via the incompatibilities of the Cygwin’s libraries) and it’s not been possible to run StarSmasher cause of some weirds errors of OpenMPI library that’s been impossible to solve. Instead with Ubuntu the installation went perfectly fine in minutes. The advantage to use Ubuntu to run StarSmasher is that, if you miss some library during installation and get error, the terminal will advice you that the “XY (generic command) command is not been found”. Then if you are going to type "XY" in the terminal you are going to know that the specific XY command is present in the Z package (a generic package) and you will be able to install it with the command 

```
Sudo apt-get install Z
```
Then after the installation of the Z package you will be able to continue the StarSmasher installation in a very clean way and without any problems. Instead, if you have windows or other OS, you can find problems, then you have to understand where is it, search where to get the specific command and other troubles that will only make you waste time unless you are not skilled to do that kind of things. This is why we suggest you to use Linux Ubuntu. If you don’t have it, just do a bipartition to your PC or Machine. You just need 30-40 GBs free and then we suggest you to use an external Hard Disk with few Tera-Bytes of space.	
If you are going to use a lot of particles (100k or more) each snapshot of your simulation is going to weight more than 30 mb. Then if you will make it run for thousands and thousands of Snapshots, well, you will need space. If you miss it StarSmasher is going to go in a weird error during the run. This is why we suggest you to make some runs on your Ubuntu Space (that part of your computer HDD where Ubuntu is installed) and then install it directly on the Hard Disk.	
A last advice. It’s impossible to run StarSmasher on VirtualBox. This is because StarSmasher need real CPU and GPU to run, not a virtual one. You will be able to install it on VB (only the CPU version) but if you’ll try to run it, mpirun will go in error. Then, please, do a bipartition!


# OBTAINING AND COMPILING

## With GitHub
Just go to the StarSmasher’s GitHub page and download the .zip file. Put it in the desired folder and extract it. You will get all the files. To travel between folders, you can simply open your terminal (Linux) in the specified folder. Once you have your files, you are ready for the installation.


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

# Importants suggestions if you are running StarSmasher on your Hard Disk
The last important advice. When your simulation is finished and you want to transfer your files in another Hard Disk or similar, never cut and paste your files. StarSmasher’s files and snapshots are very sensible, and there is the risk that during the cut and paste process you may lose or damage your files. To avoid this, just copy and paste them, be sure that every files is there and that is not been damaged (like visioning the texts or try to visualize them with SPLASH, as you will see next) and then, eventually, delete the original files once the transfer is been completed. 	
If you want instead to continue your simulation in another Hard Disk, be sure to recompile the code, or there is the risk that, after entering the command to start the simulation, mpirun will say you that’s not been possible to find the executable “test_gpu_sph” (or “test_cpu_sph”, depending on what version of the code are you using). To solve this, just open the terminal on the src folder and type

```
Make clean
```
Then

```
Make
```
Now the code is been recompiled and the mpirun command should run perfectly. If this doesn’t happen and you still get the error, you ran into a bug. To solve it just transfer back your files to the original Hard Disk. Just know that this last bug doesn’t make any sense, but just happened to me.

