Welcome to the installation guide for StarSmasher!  Please feel free to contact me (jamie.lombardi at allegheny.edu) for help or to provide feedback.

# [0] What you'll need to install StarSmasher
You'll need a few things to prepare for a clean installation of StarSmasher.

## [0.1] A Linux computer, ideally with an NVIDIA graphics card

First, you should have a computer running Linux.  An NVIDIA graphics card is worth having, since the GPU takes over the gravity calculation, but it is not essential: ``make cpu`` builds a version that runs anywhere, and sections [0.2] and [0.3] below can then be skipped.  The instructions in this section [0] below assume that you'll be setting up your computer *without* using modules.  If you are working on a research cluster of a university, for example, then you may be able to set up the necessary software with ``module load`` commands.  For example ``module load cuda/cuda-10.1.2`` would load version 10.1.2 of cuda, while ``module load mpi/openmpi-1.10.5-gcc-6.4.0`` would load version 1.10.5 of openmpi (compiled with version 6.4.0 of gcc), if available.  To see what modules are available on your system, type ``module avail``.  **If you are using modules, then load up cuda and openmpi, and move on to section [1] of this installation guide!**

## [0.2] nvcc (the NVIDIA cuda compiler)

The gravity library in StarSmasher uses the NVIDIA graphics card and is compiled using the NVIDIA cuda compiler nvcc.  To check if nvcc is installed, type
```
nvcc --version
```
**If the version information is returned, then nvcc is properly installed and you can skip to section [0.3]!**  If instead you know that your system will need cuda installed, then follow NVIDIA's installation guide to do so:

https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

--------------------
**Troubleshooting cuda compiler usage**

If ``nvcc --version`` returns a command not found error message, it's possible that ``cuda`` is installed but your $PATH environment variable is improperly set.  If you think that may be the case, try to find where nvcc is located, for example by using
```
locate nvcc
```
To check what directories are in $PATH, use
```
echo $PATH
```
Let's say for the sake of argument that you find nvcc exists in ``/usr/local/cuda-11.4/bin`` but that this directory is not in $PATH.  You should then update your $PATH and corresponding $LD_LIBRARY_PATH variables.  In bash, this can be done with
```
export PATH=/usr/local/cuda-11.4/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.4/lib64:$LD_LIBRARY_PATH
```
Such commands could be placed in your ~/.bashrc file or the equivalent so that they wouldn't need to be executed by hand in every session.

## [0.3] NVIDIA driver communication

The cuda installation should enable communications with the graphics card via the NVIDIA driver.  To test this, type
```
nvidia-smi
```
**If information about your graphics card(s) is returned, then the driver is properly functioning and you can skip to subsection [0.4]!**  

-------------------------------------------------
**Troubleshooting NVIDIA driver communications**

If cuda has been installed but ``nvidia-smi`` doesn't work, then double check that you have completed the post installation steps from the NVIDIA installation guide.  Unfortunately, a reboot after a fresh cuda installation may be the easiest way to cure driver communication issues.  If the machine stalls on reboot, make sure that in the BIOS the EFI Secure Boot is disabled.

## [0.4] openmpi and mpif90

StarSmasher is parallelized using openmpi.  To check if the openmpi fortran compiler is properly installed, type
```
mpif90 --version
```
**If version information is returned, then you can proceed to section [1]!**  If instead you know that your system will need openmpi installed, then execute one of the following two commands, depending on which type of Linux you're using.  For CentOS/RHEL/Fedora and other similar systems, use
```
sudo dnf install openmpi
```
In Ubuntu or other Debian Linux operating systems, use
```
sudo apt install openmpi-bin
```

--------------------------------
**Troubleshooting mpif90 usage**

If ``mpif90 --version`` returns a command not found error message, it's possible that ``openmpi`` is installed but your $PATH environment variable is improperly set.  If you think that may be the case, try to find where mpif90 is located, for example by using
```
locate mpif90
```
To check what directories are in $PATH, use
```
echo $PATH
```
Let's say for the sake of argument that you find mpif90 exists in ``/usr/lib64/openmpi/bin`` but that this directory is not in $PATH.  You should then update your $PATH environment variables.  In bash, this can be done with
```
export PATH=/usr/lib64/openmpi/bin:$PATH
```
This command could be placed in your ~/.bashrc file or the equivalent so that it wouldn't need to be executed by hand in every session.

# [1] Obtaining StarSmasher

## [1.1] With GitHub
Just go to the StarSmasher’s GitHub page and download the .zip file. Put it in the desired folder and extract it. You will get all the files.

## [1.2] With git

On your workstation, clone the repo.
```
git clone https://github.com/jalombar/starsmasher.git
```

------------------------------------------------------------------------
Now that you have the code, let's take a quick look in the StarSmasher folder.

```
cd starsmasher
ls
```

This main directory contains seven subdirectories:
* Blackollider: a modified StarSmasher's code optimized to make collisions Between Black Holes / Neutron Stars and Stars (Both Main Sequence, pre Main Sequence, Giants...)
* Data_visualization: a tutorial on reading and displaying the output snapshots from StarSmasher
* documentation: installation and other instructions
* example_input: contains example input files for various cases
* misc: miscellaneous files such as assorted makefiles, the StarCrash manual, and a tree gravity GPU library
* parallel_bleeding_edge: latest version of code (this is probably what you want)
* splot_directories: post-simulation code used for data extraction and analysis

# [2] Compiling StarSmasher

For the code itself, the directory you care about is parallel_bleeding_edge, and that is what section [2] below uses throughout.  There is a second version, Blackollider, which differs mainly in how smoothing lengths are dynamically evolved and is the better choice if you plan to run simulations involving compact objects treated as point masses.  Build parallel_bleeding_edge first; section [2.3] explains what changes for Blackollider.

## [2.1] Compile the gravity library

To compile the gravity library, cuda will need to be installed.  If you have no NVIDIA card, skip this section entirely and go to [2.2], where ``make cpu`` builds a version that does not use the GPU.  You can find the gravity library in both parallel_bleeding_edge/src/SPHgrav_lib2 and Blackollider/src/SPHgrav_lib2/.  The SPHgrav_lib2 subdirectory contains code written by Evghenii Gaburov (and somewhat modified by Jamie Lombardi and Sam Knarr) for calculating softened gravitational forces and potentials on NVIDIA GPUs. From within the main starsmasher folder,
```
cd parallel_bleeding_edge/src/SPHgrav_lib2
```
will navigate you to this directory, where the library will be built.  In the SPHgrav_lib2 subdirectory, there are few files including a Makefile.  This Makefile works out which GPU you have, so normally there is nothing in it to edit.

Your card's compute capability decides which machine code nvcc generates, and the Makefile asks the driver for it directly:
```
COMPUTE_CAPABILITY := $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -1 | tr -d '.')
```
You can see what it will find by running the same query yourself:
```
nvidia-smi --query-gpu=name,compute_cap --format=csv,noheader
NVIDIA GeForce RTX 5080, 12.0
```
A compute capability of 12.0 means the library is built with ``-arch=sm_120``.

If ``nvidia-smi`` is not available, on a login node without a GPU for instance, the Makefile falls back to ``-arch=sm_89``. Rather than editing the file, pass the two digits on the command line:
```
make COMPUTE_CAPABILITY=75
```
Do not override ``NVCCFLAGS`` directly for this. A value given on the command line replaces the whole variable, so ``make NVCCFLAGS="-arch=sm_75"`` quietly drops the optimisation level and the include paths along with it.

You can test the compilation of the gravity library by typing
```
make
```
Any warnings about "invalid narrowing conversion" can be ignored.  Successful compilation results in the creation of the library file libSPHgrav.a.

----------------------------------------------
**Troubleshooting compilation of the gravity library**

If compilation of the gravity library fails, the most likely explanation is that the Makefile is not identifying the correct location of the nvcc executable and/or cuda libraries.  Within the Makefile, look for the line
```
CUDAPATH       := $(shell dirname $(shell dirname $(shell which nvcc)))
```
and change it so that it identifies the main cuda directory.  For example, let's say that your system contains the files /usr/local/cuda-11.2/bin/nvcc and /usr/local/cuda-11.2/lib64/libcudart.so.  Then your main cuda directory is /usr/local/cuda-11.2, and you can change the above line in the makefile to
```
CUDAPATH       := /usr/local/cuda-11.2
```

### [2.2] Finish compiling StarSmasher

Now that your SPHgrav_lib2 is ready, go back a level to the src folder:
```
cd ..
```
If you skipped [2.1] because you have no NVIDIA card, go to that same src folder directly, with ``cd parallel_bleeding_edge/src`` from the main starsmasher folder.

The StarSmasher code lives here, along with several example makefiles. To find out if the default Makefile will work for you, type

```
make
```
for the GPU version, or
```
make cpu
```
for the version that does not need one.
If StarSmasher compiles properly, in your terminal will appear this phrase:
```
***MADE VERSION THAT USES GPUS***
```
When you complete the compilation, it will be create an executable file ending with "\_gpu\_sph." This executable is automatically moved to the directory *above* src. To run StarSmasher, you will have to move to that folder. Then just open your terminal where the .exe file is present or, after the compilation, type
```
cd ..
```
StarSmasher is ready to use! 
**You can now proceed to the optional section [3], or just start to use StarSmasher.**  To run your first simulations, follow the tutorial “How to create a MESA star” situated in the “example_input” folder!

----------------------------------------------
**Troubleshooting compilation of StarSmasher**

If the compilation is unsuccessful, then your setup is either missing or unable to find at least one library or executable.  A few other makefiles exist in the src directory.  They're unlikely to work, but you can try them with ``make -f filename``, and they may provide some inspiration as for what to change.

Start by asking the Makefile what it found:
```
make config
```
which on the machine these instructions were written on reports
```
MPI
  mpif90          /usr/lib64/openmpi/bin/mpif90
  wrapping        GNU
  fixed-form flag -ffixed-line-length-132 -fallow-argument-mismatch
CUDA
  nvcc            /usr/local/cuda/bin/nvcc
  CUDAPATH        /usr/local/cuda
  runtime libdir  /usr/local/cuda/lib64
Compilation
  FFLAGS          -O4 -mcmodel=medium
  code model      -mcmodel=medium
Output
  GPU executable  ../parallel_bleeding_edge_gpu_sph
  CPU executable  ../parallel_bleeding_edge_cpu_sph
```
Compare that against your own installation. Whichever line is wrong tells you which variable to set, and you set it on the command line rather than by editing the Makefile:

| If this is wrong | Set |
|---|---|
| the mpif90 path, or the compiler it wraps | ``make MPIF90=/full/path/to/mpif90`` |
| CUDAPATH, or the runtime libdir | ``make CUDAPATH=/opt/cuda`` |
| the optimisation flags | ``make OLEVEL=-O2`` |

If ``nvcc`` is reported as not found, the GPU build will stop and say so before it reaches the link. That is the expected result on a machine with no CUDA, and ``make cpu`` will build a working, slower executable instead.

Another possible compilation issue can occur if you've ported the source code from one system to another without removing or recreating the .o object files created on the previous system.  In this scenario, start the compilation process fresh by first doing
```
make clean
```
to remove any .o and other similar object files.

If you find a way to make the main makefile more robust, please consider sharing it with me (jamie.lombardi at allegheny.edu).

### [2.3] Building Blackollider instead

Blackollider is built exactly as above, with ``Blackollider`` in place of ``parallel_bleeding_edge`` in every path.  Its gravity library lives in ``Blackollider/src/SPHgrav_lib2`` and its source in ``Blackollider/src``, and the same ``make``, ``make cpu`` and ``make config`` commands apply.

Use it if you are running collisions involving compact objects treated as point masses.  Otherwise parallel_bleeding_edge is the one to build.

# [3] Installing SPLASH for visualization (optional but encouraged)

For visualization purposes, install Daniel Price's SPLASH.  The details are given in the installation guide

https://splash-viz.readthedocs.io/en/latest/getting-started.html

and here we summarize the main steps.  If you don't have root privileges, then install SPLASH in your home directory and set the appropriate environment variables as described in the installation guide.  If you do have root privileges and want to install SPLASH system-wide in, say, /usr/local/src/ using a local giza, then execute the following. 
```
cd /usr/local/src  
sudo git clone https://github.com/danieljprice/splash.git
cd splash
sudo git clone https://github.com/danieljprice/giza.git
sudo make SYSTEM=gfortran withgiza
sudo make install
```
If the compilation of SPLASH complains about not finding cairo.h, then try “sudo dnf install cairo-devel” on a CentOS/RHEL/Fedora system.  

Alternatively, if you want to first install cairo and giza system-wide and then install Splash after that, then execute the following as root.
```
cd /usr/local/src
svn co http://svn.code.sf.net/p/giza/code/ giza
cd giza
sudo yum install cairo-devel
./configure
make
make install
cd /usr/local/src  
sudo git clone https://github.com/danieljprice/splash.git
cd splash
sudo make SYSTEM=gfortran
sudo make install
```

To use SPLASH to look at out*.sph snapshots created by StarSmasher, use
```
splash -f starsmasher out*.sph
```

-----------------------------------------
**Troubleshooting the running of SPLASH**

If at run time SPLASH cannot find the giza library, then it will fail with an error message such as
```
splash: error while loading shared libraries: libgiza.so.0: cannot open shared object file: No such file or directory
```
(It's possible that the name of the shared library could be slightly different, but in what follows we assume the name libgiza.so.0.)
The fix is to find where libgiza.so.0 is located (perhaps with ``locate libgiza.so.0``) and update your LD_LIBRARY_PATH accordingly.  For example, if libgiza.so.0 is located in, say, /usr/local/lib then you would execute
```
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```
Placing this command in ~/.bashrc or the equivalent would help so that it wouldn't have to be executed by hand in every terminal session that uses SPLASH.

We wish you happy times smashing stars!
