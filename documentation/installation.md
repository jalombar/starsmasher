# [0] What you'll need to install StarSmasher
Before installing StarSmasher, you need a few things to make a clean installation and run the program.

## [0.1] A Linux computer with an NVIDIA graphics card

First, you should have a computer running Linux with at least one NVIDIA graphics card. While it may be possible to install StarSmasher in an operating system other than Linux, that's not yet been tested and could be challenging.

## [0.2] cuda

The gravity library in StarSmasher uses the graphics card and is compiled in cuda using the NVIDIA compiler nvcc.  To check if nvcc is installed, type
```
nvcc --version
```
**If the version information is returned, then nvcc is properly installed and you can skip to section [0.3]!**  If instead you know that your system will need cuda installed, then follow NVIDIA's installation guide to do so:

https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

### [0.2.1] Troubleshooting cuda

If ``nvcc --version`` returns a command not found error message, it's possible that ``cuda`` is installed but your $PATH environment variable is improperly set.  If you think that may be the case, try to find where nvcc is located, for example by using
```
locate nvcc
```
To check what directories are in $PATH, use
```
echo $PATH
```
Let's say for the sake of argument that you find nvcc exists in ``/usr/local/cuda-11.4/bin`` but that this directory is not in $PATH.  You should then update your $PATH and correspondingly your $LD_LIBRARY_PATH variables.  In bash, this can be done with
```
export PATH=/usr/local/cuda-11.4/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-11.4/lib64:$LD_LIBRARY_PATH
```
Such commands could be placed in your ~/.bashrc file or equivalent so that they wouldn't need to be executed by hand in every session.


Note that there are post installation steps and a reboot after installation can help.  If the machine stalls on reboot, check in the BIOS that the EFI Secure Boot is disabled.


# [1] OBTAINING

## [1.1] With GitHub
Just go to the StarSmasher’s GitHub page and download the .zip file. Put it in the desired folder and extract it. You will get all the files.

## [1.2] With git

On your workstation, clone the repo.
```
git clone https://github.com/jalombar/starsmasher.git
```

# [2] COMPILING
Now that you have the code, let's go in the StarSmasher folder!

```
cd starsmasher
```

This main directory contains seven subdirectories:
* Blackollider: a modified StarSmasher's code optimized to make collisions Between Black Holes / Neutron Stars and Stars (Both Main Sequence, pre Main Sequence, Giants...)
* Data_visualization: a tutorial on reading and displaying the output snapshots from StarSmasher
* documentation: installation and other instructions
* example_input: contains example input files for various cases
* misc: miscellaneous files such as assorted makefiles, the StarCrash manual, and a tree gravity GPU library
* parallel_bleeding_edge: latest version of code (this is probably what you want)
* splot_directories: post-simulation code used for data extraction and analysis

For the code itself, the directory you care about is either parallel_bleeding_edge or Blackollider.  These two version of the code differ primarly in how smoothing lengths are dynamically evolved.  If you plan to run simulations involving compact objects treated as point masses, then Blackollider is the better choice.

## [2.1] Prepare your machine to install the GPU's version of StarSmasher

Before installing StarSmasher, you must prepare your machine.
To compile the gravity library, cuda will need to be installed.  You can find the gravity library in either parallel_bleeding_edge/src/SPHgrav_lib2 or Blackollider/src/SPHgrav_lib2/.  For example,
```
cd Blackollider/src/SPHgrav_lib
```
will navigate you to the directory where the library can be built.



### [2.2] Compile the gravity library (CUDA) and StarSmasher

The SPHgrav_lib2 subdirectory contains code written by Evghenii Gaburov (and somewhat modified by Jamie Lombardi and Sam Knarr) for calculating softened gravitational forces and potentials on NVIDIA GPUs.
In the SPHgrav_lib2 there are few files including some makefiles.  The standard “makefile” should hopefully be sufficient after making one change, as described below. 
In particular, look for this string to edit:

```
NVCCFLAGS := -arch=sm_61
```

As written, this string is for an NVIDIA card with compute capability (version) 6.1, such as an NVIDIA GTX 1070. You can look up the compute capabilities of NVIDIA cards here:

https://en.wikipedia.org/wiki/CUDA

If your graphics card has a different compute capability, then you would want to change the "61" portion of this line accordingly.  For example, if you have a GeForce GTX 950M, which has a compute capability of 5.0, then write:

```
NVCCFLAGS := -arch=sm_50
```

Or if you have, for example, a NVIDIA TITAN RTX, then write

```
NVCCFLAGS := -arch=sm_75
```

because it has a computability version of 7.5, and so on.

Now that the Makefile of the SPHgrav_lib folder is ready, we can return to the "src" folder.

### [2.3] Setup your src's Makefile.

Now that you SPHgrav_lib2 is ready, go back a level to the src folder. In this directory there are several example makefiles. As before, the default Makefile will hopefully work for you. 

To find out, open your terminal in your src folder and type

```
make
```

The make command is going to follow all the instructions of the makefile. If make goes in error that's because you are missing one or more libraries. Then, as explained before, you need to install them and Ubuntu is the easiest way to do that.

For example, a problem that you can encounter is that your machine could not find the command "mpif90" during the compilation. To solve this problem, just type in your Ubuntu's terminal:


```
mpif90
```
It will say you that:

```
Command «mpif90» not found, but it can be installed with:

sudo apt install mpich        # version 3.3.2-2build1, or
sudo apt install openmpi-bin  # version 4.0.3-0ubuntu1
```
Then, to solve this problem, just type:

```
sudo apt install openmpi-bin
```
And Ubuntu will install the libraries you need in minutes. Then, if you had this little trouble, to install StarSmasher correctly type again:

```
make
```
Now the terminall will compile StarSmasher in tens of seconds.


If the GPU version of StarSmasher is been installed correctly, in your terminal will appear this phrase:

```
***MADE VERSION THAT USES GPUS***
```

Now StarSmasher is ready to use! When you complete the installation, it will be created an executable file called “test_gpu_sph. This executable is situated in the upper folder of src. To run StarSmasher, you will have to move to that folder. Then just open your terminal where the .exe file is present or, after the installation, just type:

```
cd ..
```

Now, to run your first simulations, follow the tutorial “How to create a MESA star” situated in the “example_imput” folder!


# [3] Importants suggestions if changing computer systems

If you want to do a new install, start the compilation process fresh by first doing

```
make clean
```
to remove any .o and other similar object files.

We wish to the user a good use of StarSmasher!
