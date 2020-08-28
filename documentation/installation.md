# Before installing StarSmasher
Before installing StarSmasher, you need few things to make a clean installation and run the program.  
The 1st thing is that you must have one or more NVIDIA Graphic Cards. Why? 
There are 2 version of the code. The GPU version and the CPU version. The GPU version of the code uses Nvidia graphic card to calculate gravity and let to the CPU the rest (and less) of the work. The CPU version of the code make the CPU handle all the calculations.	
Depending on your machine, in most cases (unless you have a very high number of CPU per node) the GPU version of the code is much faster than the CPU version of the code. To have a comparison, my laptop has an NVIDIA GeForce GTX 950M, and I can handle 150000 Particles with the GPU version of the code. Instead, with the CPU version (and 2 quad-core intel i7), I can only get 15000 particles at the same amount of time respect to the GPU version.	
Then you must need an NVIDIA Graphic Card to run StarSmaser, not others. If you don't have get a computer or a machine with NVIDIA Graphic card.  
You must also be sure that if you have more than one NVIDIA Graphic Card they MUST BE of the same version. What i mean?
In this wikipedia page
https://en.wikipedia.org/wiki/CUDA
there is written the “compute capability version” of each NVIDIA’s Graphic Card. You must be sure that if you have differents NVIDIA Graphic cards in your machine, they must be of the same compute capability version. For example, in your machine you can have a Nvidia TITAN Xp, a Titan X and a GeForce GTX 1080 Ti, because they are all 6.1 compute capability version and then they will calculate gravity correctly. You can’t have, however, a Titan X and a NVIDIA TITAN V, because one is a 6.1, the other is a 7.0. This means that they will calculate gravity both in a different way and StarSmasher will give you a very bad simulation where eventually the stars will explode. Then, please, be sure that all of your graphic cards are of the same compute capability version if not the same.

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

### Prepare your machine to compile the StarSmasher's GPU version
Before installing the GPU version of StarSmasher, you must prepare your machine/PC. If you are running on Linux (Ubuntu especially), the first thing that you have to do is to install the correct driver for your Graphic card (even more than one if you have different graphic cards). To do that, I recommend you to follow this guide (Automatic Install using standard Ubuntu Repository, if you have Ubuntu):

https://linuxconfig.org/how-to-install-the-nvidia-drivers-on-ubuntu-18-04-bionic-beaver-linux

There is an error, however. The guide says that the system detected a NVIDIA GeForce GTX 680 and that the recommended driver to install is the nvidia-384. This is not correct. The system detected a GeForce GTX 1060 6GB and then the correct driver to install was the nvidia-390. Then, please, keep attention!
Now you have to install the NVIDIA nvcc toolkit. There are numerous ways to install it, you can check the NVIDIA official guide:

https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

Or another NVIDIA official’s guide if you are not installing StarSmasher in Ubuntu. However, the easiest way to install nvcc in Ubuntu is just this. Open your terminal and run this command: 
```
Sudo apt get install cuda-toolkit
```
But keep attention. This command is going to install you the latest nvcc version that NVIDIA developers installed in the Ubuntu's libraries. If you have, for example, the latest NVIDIA graphic card, there is the possibility that NVIDIA developers didn’t updated nvcc to the latest version for your graphic card. This means that you have to install nvcc with a specified version, as there is written in the NVIDIA official guide. There is the possibility, then, that you could get in trouble during installation. My personal suggestion is to get an NVIDIA graphic card that has a high gravitation’s compute performance. A very new graphic card doesn’t mean that is less powerful to calculate gravity respect to the old one. You can check it there:

https://gpu.userbenchmark.com/Compare/Nvidia-Titan-RTX-vs-Nvidia-RTX-2080-Ti/m664199vs4027

As you can see the NVIDIA Titan RTX is newer than the NVIDIA RTX 2080 Ti, However, about gravity calculation, is just a 10% more performant, there is like no difference then. In this case the NVIDIA RTX 2080 Ti is a better choice because it’s older and there surely is the nvcc version ready for that graphic card.
Once your installation is complete, you now have a nvcc file in your Ubuntu OS. Search it and keep in mind its location.


### Compile the gravity library (CUDA) and StarSmasher

If you want to use the CPU version of StarSmasher, avoid this part of the tutorial and just go to the next paragraph.

The SPHgrav_lib subdirectory contains code written by Evghenii Gaburov (and somewhat modified by Jamie Lombardi and Sam Knarr) for calculating softened gravitational forces and potentials on NVIDIA GPUs.
Now you must do some important changes in this. All depends on your OS. In the SPHgrav_lib there are few files including some makefiles. Every makefile is written of an OS. The Standard “makefile” is written for Keeneland. If you are using Ubuntu, you will use the makefile.ubuntu makefile, if you are using Quest, you will the makefile.quest and so on. You will find every makefile in the “misc” folder if there isn’t in the SPHgrav_lib. If there isn’t even there, you must substitute in the standard makefile the path to nvcc. Once you have chosen your makefile (i did a makefile called "MakefileUbuntuCUDAToolkit" if you are following this guide spet by step), just change the name in “makefile” and delete the others, you don’t need them.	
In the makefile (Ubuntu case) there is written this string:

```
CUDAPATH       := /usr/lib/nvidia-cuda-toolkit/
```

In Ubuntu this is the path where there is the nvcc file if you installed it withe the sudo apt-get command. Be sure that there is! If not because you isntalled the CUDA toolkit in anoher way, edit this path!
In the same file there is also this string to edit:

```
NVCCFLAGS := -arch=sm_61
```

As there is written, this string is for the NVIDIA GTX 1070. The number 61 indicate the computability version of the NVIDIA Graphic card. If you check there:

https://en.wikipedia.org/wiki/CUDA

the NVIDIA GTX 1070 has a computability version of 6.1. then in the string there is written “NVCCFLAGS := -arch=sm_61”. If you have a different graphic card, this number must be change according to the computability version as the Wikipedia page says, or the gravity won’t be calculated and your stars will explode. Then, if you have a GeForce GTX 950m, like me, write:

```
NVCCFLAGS := -arch=sm_50
```

Because my graphic card is a 5.0, if you have a NVIDIA TITAN RTX, write:

```
NVCCFLAGS := -arch=sm_75
```

because it has a computability version of 7.5, and so on. This is why all of your graphic cards must be of the same computability version, or Gravity won’t be calculated well and your stars will eventually explode.


Now that the Makefile of the SPHgrav_lib folder is ready, we can return to the "src" folder.

### Go there if you want to install the CPU version of the code

Now that you SPHgrav_lib is complete (only if you are installing the GPU version of StarSmasher), go to the src folder . In this folder there are a ton of makefiles. As before, choose the correct one for you. For Ubuntu users that want to install StarSmasher with OpenMPI  (this is the easiest way, you can also install it with IFort) there is an already written file called "MakefileGPUubuntu", if you want to use the GPU version (recommended), and "MakefileCPUubuntu" if you want to use the CPU version (unrecommended as explained earlier).  
Then, if you have Ubuntu, just choose one of these two files and rename it as “Makefile”. Be also sure that there are not other files called only “Makefile”. Delete the others makefiles that you don’t need if you want. If you are not an Ubuntu user, be sure to edit your makefile in way that it contains the installation path only for CPU/GPU and OpenMPI/IFort. I suggest you to follow the structure of "MakefileGPUubuntu" or "MakefileCPUubuntu" to do a clean EDIT.	
Once your Makefile is done, everything is ready. Open your terminal in your src folder and type:

```
Make
```

The make command is going to follow all the instructions of the makefile. If make goes in error that's because you are missing one or more libraries. Then, as explained before, you need to install them and Ubuntu is the easiest way to do that.

If the GPU version of StarSmasher is been installed correctly, in your terminal will appear this phrase:

```
***MADE VERSION THAT USES GPUS***
```

If the CPU version of StarSmasher is been installed correctly, in your terminal will appear this phrase:

```
***MADE VERSION THAT DOES NOT NEED GPUS***
```

Now StarSmasher is ready to use! When you complete the installation, it will be created a .exe file called “test_gpu_sph if you installed the GPU version of the program, or “test_cpu_sph” if you installed the CPU version of the code. This .exe is situated in the upper folder of src. To run StarSmasher, you will have to move on that folder. Then just open your terminal where the .exe file is present or, after the installation, just type
"cd .."

Now, to run your first simulations, follow the tutorial “How to create a MESA star” situated in the “example_imput” folder!


# Importants suggestions if you are running StarSmasher on your Hard Disk
A last important advice. When your simulation is finished and you want to transfer your files in another Hard Disk or similar, never cut and paste your files. StarSmasher’s files and snapshots are very sensible, and there is the risk that during the cut and paste process you may lose or damage your files. To avoid this, just copy and paste them, be sure that every files is there and that is not been damaged (like visioning the texts or try to visualize them with SPLASH, as you will see next) and then, eventually, delete the original files once the transfer is been completed. 	
If you want instead to continue your simulation in another Hard Disk, be sure to recompile the code.

```
Make clean
```
Then

```
Make
```
Now the code is been recompiled and the mpirun command should run perfectly. 

We suggest to the user to run StarSmasher into an Hard Disk. If you are going to use a lot of particles (100000 or more) each snapshot of your simulation is going to be around 40 mb. Then, if you have thousands of them, there is the possibility that all your space is going to be occupied. To avoid this, we suggest a NEW Hard Disk of 1 TeraByte or more. Why a New? There is not any known risk (actually) that StarSmasher is going to delete your personal files on your Hard Disk, but there is a problem with Linux and the file systems of the Hard Disk itself. If you are going to run StarSmasher in Ubuntu and try to install StarSmasher in an Hard Disk that has FAT32 file system, when you'll try "mpirun -np N test_gpu_sph", OpenMPI is going to say you that it is not possible to find the executable test_gpu_sph. This is because Linux won't give the permission to run the executable in a space that has not the same file system of linux. However, everything will run perfectly if your Hard Disk is going to be have NTFS's file system. Then, why a new Hard Disk? Because so you can format it with the correct file system of your linux OS and you won't lose any of your personal data into the Hard Disk if you can't do a backuo or for some other reason.

The solution of having a NTFS's Hard Disk is not universal and it may change in base of your OS or other factors.


We wish to the user a good use of StarSmasher, actually best star collider in the planet!

