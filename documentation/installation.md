# [0] Before installing StarSmasher
Before installing StarSmasher, you need a few things to make a clean installation and run the program.

First, you must have at least one NVIDIA graphics card. If you have more than one card, they should be of the same compute capability (version): this wikipedia [page](https://en.wikipedia.org/wiki/CUDA) gives the compute capability of each NVIDIA graphic card.

[//]: # "Why? There are 2 version of the code. The GPU version and the CPU version. The GPU version of the code uses Nvidia graphic cards to calculate gravity and let the CPU do the rest (and less) of the work. The CPU version of the code has the CPU handle all the calculations."

[//]: # "Depending on your machine, in most cases (unless you have a very high number of CPU per node) the GPU version of the code is much faster than the CPU version of the code. For example, an NVIDIA GeForce GTX 950M can handle 150,000 Particles with the GPU version of the code in the same amount of time that the CPU version handles 15,000 particles on a dual quad-core intel i7).  Then you must need an NVIDIA graphics card to run StarSmasher, not others."

[//]: # "If, for example, your machine has an Nvidia TITAN Xp, a Titan X and a GeForce GTX 1080 Ti, then they will calculate gravity correctly because they are all a 6.1 compute capability version. You can’t have, however, a Titan X and a NVIDIA TITAN V, because one is a 6.1, the other is a 7.0. This means that they will calculate gravity both in a different way and StarSmasher will fail to run properly. So please be sure that all of your graphic cards are of the same compute capability version if not the same."

[//]: # "The 2nd thing is that the CPU version is actually bugged and you won’t be able to create a star. During the relaxation of it from a MESA profile (as you will see in the tutorial “How to create a Star using MESA”), the Star will just explode if you use more than 6000 particles. Unless the bug will be solved, we suggest you to use the GPU code of the program, because it’s much faster and safe."

Second, we suggest Linux as the operating system (OS) to run StarSmasher. This doesn’t mean that you won’t be able to install StarSmasher in another OS: it's just not supported and could be challenging.

[//]: # "Just know that you can met difficulties and waste a lot of time if you are not skilled and experienced on installing libraries.
For example, one of us tried (for 2 months) to install StarSmasher in Windows via Cygwin with a very bad result. It proved impossible to install the GPU version of the code (only the CPU after lot of difficulties via the incompatibilities of the Cygwin’s libraries), and it was not possible to run StarSmasher because of some unresolved OpenMPI library errors. Instead, with Ubuntu the installation went perfectly fine in minutes. The advantage to use Ubuntu to run StarSmasher is that, if you miss some library during installation and get an error, the terminal will advice you that the 'XY (generic command) command is not been found'. Then if you are going to type 'XY' in the terminal you are going to know that the specific XY command is present in the Z package (a generic package) and you will be able to install it with the command"

[//]: # "```"
[//]: # "sudo apt-get install Z"
[//]: # "```"
[//]: # "Then, after the installation of the Z (generic package) package you will be able to continue the StarSmasher installation in a very clean way and without any problems. Instead, if you have windows or other OS, you can find problems, then you have to understand where the problem is, search where to get the specific command and other troubles that will only make you waste time unless you are skilled to do that kind of thing. This is why we suggest you to use Linux. If you don’t have it, just do a bipartition to your PC or machine. You just need 30-40 GBs free and then we suggest you to use an external Hard Disk with few Tera-Bytes of space to store your simulation."

[//]: # "If you are going to use a lot of particles (100k or more) each snapshot of your simulation is going to weight more than 30 Mb. Then if you will make it run for thousands and thousands of snapshots, well, you will need space. <!--- If you miss it StarSmasher is going to go in a weird error during the run. This is why w --> We suggest you make runs on your Ubuntu partition (that part of your computer HDD where Ubuntu is installed)<!--- to make practice and then install it directly on the Hard Disk-->.	
Be also sure, if you have a dipartition with Windows, to do a backup of your LINUX OS before installing Windows Updates. Why? Because it can happen that Windows Update will delete your LINUX dipartition, losing everything you installed. So we suggest you or to a backup or install LINUX on a totally different Hard Disk of your Pc. In this way you won't have any problem. We also suggest you an Hard Disk of some TBs of memory and make all your simulations there. THat's because StarSmasher's snapshot need space and you won't lose any data if something goes bad (in general)."

[//]: # "It’s impossible to run StarSmasher on VirtualBox. This is because StarSmasher need real CPU and GPU to run, not a virtual one. You will be able to install it on VB (only the CPU version) but if you’ll try to run it, mpirun will go in error. Then, please, do a bipartition!
We don't know if it's possible to run StarSmasher on a GPU cloud or to another Server with 8 GPU that you can rent. We will write there when we will try that. We hope of yes so you can run simulations with millions of particle without spending thousands of dollars/euros to build a very powerfull computer for only some definitives runs."

# [1] OBTAINING

## [1.1] With GitHub
Just go to the StarSmasher’s GitHub page and download the .zip file. Put it in the desired folder and extract it. You will get all the files. <!--- To travel between folders, you can simply open your terminal (Linux) in the specified folder. Once you have your files, you are ready for the installation.-->


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

[//]: # "--------------------------------------------------------"
[//]: # "PS: If you want to install the CPU version of the code, go to the paragraph [2.3]."
[//]: # "--------------------------------------------------------"

[//]: # "If you are running on Linux (Ubuntu especially), the first thing that you have to do is to install the correct driver for your Graphic card (even more than one if you have different graphic cards). To do that, type in your Ubuntu's terminal:"

[//]: # "```"
[//]: # "ubuntu-drivers devices"
[//]: # "```"
[//]: # "Your Ubuntu terminal will write (for my case):"

[//]: # " ```
== /sys/devices/pci0000:00/0000:00:01.0/0000:01:00.0 ==
modalias : pci:v000010DEd0000139Asv00001043sd00001C9Dbc03sc02i00
vendor   : NVIDIA Corporation
model    : GM107M [GeForce GTX 950M]
driver   : nvidia-driver-390 - distro non-free
driver   : nvidia-driver-418-server - distro non-free
driver   : nvidia-driver-450 - distro non-free
driver   : nvidia-driver-440-server - distro non-free
driver   : nvidia-driver-450-server - distro non-free
driver   : nvidia-driver-460 - distro non-free recommended
driver   : xserver-xorg-video-nouveau - distro free builtin"
[//]: # "```"
[//]: # "Then i will install the driver that is good for ME (460). So:"

[//]: # "```"
[//]: # "sudo apt install nvidia-driver-460"
[//]: # "```"
[//]: # "Obviously you have to install the driver that you need and that it's good for your NVIDIA Graphic card. It's reccomended to install the driver that Ubuntu suggest to you, the 460 in my case."

[//]: # "now reboot"

[//]: # "```"
[//]: # "sudo reboot"

[//]: # "```"

[//]: # "After rebooting, you have to install the NVIDIA nvcc toolkit, then just type this command in your Ubuntu terminal:"

[//]: # "```"
[//]: # "sudo apt install nvidia-cuda-toolkit"

[//]: # "```"
[//]: # "This last command is reccomended to install nvcc because we have already done a Makefile (MakefileGPUubuntu, MakefileGPUUbuntu10.04 and MakefileUbuntuCudaToolkit) that contain the path already written for this installation (as you will see later)."

[//]: # "But keep attention. If you have, for example, the latest NVIDIA graphic card or another OS, there is the possibility that NVIDIA developers didn’t updated 'nvidia-cuda-toolkit' to the latest version for your graphic card. There is the possibility, then, that you could get in trouble during installation. My personal suggestion is, if you won't be able to run the GPU version of StarSmasher following the installation via the last command, to"

Install cuda following the NVIDIA official [guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html), which is constantly updated.

[//]: # "Then, choose your Linux OS and follow step by step what the guide says to do. When you will choose your NVIDIA toolkit package, it's reccomended to follow the installation via deb(network) because it's the easier and it avoid many problems that could occur with other installation methods."

[//]: # "We also reccomend you to get an NVIDIA graphic card that has a high gravitation’s compute performance. A very new graphic card doesn’t mean that is less powerful to calculate gravity respect to the old one. You can check it [there](https://gpu.userbenchmark.com/Compare/Nvidia-Titan-RTX-vs-Nvidia-RTX-2080-Ti/m664199vs4027)"

[//]: # "As you can see the NVIDIA Titan RTX is newer than the NVIDIA RTX 2080 Ti, However, about gravity calculation, is just a 10% more performant, there is like no difference then. In this case the NVIDIA RTX 2080 Ti is a better choice because it’s older and there surely is the nvcc version ready for that graphic card."
Once your installation is complete, you now have a nvcc executable that will be used to compile the gravity library of StarSmasher.
[//]: # "Search it and keep in mind its location."


### [2.2] Compile the gravity library (CUDA) and StarSmasher

The SPHgrav_lib subdirectory contains code written by Evghenii Gaburov (and somewhat modified by Jamie Lombardi and Sam Knarr) for calculating softened gravitational forces and potentials on NVIDIA GPUs.
Now you must do some important changes in this. All depends on your OS. In the SPHgrav_lib there are few files including some makefiles. Every makefile is written for an OS. The Standard “makefile” is written for Keeneland. If you are using Ubuntu, you will use the makefile.ubuntu makefile, if you are using Quest, you will the makefile.quest and so on. You will find every makefile in the “misc” folder if there isn’t in the SPHgrav_lib. If there isn’t even there, you must substitute in the standard makefile the path to nvcc (as we told you some lines ago...). Once you have chosen your makefile (i did a makefile called "MakefileUbuntuCUDAToolkit" if you are following this guide spet by step), just change the name in “makefile” and delete the others, you don’t need them.	
In the makefile (Ubuntu case) there is written this string:

```
CUDAPATH       := /usr/lib/nvidia-cuda-toolkit/
```

In Ubuntu this is the path where the nvcc file is located if you installed it with the sudo apt-get command. Be sure that there is! If you installed the NVIDIA toolkit in another way, just go to the computer, find a file called "nvcc" and edit this path up to the folder "bin" where nvcc is located.

In the same file there is also this string to edit:

```
NVCCFLAGS := -arch=sm_61
```

As there is written, this string is for the NVIDIA GTX 1070. The number 61 indicate the computability version of the NVIDIA Graphic card. If you check there:

https://en.wikipedia.org/wiki/CUDA

the NVIDIA GTX 1070 has a computability version of 6.1. Then in the string there is written “NVCCFLAGS := -arch=sm_61”. If you have a different graphic card, this number must be change according to the computability version as the Wikipedia page says, or the gravity won’t be calculated and your stars will explode. Then, if you have a GeForce GTX 950M, like me, write:

```
NVCCFLAGS := -arch=sm_50
```

Because my graphic card is a 5.0, if you have, for example, a NVIDIA TITAN RTX, write:

```
NVCCFLAGS := -arch=sm_75
```

because it has a computability version of 7.5, and so on. This is why all of your graphic cards must be of the same computability version, or gravity won’t be calculated well and your stars will eventually explode.


Now that the Makefile of the SPHgrav_lib folder is ready, we can return to the "src" folder.

### [2.3] Setup your src's Makefile. - Go there if you want to install the CPU version of the code

Now that you SPHgrav_lib is complete (only if you are installing the GPU version of StarSmasher), go to the src folder. In this folder there are a tons of makefiles. As before, choose the correct one for you. For Ubuntu users that want to install StarSmasher with OpenMPI (this is the easiest way, you can also install it with IFort) there is an already written file called "MakefileGPUubuntu", if you want to use the GPU version (recommended), and "MakefileCPUubuntu" if you want to use the CPU version (unrecommended as explained earlier).  
But keep attention! For users that are using a version of Ubuntu 19.04 or older, they have to use the makefile called "MakefileGPUubuntu". If you are using Ubuntu 20.04 or newest, you have to use the makefile called "MakefileGPUubuntu20.04", or your compilation will fail with another makefile.
Then, if you have Ubuntu, just choose one of these three files and rename it as “Makefile”. 

Not makefile. Makefile. with the first capital letter M. Why? Because inside the makefile of the folder called "src" that you choosed (GPU version), there is written wich makefile to look. Usually the makefile to look is ever Makefile or Makefile.Ubuntu or Makefile.quest. Then be sure that the makefile inside SPHgrav_lib is called as the one written in the src's Makefile. So be careful, or the compilation will fail. To recap, if i'm running Ubuntu 20.04 I'll have to use the makefile called "MakefileUbuntu10.04". SO i will rename this makefile inside src as makefile or Makefile. Then, i'll go inside SPHgrav_lib to choose the makefile which will compile the GPU version of the code (MakefileGPUCudaToolkit for example). I must rename this makefile as Makefile. WIth capital M.

Be also sure that there are not other files called only “Makefile”. Delete the others makefiles that you don’t need if you want. If you are not an Ubuntu user, be sure to edit your makefile in way that it contains the installation path only for CPU/GPU and OpenMPI/IFort. I suggest you to follow the structure of "MakefileGPUubuntu" or "MakefileCPUubuntu" to do a clean EDIT.	
Once your Makefile is done, everything is ready. Open your terminal in your src folder and type:

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

If the CPU version of StarSmasher is been installed correctly, in your terminal will appear this phrase:

```
***MADE VERSION THAT DOES NOT NEED GPUS***
```

Now StarSmasher is ready to use! When you complete the installation, it will be created a .exe file called “test_gpu_sph if you installed the GPU version of the program, or “test_cpu_sph” if you installed the CPU version of the code. This .exe is situated in the upper folder of src. To run StarSmasher, you will have to move on that folder. Then just open your terminal where the .exe file is present or, after the installation, just type:

```
cd ..
```

Now, to run your first simulations, follow the tutorial “How to create a MESA star” situated in the “example_imput” folder!


# [3] Importants suggestions if you are running StarSmasher on your Hard Disk
A last important advice. When your simulation is finished and you want to transfer your files in another Hard Disk or similar, never cut and paste your files. StarSmasher’s files and snapshots are very sensible, and there is the risk that during the cut and paste process you may lose or damage your files. To avoid this, just copy and paste them, be sure that every files is there and that they aren't damaged (like visioning the texts or try to visualize them with SPLASH, as you will see next) and then, eventually, delete the original files once the transfer is been completed. 	
If you want instead to continue your simulation in another Hard Disk, be sure to recompile the code in your Hard Disk:

```
make clean
```
Then

```
make
```
Now the code is been recompiled and the mpirun command should run perfectly. 

We suggest to the user to run StarSmasher into an Hard Disk. If you are going to use a lot of particles (100000 or more) each snapshot of your simulation is going to be around 40 mb. Then, if you have thousands of them, there is the possibility that all your OS's space is going to be occupied. To avoid this, we suggest a NEW Hard Disk of 1 TeraByte or more. Why a New? There is not any known risk (actually) that StarSmasher is going to delete your personal files on your Hard Disk, but there is a problem with Linux and the file systems of the Hard Disk itself. If you are going to run StarSmasher in Ubuntu and try to install StarSmasher in an Hard Disk that has FAT32 file system, when you'll try "mpirun -np N test_gpu_sph", OpenMPI is going to say you that it is not possible to find the executable test_gpu_sph. This is because Linux won't give the permission to run the executable in a space that has not the same file system of Linux. However, everything will run perfectly if your Hard Disk is going to have NTFS's file system. Then, why a new Hard Disk? Because so you can format it with the correct file system of your linux OS and you won't lose any of your personal data into the Hard Disk, if you can't do a backup or for some other reason.

The solution of having a NTFS's Hard Disk is not universal and it may change in base of your OS or other factors. If your Linux OS is a FAT or FAT32, you'll need an Hard Disk with FAT or FAT32 file system and so on.


We wish to the user a good use of StarSmasher, actually best star collider simulator in the planet!
