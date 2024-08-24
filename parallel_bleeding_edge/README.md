Note
=====

These notes are being developed as a stand-alone in the
`parallel-bleeding-edge` folder by Pau Amaro Seoane. Any mistake should be
attributed to me. If you spot any, please let me know.

I am focusing on:

- Initial conditions created with MESA.
- A parabolic collision between two stars.

But StarSmasher is capable of many more things. Have a look at the upper
folders and the `README.md` file there.


Compile 
========

Compile the source and place the executable in the `bin/` directory, after you
have modified the makefile for your system. 

*As an important advice:* 

Please make sure that you do not have more than one `cuda` library in your system. If you do, then make sure that the makefile in the `SPHgrav_lib2` and in the main `src` directories are (1) using the same one and (2) that you choose the more recent release.


```
   $ cd ./src/SPHgrav_lib2
   $ make
   $ cd ..
   $ make
   $ mv ./*_sph ../bin/
```

*Another important note:*

You might run into difficulties when trying to set up the `Makefile` in the `SPHgrav_lib2`
folder. My advice is that you

First find out what CUDA card you have. The best is to run `$ lspci | grep VGA`

```
$ lspci | grep VGA
45:00.0 VGA compatible controller: NVIDIA Corporation TU102GL [Quadro
RTX 6000/8000] (rev a1)
61:00.0 VGA compatible controller: NVIDIA Corporation TU102GL [Quadro
RTX 6000/8000] (rev a1)
```

In my case, I have two cards, `TU102GL`. Thanks to the wikipedia, 

<a href="https://en.wikipedia.org/wiki/CUDA">https://en.wikipedia.org/wiki/CUDA</a>

(look for "Compute Capability, GPU semiconductors and Nvidia GPU board products"),
I know that I have to use 

```
NVCCFLAGS := -arch=sm_75
```

In my case, I had to remove the `-abi=no` flag from the makefile to make it compile.

Adding the `bin/`directory to your path: Why you should do it
=============================================================

In `bin/` you will find not just the executable, but other interesting scripts
useful for analysing the data, cleaning, etc.  All of them start with `StSm_`
so that you can easily find them from the command line.

You should add the `bin/`directory to your shell path. In `zsh` (but also
`bash`) the syntax is

```
export PATH="$PATH:/bin:/usr/bin:/usr/local/bin:/sbin:/home/pau/StarSmasher/parallel-bleeding-edge/bin/"
```

This way you can start a simulation anywhere you want and analyse the datafiles, etc.


Create 3D models from MESA 1D models
=====================================

Get into MESA/ and use your 1D profile file created with MESA
as input for the sph.input (last rows). This will create an SPH
profile from the 1D MESA file, say profile9.data. It needs to relax.

```
   $ cp MESA_profiles/LOGS/profile9.data ./MESA/
   $ sed -i.bak "s|^profilefile=.*|profilefile='profile9.data'|" sph.input
   $ StSm_run.sh
```

Then wait for it to finish.

Take the last snapshot as your input 3D file for the
stellar collision. Rename it

```
   $ mv out0300.sph sph.start1u
```    

Repeat for the second profile, corresponding to the second star.

```
   $ mv out0300.sph sph.start2u
```

Brace for the collision
========================

Use the script `StSm_PrepColl.sh`

```
   $ StSm_PrepColl.sh
```

Or do it step by step. Create a collision directory,

```
   $ mkdir Collision
```

Copy all necessary collision files from `tools/` into that directory,

```
   $ cp tools/* Collision
```


Move the 3D MESA SPH initial profiles into the Collision directory

```
   $ mv MESA/sph.start* Collision/
```


Get into the directory and start the simulation,

```
   $ cd Collision/
   $ StSm_run.sh
```
