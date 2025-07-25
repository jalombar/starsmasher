"Bleeding edge"
=================

As the name points out in a subtle way, this is the most updated 
version of the code. These notes are being developed as a stand-alone in the
`parallel-bleeding-edge` folder by Pau Amaro Seoane. Any mistake should be
attributed to me. Please do not use this version of starsmasher unless you are
willing to run intro trouble. If you spot any error, please let me know.

I am focusing on:

- Initial conditions created with MESA.
- A parabolic collision between two stars.

But StarSmasher is capable of many more things. Have a look at the upper
folders and the `README.md` file there.

WARNING
=======
Unless you know what you are doing and have a good overview of what StSm does and
how it works, *you should NOT be using bleeding edge*. You will be running into bugs
and problems, since this is meant to be a rolling version of StSm. Use the default
one instead.

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

A note on units
===============

StarSmasher assumes that `R_sun = M_sun = G = 1`

where

```
R_sun = 6.9599 * 10^10 cm
M_sun = 1.9891 * 10^33 g
G = 6.6739 * 10^-8 cm^3 g^-1 s^-2 # Newton's constant in cgs units
```

You can combine these three values to get the StarSmasher units for other quantities. The idea for any 
quantity is that its units are whatever combination of G, Msun, and Rsun give you that quantity.  
For example, the unit of density is Msun/Rsun^3.

A StarSmasher unit of time is

``` 
t_sun = sqrt(R_sun^3 / (G*M_sun)) = 1593.6 seconds
``` 

You can use the same approach for other physical quantities such as to get the pressure unit
that was also included in the output you provided, where

``` 
P_sun = G * M_sun^2 / R_sun^4 = 1.1253*10^16 g cm^-1 s^-2
```

A few more examples: You have to multiply the code output by these factors to
get physical units:

```
ConversionFactorLength   = 1           # To get lengths in Rsun
ConversionFactorLength   = 2.255e-8    # To get lengths in parsecs
ConversionFactorTime     = 1593.6      # To get time in seconds
ConversionFactorTime     = 5.053e-5    # To get time in years
ConversionFactorMass     = 1           # To get mass in Msun
ConversionFactorMass     = 1.9891e33   # To get mass in grams
ConversionFactorPressure = 1.1253e16   # To get g/(cm sec)
ConversionFactorDensity  = 4.1063e+11  # To get g/cm^2 (column density!)
```

Create 3D models from MESA 1D models
=====================================

I am going to assume that you have a directory in the main level named `MESA/`
with two subdirectories:

- `MESA_initial_1D_models/`
- `MESA_initial_3D_models/`

The folder `MESA_initial_1D_models/` is the usual working folder of MESA, i.e.

```
$ cp -r $MESA_DIR/star/work ./MESA_initial_1D_models/
$ cd ./MESA_initial_1D_models/
$ cp $MESA_DIR/star/defaults/profile_columns.list .
$ vim $MESA_DIR/star/defaults/profile_columns.list
```
These last two steps is to have information about the abundances, which might
be useful. Edit `profile_columns.list` and add just below `z_mass_fraction_metals`
the following:
```
   zone       ! numbers start with 1 at the surface
   mass       ! m/Msun. mass coordinate of outer boundary of cell.
   logR       ! log10(radius/Rsun) at outer boundary of zone
   logT       ! log10(temperature) at center of zone
   logRho     ! log10(density) at center of zone
   logP       ! log10(pressure) at center of zone
   x_mass_fraction_H
   y_mass_fraction_He
   z_mass_fraction_metals
   h1
   he3
   he4
   c12
   n14
   o16
   ne20
   mg24
   q
```
Save the file and modify `inlist_project` to make sure MESA reads in what you want:
```
&star_job
  ! see star/defaults/star_job.defaults

  ! begin with a pre-main sequence model
    create_pre_main_sequence_model = .true.

  ! save a model at the end of the run
    save_model_when_terminate = .false.
    !save_model_filename = '15M_at_TAMS.mod'

  ! display on-screen plots
    pgstar_flag = .true.

    profile_columns_file = 'profile_columns.list' ! < ---- HERE

/ ! end of star_job namelist

```
You might also be interested in having the evolutionary plots in pdf
format with a white background instead of the default one. For that,
change your pgplot configuration file to this
```
$ cat cat inlist_pgstar 
&pgstar
  ! see star/defaults/pgstar.defaults
    
    file_device = 'vcps'
    file_extension = 'eps'

  ! MESA uses PGPLOT for live plotting and gives the user a tremendous
  ! amount of control of the presentation of the information.

  ! show HR diagram
  ! this plots the history of L,Teff over many timesteps
    HR_win_flag = .true.

  ! set static plot bounds
    HR_logT_min = 3.5
    HR_logT_max = 4.0
    HR_logL_min = -1.0
    HR_logL_max = 4.0

  ! set window size (aspect_ratio = height/width)
    HR_win_width = 6
    HR_win_aspect_ratio = 1.0

  ! show temperature/density profile
  ! this plots the internal structure at single timestep
    TRho_Profile_file_flag = .true.
    TRho_Profile_file_dir = 'png'
    TRho_Profile_file_prefix = 'trho_profile_'
    TRho_Profile_file_interval = 5
    TRho_Profile_file_width = -1
    TRho_Profile_file_aspect_ratio = -1

  ! add legend explaining colors
    show_TRho_Profile_legend = .true.

  ! display numerical info about the star
    show_TRho_Profile_text_info = .true.

  ! set window size (aspect_ratio = height/width)
    TRho_Profile_win_width = 8
    TRho_Profile_win_aspect_ratio = 0.75


    TmaxRho_file_flag = .true.
    TmaxRho_file_dir = 'png'
    TmaxRho_file_prefix = 'trho_'
    TmaxRho_file_interval = 5
    TmaxRho_file_width = -1
    TmaxRho_file_aspect_ratio = -1

    logL_Teff_file_flag = .false.
    logL_Teff_file_dir = 'png'
    logL_Teff_file_prefix = 'logL_Teff_'
    logL_Teff_file_interval = 5
    logL_Teff_file_width = -1
    logL_Teff_file_aspect_ratio = -1

    L_Teff_file_flag = .false.
    L_Teff_file_dir = 'png'
    L_Teff_file_prefix = 'L_Teff_'
    L_Teff_file_interval = 5
    L_Teff_file_width = -1
    L_Teff_file_aspect_ratio = -1

    file_white_on_black_flag = .false.
    win_white_on_black_flag = .false.



/ ! end of pgstar namelist
```


Make and run
```
$ ./mk
$ ./rn
```

Use your 1D profile file created in a directory with MESA as input for the
`sph.input` (last rows). You will also need `sph.init`. You can find both
files in the `tools/` directory:

```
$ cp tools/sph.init_MESA  MESA_initial_3D_models/sph.init
$ cp tools/sph.input_MESA MESA_initial_3D_models/sph.input
```

After you have your profile from the 1D MESA file, say profile9.data, copy it
to the 3D folder to create an SPH profile which needs to relax before we start
the simulation,

```
   $ cp MESA_initial_1D_models/LOGS/profile9.data MESA_initial_3D_models/
   $ cd MESA_initial_3D_models/
   $ sed -i.bak "s|^profilefile=.*|profilefile='profile9.data'|" sph.input
```

Edit `sph.input` according to your needs and, then, run

```
   $ StSm_run.sh
```

and wait for it to finish. After that, take the last snapshot as your input 3D
file for the stellar collision. Rename it

```
   $ mv out0300.sph sph.start1u
```    

Repeat for the second profile, corresponding to the second star.

```
   $ mv out0300.sph sph.start2u
```

After that, you can clean the directory to save disc space,

```
   $ StSm_clean.sh
```

## Kinds of Collisions

### To initiate a collision

The three-letter code word in `sph.init` should be `hyp`, which directs the code to the initialization routine `initialize_hyperbolic.f`. In `sph.input`, you can use any of the following pairs of data to initialize the orbital parameters: `(e0, vinf2)`, `(e0, rp)`, `(semimajoraxis, rp)`, `(semimajoraxis, e0)`, or `(rp, vinf2)`. Here, `e0` is the eccentricity of the initial orbit, and `rp` is the periastron separation.

### Hyperbolic, parabolic or elliptical: `vinf2`

Whether the encounter is hyperbolic, parabolic, or elliptical can be controlled by the parameter `vinf2` (the relative velocity at infinity squared in code units) or by the parameter `semimajoraxis` (the semimajor axis). The code unit of velocity is approximately 436.73 km/s. The value of `vinf2` is negative for elliptical encounters, zero for parabolic encounters, and positive for hyperbolic encounters.

### Examples

For example, if you'd like a relative velocity at infinity of 1000 km/s, then you would set `vinf2` to approximately 5.243. As another example, if `vinf2` is set to approximately 11.796, then the relative velocity at infinity would be around 1500 km/s.

Note that the orbital energy (which relates to the semimajor axis) has a simple relationship with `vinf2`: the semimajor axis `a` is equal to `-G * M / vinf2`, where `M` is the total mass of the two stars. 

I.e. example, if you are colliding a 0.3 solar mass star with an 8 solar mass star, and `vinf2` is set to `-0.07`, then the semimajor axis `a` would be approximately 118.57, assuming the units are the default ones where `G`, the mass of the Sun, and the radius of the Sun are all equal to 1.


Brace for collision
====================

You can either (blindly) use the script `StSm_PrepColl.sh`

```
   $ StSm_PrepColl.sh
```

Or do it step by step. 

*(1) Create a collision directory,*

```
   $ mkdir Collision
```

*(2) Copy all necessary collision files* 

from `tools/` into that directory,

```
   $ cp tools/* Collision/
```

Move the 3D MESA SPH initial profiles into the Collision directory

```
   $ mv MESA_initial_3D_models/sph.start* Collision/
```

*(3) Start the simulation*

```
   $ cd Collision/
   $ StSm_run.sh
```

### Analysing the collision

You can use the scripts in the `bin/` folder. There's another `README.md` file
there.

Besides that, plotting quantities `energy?.sph` can help you to understand for instance
how many passages there are before merger. The columns in that file are as follows, with 
all quantities in code units (G=Msun=Rsun=1):

1. time t
2. potential energy W
3. kinetic energy T
4. internal energy U
5. total energy E =  W + T + U
6. total entropy
7. angular momentum

For certain types of runs, there will be additional columns as well.

For instance, ounting the number of major "dips" in W versus t gives the number of passages. 

But... what is a major dip from a passage, and what is a minor dip from an oscillation?

Dips in the potential energy W correspond either to a collision of the two components before 
final merging or to a maximum contraction during the subsequent oscillations of the merger remnant. 
The criterion we use to distinguish collisions (which should be included in the number of interactions 
`n_p` before the stars merge) from oscillations is that the first local maximum of W which is lower than 
the previous local maximum occurs immediately after the final merging. The idea behind this criterion 
is that a collision without merger ultimately tends to increase the system's gravitational potential 
energy, whereas a merger will decrease the potential energy.


Beyond the collision
=====================

The SPH collision will give you an outcome which is interesting in many
different ways, but maybe you are interested in the long term evolution of the
object, and that could take a long time to integrate, and a large accumulated
numerical error.

Fortunately, there is a tool to convert the three-dimensional output of
StarSmasher into a 1D file that MESA can read. This code has been developed by
Charles Gibson and can be found here:

<a href="https://github.com/charlie-gibson/sph-to-mesa/">https://github.com/charlie-gibson/sph-to-mesa/tree/main</a>

You can find a copy of this code in the `tools2/` directory. Please read the
detailed `README.md` file there.
