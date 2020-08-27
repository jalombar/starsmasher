Before reading this guide, we reccomend you to read the installation tutorial first and then "Creating a MESA star" tutorial to get confidence with the code and the commands.

## Walkthrough: fly-by of two stars

There are two steps to modelling a collision or fly-by of two stars.
First, you have to model the individual star: a relaxation run.
Second, you simulate the interaction: a dynamical calculation.

#### A note about Running the code
If you are running the code on a machine that uses a pbs scheduler, then I recommend that you do the following to prepare to run the code later on:
```
cd ~/starsmasher/parallel_bleeding_edge  ####  depending on where your starsmasher installation is located
cp ../misc/sph.pbs .
```

### Step 1: Relaxation runs

### Step 1(a): Relaxing a n=1.5 polytrope

Let's have you do your own simulation for practice.  We'll use N=1000, so that the code runs very quickly.  Let's try an equation of state
governed by the adiabatic index gamma=5/3, and a structure specified
by the polytropic index n=1.5.

First, get the directories ready:
```
cd ..  ####   or maybe "cd ~/starsmasher" depending on where your installation is
cp -r parallel_bleeding_edge GAM1.667_n1.5
cd GAM1.667_n1.5
cp ../example_input/relaxation_preMS/sph.in* .
```
The file sph.init controls what type of calculation we will do.
Confirm that the three letter code word inside sph.init is '1es' for this polytrope relaxation (and update the text if necessary):
```
emacs -nw sph.init
```
You should look at sph.input:
```
emacs -nw sph.input
```
Change line 4 of sph.input to be N=1000 and save the file.
The radius and mass of the polytrope can be set within the sph.input file using the variable names starradius and starmass, respectively.
For example "starmass=0.2" and "starradius=0.2"
These values are in physical units of runit and munit, which by default are the radius and mass of the sun in cgs units.
There are several other variables that could be changed or set, but the defaults should be fine for now.
See the Starcrash documentation for an explanation of what some of the other variables are or visualize the file "tips and tricks to improve your simulation" This file is going to be updated during time.

Now let's compile the source code:
```
cd src
make
```
Should you want to try a slow cpu-only version for test purposes, use "make cpu" instead.
Once everything goes well with the
compilation, the last line you should see is "mv GAM1.667_n1.5_sph
..".  Get positioned to run the code:
```
cd ..
```
If you are using a pbs scheduler, run the code using the following:
```
qsub sph.pbs 
```
Otherwise, use this command: 
```
mpirun -np 8 ./*_sph 
```
Where '8' can me changed to the number of parallel processes (often equal to the number of cores on your machine).

This creates, among other things, a file log0.sph that collects
ascii tex about the run.

### Step 1(b): Relaxing a TWIN model

Get the directories ready:
```
cd ~/starsmasher   ### or cd to wherever your main installation is
cp -r parallel_bleeding_edge TWIN
cd TWIN
cp ../example_input/relaxation_TWIN_model/* .
```

To relax a TWIN model instead, the three letter code word in sph.init should be 'erg'.
The SPH code will look for a model named 'eg.last1.muse_s2mm' by default.
You can create a symbolic link from an appropriate output file of TWIN if you like (e.g., "cp -s m_8_t_5.7.last1.muse_s2mm eg.last1.muse_s2mm").
The file m_8_t_5.7.last1.muse_s2mm is the default in example_input/relaxation_TWIN_model.
It represents an 8M_sun MS star at t=5.7Myr.

The sph.input file, also available within example_input/relaxation_TWIN_model/, should be similar to the sph.input file used above for the polytrope (probably just change N and nothing else unless necessary).
The file on the svn site does have TRELOFF to 200 in this case, to make sure there was enough time for oscillations to die off, and TF is 300 so the star could be monitored for 100 time units after the relaxation was turned off.

As discussed before, if you are using a pbs scheduler, run the code using the following:
```
qsub sph.pbs 
```
Otherwise, use this command: 
```
mpirun -np 8 ./*_sph 
```
Where '8' can me changed to the number of parallel processes (often equal to the number of cores on your machine).

### Step 2: Dynamical calculation

If the last model from the relaxation run looks good, then you can use the t=TRELOFF output file as an input to a dynamical calculation.
Start by get the directories ready:
```
cd ~/starsmasher   ### or cd to wherever your main installation is
cp -r parallel_bleeding_edge/ collision_rp3.9_a119
cd collision_rp3.9_a119
cp ../example_input/collision/sph.in* .
emacs -nw sph.init
```
Make sure the second line of sph.init says INAME='hyp'.

The new sph.input file is set up for doing dynamical calculations.
Take a look at it:
```
emacs -nw sph.input
```
While you're in the file, check that DTOUT=1 on line 3, which means that an out*.sph file will be dumped every 1 time unit.
Also, on line 16 or so, make BIMPACT=3.9d0, which means the periastron separation of the initial orbit is 3.9 stellar radii.

In sph.input you can use any of the following pairs of data to initialize a collision in hyperbolic.f: (e0,vinf2), (e0,bimpact), (semimajoraxis,bimpact), (semimajoraxis,e0), (bimpact,vinf2).  Here e0 is the eccentricity of the initial orbit, and bimpact is the periastron separation (not really the impact parameter).
For example, sph.input could contain

```
BIMPACT=4.d0,
semimajoraxis=118.57d0,
```

and e0 and vinf2 could be unspecified.
Or you could specify BIMPACT and e0, and the same code works, etc.
The code in initialize_hyperbolic.f figures out what it needs to solve for.

Whether the encounter is hyperbolic, parabolic, or elliptical is controlled by the velocity at infinity squared, vinf2, or the semimajoraxis a.
The value of vinf2 is negative for elliptical encounters, zero for parabolic encounters, and positive for hyperbolic encounters.
Note that the orbital energy G*M*mu/(-2a) just equals mu*vinf2/2.
Therefore, there is a simple relation between the semimahor axis a and vinf2, namely.  
a=-G*M/vinf2, where M is the total mass of the two stars.  
For example, if you are colliding a 0.3 Msun and 8Msun star and if vinf2=-0.07, then the semimajor axis a=(0.3+8)/0.07=118.57, where I assume the units chosen are the usual G=Msun=Rsun=1.

Note that ngravprocs=-6 and ppn=16... this means there are abs(ngravprocs)=6 GPU units per 16 CPU cores, as on one node of the supercomputer forge.
You should change these values to be appropriate for your machine.
For example, on Keeneland KFS, there are 3 gpus and 12 cores per node.

Let's get the star models:
```
cp ../TWIN/*/out0200.sph ./sph.start1u
cp ../GAM1.667_n1.5/*/out0200.sph ./sph.start2u
```
where we assume the 8 solar mass zero age main sequence star was relaxed in the directory ../TWIN/.
Now let's compile the code:
```
cd src
make
```
To run the code, get to the same directory as the executable:
```
cd ..
```
As discussed before, if you are using a pbs scheduler, run the code using the following:
```
qsub sph.pbs 
```
Otherwise, use this command: 
```
mpirun -np 8 ./*_sph 
```
Where '8' can me changed to the number of parallel processes (often equal to the number of cores on your machine).


## RESTARTING RUNS

Every few iterations, a restartrad.sph file is dumped (overwriting any previously existing restartrad.sph file).
If a restartrad.sph file exists in the directory when a new run is launched, then the code will automatically use that file to initiate the calculation.
This is useful for restarting in the middle of a long simulation.
You can also restart from any out*.sph file simply by renaming it to restartrad.sph.

