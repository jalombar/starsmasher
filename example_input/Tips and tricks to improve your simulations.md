sph.input can be edited to add features for you stars, for example:

# How to make a rotating star 


You can make your star rotate. To do that, add inside your sph.input file the line:
omega_spin=N

Where N is the angular velocity of your star. To decide the correct number to put, we have to do a little calculation.
If you want, for example, that you star do exactly ONE rotation (2π) in ONE day of simulation (54 times of unit in Starsmasher,
that equals to 1 day)
N will be:

2π/54 =0.116             then  

omega_spin=0.116


If you want 2 rotations in 1 day, then the equation is 4π/54=0.29 ;     then omega_spin=0.29

If you want 1 rotations in 2 days the the equation is π/54=0.058;      then omega_spin=0.058

At the end just relaxe your star, after 200 time units it's ready to use it in a collision! If you notice, your star during the relaxation
is rotating!

# How to make High framerate simulations

Every simulation will release a snapshot. If your simulation looks laggy, that’s because the simulation releases a snap file at a low frame rate. 
To solve this problem, just open your sph.input file. There is a string called “dtout=1”. This means that every 1 t of simulation, StarSmasher will release an out*.sph,
however this means a laggy visualization. To solve, do that. Change this number with 0.1 and also change the tf=N to tf=200. This means that StarSmasher will release
a Snapshot every 0.1 t, then 10 snapshot every t. Now your simulation will have a high frame rate during visualization with SPLASH, however you will likely have 2000
snapshots files after 200 t. Now that your stars smashed and you want to continue to visualize how the collision evolve, you don’t need any more a high frame rate. 
Then, once your simulation is finished, for example at a t=200, or whatever you want, just change “tf=N” (until you want to continue, obviously a number higher than 
200) and dtout=1 or more. Now your simulation will release 1 snapshot every 1 t, and the snapshot numeration is going to continue (then if it stopped to 2000, it will
continue with 2001, 2002 and so on).






----------------------------

#Computeexclusivemode; what is it and what does it do?

The safest thing to do is not to set computeexclusivemode in sph.input so that it's assigned its default value of 0.  (computeexclusivemode=0).
If a user is going to set it, they should set it to 0 unless (a) they want more than one job is to be run on the same compute node and 
(b) the GPUs have been configured by the system administrator to be in exclusive mode.  
For example, if a compute node has 4 GPUs and the user wants one job to use two of the GPUs and another job to use the other two GPUs, then, depending on how 
the hardware was configured, it might be appropriate to set computeexclusivemode to 1.  By setting computeexclusivemode to 1 (computeexclusivemode=1), the user
is telling StarSmasher that the GPUs have been configure such that only one COMPUTE thread is allowed to run on each GPU.  
Setting the value to 1 does *not* actually change the configuration of the GPUs (that would be done by the system administrator with commands 
such as "nvidia-smi --id=0 --compute-mode=EXCLUSIVE_PROCESS".  If the user sets computeexclusivemode equal to 1 when the GPUs are in their default configuration, then 
the code will run slower because a single GPU will do all of the calculations which could have been spread out over multiple GPUs

Then, if you want a very fast simulation, computeexclusivemode=1 in sph.input
## THIS TUTORIAL WILL BE UPDATED DURING TIME

For any questions contact Francesco Radica at:  francyrad.info@gmail.com
