# How to install splot routines

What are splot routines? And what is their utility in StarSmasher?
Splot routines are tiny programs useful to calculate data that are not available with StarSmasher alone. For example, the simulation alone cannot say you how much mass it’s been loss by both star during a collision, or it cannot say you how the eccentricity of the orbit changed during time… and so on. However, splot routines can! Splot routines are able to calculate these kinds of data writing all you need in a file that you can use to make plots with SPLASH. They are essentials to make research and write a paper! 
How do they work? How do you install it? The procedure is very simple! In the folder where this “Readme” is located there are other folders where splot routines are present. Each folder contains a splot routine (files and a Makefile). To install the splot routine and make the program, just open your terminal in the folder and type:

```
make
```

Be sure that ALL the folders of the path where the splot routine is located doesn’t have spaces between the name. For example, your folder called “splotroutine”, where the files are present, cannot have this path:
Pc: Desktop/StarSmasher/splot routines test/splotroutine
Because the compilation will fail. To solve, the folders must not have spaces in the name. Then you can adjust your folder path, according to the preceding example, in this way:
Pc: Desktop/StarSmasher/splot_routines_test/splotroutine
Now “make” should work perfectly and it should create a program called “splotroutine”. The program name will ever be like the name of the folder. Now, how to use it?
To use the program, just open the terminal where the program is, and type (according to the last example):
./splotroutine
The terminal will open it and it will ask you to type 2 numbers in this case: 0 or 71. Where 0 is to exit, 71 is to execute the program.
To read your snapshots files of your simulation, you must have all your out*.sph file in the folder where the splot routine is. So just cut and paste all your file in the folder of the splot routine or, better, just install it where your snapshots files have been created during the simulation. I suggest this last one because it’s safer, easier and quicker.
Then, once your out*.sph files are present, type in the terminal “71”; it will ask you:
(1) The number of the first output to read, (2) the last number of the output file to read, (3) the step to read the output files.
For example, if your simulation is composed of 501 snapshot files, where the first is out0000.sph and the last is out0500.sph and you want to read each of them, you have to type:
0000 0500 1
Then press enter. The program will start to compile a file called “MassAndMore.out”.     The file will contain the data that you need, calculated time by time.
To plot your data, type in the terminal:
splash MassAndMore.out     
SPLASH will then show you a list of numbers, each number correspond to a number related to the column of your file (for example “1” is time, “2” is the mass lost form the first star and so on; according to the standard splot routine) and will ask you what do you want to plot in the y axis and in the x axis. Type your desired numbers and then check the graph typing	

```
/xw
```

SPLASH will show you the graph. Adjust it as you want (you will find all the instructions on how to use SPLASH in the SPLASH’s user manual), when you like it press:

```
s
```

Then:

```
q
```

To exit.

 Then again. Plot the y axis as before, the x axis as before and then write:
 
 ```
NameYouWant.png
```

SPLASH will create your picture ready to use for your paper!



You can also edit the numbers as you want with further post processing with Excel or another editor, just remeber to insert columns between the numbers!

# Other post processing files
StarSmasher is able to produce some files that contains useful data. 

## Energyx.sph
This file is described in details in the documentation.
The file energy.sph lists:
```
column 1: the current time,
column 2: the total system potential energy, 
column 3: kinetic energy, 
column 4: internal energy, 
column 5: total energy, 
column 6: total entropy, 
column 7: total angular momentum, 
column 8:the maximum SPH density for any particle.
```

## parent.sph
That file is generated during the initialization of a relaxation run by initialize_parent.f.  The file isn't used for the operation of the code but rather is provided in case the user wants to check the splined profiles.  Unless otherwise stated, quantities are in code units.
This file contains:
```
Column 1: radius of shell in solar radius,
Column 2: pressure at that shell,
Column 3: density at that shell,
Column 4: temperature of that shell in Kelvin,
Column 5: mean molecular weight of that shell in grams,
Column 6: specific internal energy of that shell.
```
## Cool****.sph
These files are created in output.f and have much of the same data as in the out*.sph files but are in ascii format.  They are not really need but can be convenient for someone who wants to check if the relaxation is going well.  The particle data within them could be compared against the profiles in the parent.sph file to make sure the SPH model is close to the desired model.  Each line in the col*.sph files corresponds to data about one particle.  The quantities are in code units unless otherwise stated:
```
Column 1 is radius
Column 2 is pressure
Column 3 is density
Column 4 is temperature in Kelvin
Column 5 is mean molecular mass in grams
Column 6 is particle mass
Column 7 is smoothing length
Column 8 is neighbor number
Column 9 is the radial component of the gravitational acceleration
Column 10 is the radial component of the hydrodynamic acceleration
Column 11 is the x coordinate
Column 12 is the y coordinate
Column 13 is the z coordinate
Column 14 is the gravitational potential
Column 15 is the specific internal energy
Column 16 is the velocity squared
```

# Unit of measure
Those quantities  are in code units: G=M=R=1, where the default M and R are the mass and radius of the Sun.   You can change the mass and radius unit by setting values of munit and ruint in sph.input.  The usual thing to do is not set these values and then the default values, set in init.f, will be used:

```
runit=6.9599d10          ! number of cm in the unit of length.  use 6.9599d10 if want solar radius.
munit=1.9891d33          ! number of g in unit of mass.  use 1.9891d33 if want solar mass.
```

So, anyway, let's say we're using the default values where M=Msun and R=Rsun.  So all masses are in units of Msun and all distances are in units of Rsun.  Then The unit of density is Msun/Rsun^3=5.8999 g/cm^3.  The unit of time is sqrt(Rsun^3/(GMsun))=0.018445 days. The unit of pressure is G*Msun^2/Rsun^4.  The unit of energy is G*Msun^2/Rsun.  The unit of specific energy (like the specific internal energy in parent.sph for example) is GMsun/Rsun.  Etc.  Basically you take whatever combination of G, Msun, and Rsun that gives the dimensions of the quantity in question

# This  tutorial is in continuous updating. We are collecting all the splot routines around to make all them available! 
