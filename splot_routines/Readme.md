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
1The number of the first output to read, (2) the last number of the output file to read, (3) the step to read the output files.
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

a last suggestion: To have a better resolution data, we suggest to make the first snapshot of the simulations with an output rate of 0.1 t. Then stop the simulation when all the major events are done and then change it to 1.

# This  tutorial is in continuous updating. We are collecting all the splot routines around to make all them available! 
