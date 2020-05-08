# How to see the snapshots of your simulation

Every simulation will create many files called "out0000.sph; out0001.sph; out00002.sph" and so on.
Each one of these files is a snapshot of your simulation. To visualize them we suggest you to use the program "SPLASH" of Daniel Price: 

https://github.com/danieljprice/splash


YOU MUST install the version 2.5.0 of the 1st february 2020, but with some changes. When you are going to install SPLASH, you must substitute the "Makefile" of the SPLASH folder called "build" with the Makefile of this folder. You have also to add the file "read_data_jamiesph.f90" in the folder of SPLASH repository called "src".

Be sure to follow the Daniel Price's guide to install SPLASH. If everything went fine, it should be added a command called "jsplash" to your terminal.

Be also sure to have the files "splash.defaults & splash.units" in the folder of your simualtions where the files outputNNNN.sph are present. These files are needed to visualize the correct space units in solar radius. Then just copy and past "splash.defaults & splash.units" in your simulation's folder. Be also sure to open the terminal in your folder's simulation obviously.

At the end, to read your data. just type "jsplash out0000.sph" in your terminal. If everything went well, SPLASH will ask you to choose the axis of the simulation to visualize in the y and x axis of SPLASH, it will also ask you if you want to visualize the density, the temperature or the mass particle. For a standard visualization, just type 2 (y axis), 1 (x axis), 8 (density) and 0 (no speed vector). Then type "/xw". It will appear you the snapshot of the simulation. Typing "h" in the virtual window make appears a list of commands in your terminal that you can use to visualize better your snapshot. For example you can zoom, scale down, rotate the picture and so on. then press "s" to save the current visual position of your snapshot and then q to close it.  
At this point press again, for a standard visualization: 2, 1, 8, 0. At this point, to get your picture, instead of "/xw" type /NAME.png. In this way SPLASH will create a picture on your folder.

If you want to visualize all the data, just type "jsplash out*.sph". This command will read you all the files. When you will type "/xw", just press the space bar to visualize the next snapshots and "b" the previous one. If you want to quit, just press "q; esc or ctrl+c". We suggest also to read all the Daniel Price's guidelines for insights and tricks about visualizing datas and creating other kinds of pictures or videos:    

http://users.monash.edu.au/~dprice/splash/userguide/index.html


For any problems or question, contact Francesco Radica at the email "francyrad97@gmail.com"


Edit:   Daniel Price is working to implement the command "jsplash" in splash. In this way it will be already present when downloading splash in the classical way with the command:

sudo apt-get install splash


This mean that you won't need anymore in the future to install splash from the github repository (in Ubuntu).
