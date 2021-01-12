# How to see the snapshots of your simulation

Every simulation will create many files called "out0000.sph; out0001.sph; out00002.sph" and so on.
Each one of these files is a snapshot of your simulation. To visualize them we suggest you to use the program "SPLASH" of Daniel Price: 

https://github.com/danieljprice/splash


# Install SPLASH (new way)
The upper web page has a full documentation about the installation of SPLASH, however we want to do there a little guide that is setted up just for StarSmasher.

First of all, download Splash 3.0 from the github repository and extract it to your HOME directory. It must be there, not elsewhere, as you will see later.

Now go to the splash folder, clicking on it or with:

```
cd splash
```
Inside there will be many folder, but the program is not complete! You have to add another one called "Giza". To do so, inside the splash folder, type:

```
git clone https://github.com/danieljprice/giza.git
```
This will create a folder called Giza with many files inside it. Now go to the Giza's direcotry clcking on it or typing:

```
cd giza
```
Now that you are inside the giza's folder, type in your terminal:

```
./configure
```

It will appear many lines of installations

Then type:

```
make
```

still lines; then

```
sudo make install
```

It can happen that the installation of giza ends in errors, but don't worry, it's normal.

Now go back to the splash folder:

```
cd ..
```

We have to install another library now called "Cairo". To do so, we have to use a script inside the folder, then type in your terminal:

```
./install-cairo.sh
```

This will download the library and it could take some time. Once done, we are ready to install Splash. Now type:

```
make SYSTEM=gfortran withgiza
```

If everything went well, after many lines it will say you that SPLASH is been installed and to complete the installation you have to type the command "sudo make install"; then type it!

```
sudo make install
```
------------------------------------------------------------------------------------------------------
1S POSSIBLE MISHAP

Depending on your Machine/OS, this installation can go in error because it couldn't find the path /usr/local/bin of you computer. If that happen, type this command to solve the problem:

```
sudo mkdir usr/local/bin
```
It will create the bin folder, now type again

```
sudo make install
```

Now it should go fine.

---------------------------------------------------------------------------------------------
To end the installation, you have to add some enviromental variables to make splash avaiable every time that you open the terminal. To do so, open the terminal in your HOME and type:

```
ls -a
```
This will list many files including the hidden file. We have to look to a file called .bashrc. We have to open it. Now type:

```
gedit .bashrc
```

This will open thE hidden file. Inside the file, you will see many lines. Go to the end of the file, add some spaces to separate the last lines from what we are going to write and copy and past these 3:

```
export SPLASH_DIR=$HOME/splash
export PATH=$PATH:$SPLASH_DIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$SPLASH_DIR/giza/lib
```

Now save the file. 
If everything went well, close your terminal and open it again. Now type:

```
splash
```

If everything went fine, the Splash home should appear. 

```
      _(_)     ___ _ __ | | __ _ ___| |__   (_)     _ (_)
   _ (_)  _   / __| '_ \| |/ _` / __| '_ \      _  (_)   
  (_)  _ (_)  \__ \ |_) | | (_| \__ \ | | |  _ (_)  _    
      (_)  _  |___/ .__/|_|\__,_|___/_| |_| (_)  _ (_)   
          (_)  (_)|_| (_) (_)  (_)(_) (_)(_) (_)(_)      

  ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e )

  ( v3.0.1 [27th Aug 2020] Copyright (C) 2005-2020 )

```
--------------------------------------------
2ND POSSIBLE MISHAP

If instead to appear the splash home it gives to you this error:
```
splash: error while loading shared libraries: libgiza.so.0: cannot open shared object file: No such file or directory
```
add this enviromental variable to your bash file

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
```

PS: This path is been setted because the file libgiza.so.0 is in /usr/local/lib most of the times (like in Ubuntu), but it can change depending on your machine.

We suggest you, then, to find this file first and then edit the path xxx/yyy/zzz depending where your giza file is.

Close your terminal and open it again. Now Splash should work.

----------------------------------------------------------------------------

The command to read the Starsmasher's snapshots is:

```
splash -f starsmasher out*.sph
```
However this command is too long and it may be annoying to write it every single time to read our data. Before the update 3.0 of Splash, Starsmasher's user used to use the command "jsplash" (from jamie splash). It's quick to write and it will permit you to check the datas fast during your simulation. To implement the command jsplash along "splash -f starsmasher", open again your .bashcr as you did to write the paths during the installation and, at the end of the file, write:

```
alias jsplash='splash -f starsmasher'
```
Save .bashrc now. Close your terminal and open it again. now type:

```
jsplash
```
If everything went well, it should now appear the splash's home and you'll be able to use the easier command "jsplash" to read the data of StarSmasher.


# Install SPLASH (old way)
You shouldn't need this nowday, but we keep the files there to prevent any inconvenience
With the old way, you have to install the version 2.5.0 of the 1st february 2020, but with some changes. When you are going to install SPLASH (as described in the last paragraph), you must substitute the "Makefile" of the SPLASH folder called "build" with the Makefile of this folder. You have also to add the file "read_data_jamiesph.f90" in the folder of SPLASH repository called "src".

If everything went fine, it should be added a command called "jsplash" to your terminal.


# How to use SPLASH to read Starsmasher's files.

Be also sure to have the files "splash.defaults & splash.units" in the folder of your simualtions where the files outputNNNN.sph are present. These files are needed to visualize the correct space units in solar radius. Then just copy and past "splash.defaults & splash.units" in your simulation's folder. You can find them on the folder called "misc". Be also sure to open the terminal in your folder's simulation obviously.

At the end, to read your data. just type in your terminal:
```
"jsplash out0000.sph" 
```
If everything went well, SPLASH will ask you to choose the axis of the simulation to visualize in the y and x axis of SPLASH, it will also ask you if you want to visualize the density, the temperature or the mass particle. For a standard visualization, just type 2 (y axis), 1 (x axis), 8 (density) and 0 (no speed vector). 
Then type:

```
/xw
```

 It will appear you the snapshot of the simulation. Typing "h" in the virtual window make appears a list of commands in your terminal that you can use to visualize better your snapshot. For example you can zoom, scale down, rotate the picture and so on. then press "s" to save the current visual position of your snapshot and then q to close it.  
At this point press again, for a standard visualization: 2, 1, 8, 0. At this point, to get your picture, instead of "/xw" type:

```
/NAMEYOUWANT.png
```
In this way SPLASH will create a picture on your folder.

If you want to visualize all the data, just type:

```
jsplash out*.sph
```

This command will read you all the files. When you will type "/xw", just press the space bar to visualize the next snapshots and "b" the previous one. If you want to quit, just press "q; esc or ctrl+c". We suggest also to read all the Daniel Price's guidelines for insights and tricks about visualizing datas and creating other kinds of pictures or videos:    

https://splash-viz.readthedocs.io/en/latest/


For any problems or question during the installation, visualization and other basic functions, contact Francesco Radica at the email "francyrad.info@gmail.com"


Edit:   If you are an Ubuntu user, there is this command to install splash

```
sudo apt-get install splash
```
However it stil doesn't have the command to read the Starsmasher's file. This is because you need to install directly the repository of SPLASH. When it will be avaiable, we will move this command upward and left the rest of the guide for other OSs users.
