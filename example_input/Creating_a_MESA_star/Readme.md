# How to create a star in Starsmasher and use it to make a collision

To create a Star in Starsmasher and use it, at the moment the best way to do that is using MESA stellar evolution code.
MESA creates star models after you specify their mass and the metallicty.


## Installing MESA

If you are interested, there is a Youtube video(https://www.youtube.com/watch?v=NmaLHFxpALg) that explains how to install MESA. Here we summarize the main steps, but you should refer to the official documentation for the latest instructions and for more details.

First of all we have to install MESASDK, as it contains the libraries needed to install MESA
Download(http://www.astro.wisc.edu/~townsend/static.php?ref=mesasdk) it for your own OS.

Now, extract the folder in your HOME. It will create a directory called mesasdk-somenumbers, delete the upperfolder if there is one, and rename the direcory "mesasdk".

Now you can download MESA from the website(http://mesa.sourceforge.net/)

Once you have your MESA zip, extract it also in your HOME. It will create a folder called mesa-somenumbers. Just rename it as Mesa.
Now open the terminal in your HOME and type
```
ls -a
```
It will do a list of files, even the hidden files. We have to look at one file called .bashrc. We have to open it with:
```
gedit .bashrc
```
It will open a file with tons of line. Just go to the last lines and write these lines:


```
LINE TO WRITE                    COMMENTS

export MESA_DIR=/~Mesa           Here we have to write the path where the mesa folder is. In my case we are in the Home.
export OMP_NUM_THREADS=N        "N" is the number of CPUs of your machine per node. I've a quad-core, so i'll write 4.
export MESASDK_ROOT=~/mesasdk
source $MESASDK_ROOT/bin/mesasdk_init.sh
```

Now save the file and exit. CLOSE your terminal and open it again. In this way enviromentals paths will be applicated. Now, let's go to our Mesa directory:
```
cd Mesa
./install
```
Now many lines appear on screen. The whole process can take several minutes. Just take a coffe and come back when all is finished, which is when this phrase appears in your terminal:


```
mesa/astero has been built, tested, and exported.

************************************************







************************************************
************************************************
************************************************

MESA installation was successful

************************************************
************************************************
************************************************
```

# Create your star profile

Now that MESA is installed, we have to look at a folder called "work". You can find it in /Mesa/star/. The folder work is where we are going to create our star. Then copy this folder and put it where you want to do your simulation. As suggested in the tutorial isntallation, I suggest to put it in your hard disk.

Now let's look the work folder. Inside of it there is a file called inlist_project. Open it.

Here there are the options to create the star you want, the most important are as follows.

```
save_model_when_terminate = .false.         Set it as true
save_model_filename = '15M_at_TAMS.mod'     The name of the model
Zbase = 0.02
initial_z = 0.02                            These two should be the same. They controll the metallicity of the star.
initial_mass = 15 ! in Msun units           The mass of the star that you want
stop_near_zams = .true.                     This is when the simulation will end. If it's setted on true, then the star will stop to grow when it reaches the Main Sequence. If false, it will continue to grow until the simulation stop. If you prefer a precise age, edit the phrase:

create_pre_main_sequence_model = .true.
```
in
```
create_pre_main_sequence_model = .false.
```
And under it write:

```
max_age = 1.3d10    A number in years that you prefer, this is only an example
```

Now we are ready to start the simulation. Open the terminal in your work folder and type:

```
./mk
./rn
```

Now the simulation will start, and your star profiles will be created. 

Remember that you can run this folder where you want. If you want to make another different star, you can simply copy the folder from the original one situated in Mesa.

## Suggestion: Near the original "work" direcotry in Mesa, there is a directory called "test_suite". Here there are a lot of differents simulations that you can run to produce a variety of different star models. Have fun using these in StarSmasher!

# Transorm your MESA's star in a SPH star

At the end of the Mesa's simulation we have to look at the folder called LOGS and then at the files called "ProfileN.data". There will be many, each numbered "1", "2", "3", and so on. Each is a "snapshot" of your star. You have to look obviously at the profile you are interested in (that differ in radius, age, helium concentration and so on!). It is likely that the file that you want is the one with the biggest number, and the last.

Once you have your file, copy/cut it and paste it inside the folder where you want to do you simulation, besides the src folder and besides "sph.init" and "sph.input". You will have to edit your sph.input file. Be sure that, if your profile file is called, for example, "start1" (you can rename it as you want), Inside the sph.input you will have to write:

profilefile='start1'

Be sure then that the name of the profile and the name inside sph.input are the same!

Be sure that the nrelax=1. In this way Starsmasher will convert your mesa stars in SPH particles. The number N of SPH particles is set in the sph.input file.

After compiling the code and ensuring that all necessary files are in your folder (as shown in the example folder) and that everything is properly set in sph.input (such as the number of cores and GPUs per node on your machine), type the following command in your terminal:

```
mpirun -np N test_gpu_sph
```
Where N are the numbers of cores (cpu) of your machine.

At the end of the simulation, you will have many files called out0000.sph, out0001.sph and so on. The last file should be a "out0300.sph"; it obviously depend when you decide to finish your simulation (specified in tf=N inside sph.input).

# Use the star in a collision

To use your star in a simulation, create a new folder, copying the "collision" folder that is in "example_imput" folder. In that folder you have to delete the files called "sph.start1u and sph.start2u" replacing them with the last files of your relaxation runs. These files contain you star models!

If you want that your star is the target of the collision, rename it "sph.start1u". If you want that is the impactor, call it "sph.start2u. Often, you want two differents stars there.  Alternatively, you can use two identical models: if so, copy the last output of the relaxation and paste it two times. Call one sph.start1u and the other sph.start2u.

At this point set the collision parameters in the sph.input for the that collision that you want.

Then start the simulation:

```
mpirun -np N test_gpu_sph
```

Where N is the number of cores (cpu) of your machine.

The simulation will create many out*.sph files. You can see them following the tutorial inside the "visualizing data" folder of the Starsmasher repository!

In this folder there are also 2 already done MESA stars!    Start1 is a SUN (1 solar mass and 1 solar radius), Start2 is a red dwarf (0.5 solar mass and 0.5 solar radius). You can use them to practice with Starsmasher! 

There is also another file in the previous folder called "tips and tricks to improve your simulation". Here is a list of tricks that you can add inside your sph.imput before relaxing the star, for example, to make it rotate!

For problems please contact Francesco Radica at "francyrad.info@gmail.com".

