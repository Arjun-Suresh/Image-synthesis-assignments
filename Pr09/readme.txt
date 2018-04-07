The cpp files runs using c++ and employs freeglut for creating a window. In order to set the environment for running this app, please follow these steps for linux:
1. sudo apt-get update
2. sudo apt-get install libglu1-mesa-dev
3. sudo apt-get install libxi-dev
4. Download freeglut package for the appropriate OS version (I used Ubuntu 16.0.4 and took the package from https://launchpad.net/ubuntu/+archive/primary/+files/freeglut_2.8.1.orig.tar.gz)
5. Extract the files from the package and cd to the extracted folder in terminal.
6. Run ./configure, make and sudo make install commands on the extracted folder.
7. The compiled binaries need to be copied to /usr/lib. Run cp /usr/local/lib/* /usr/lib

Now, the environment is set up in the system and the source code can be compiled. Each application creates an output ppm files in the same location as the output file. In my submission, I am attaching this readme and 2 cpp files. The 2 description of the 2 files and command to compile each is given below:
1. pr04.cpp - This file contains the code for point, area, spot and direction light illumination and shadows.
   To compile, please run: g++  -o pr04 pr04.cpp -lglut -lGL -lm
2. subSurface.cpp - This file contains the code to get shadows using d/r method.
   To compile, please run: g++  -o subSurface subSurface.cpp -lglut -lGL -lm

These 2 commands create the output files which can then be run seperately.
