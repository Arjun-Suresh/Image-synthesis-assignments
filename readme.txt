The cpp files runs using c++ and employs freeglut for creating a window. In order to set the environment for running this app, please follow these steps for linux:
1. sudo apt-get update
2. sudo apt-get install libglu1-mesa-dev
3. sudo apt-get install libxi-dev
4. Download freeglut package for the appropriate OS version (I used Ubuntu 16.0.4 and took the package from https://launchpad.net/ubuntu/+archive/primary/+files/freeglut_2.8.1.orig.tar.gz)
5. Extract the files from the package and cd to the extracted folder in terminal.
6. Run ./configure, make and sudo make install commands on the extracted folder.
7. The compiled binaries need to be copied to /usr/lib. Run cp /usr/local/lib/* /usr/lib

Now, the environment is set up in the system and the source code can be compiled. Each application creates an output ppm files in the same location as the output file.
Each of the individual folders contain different projects from basic camera rendering to illumination, shadows etc.
