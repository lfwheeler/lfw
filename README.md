Simultaneous correlation of multiple well logs code
---------------------------------------------------


#########Software needed to run our code################################################
To use this code, you will need to install jython (to use our .py scripts) and the java jdk (to run the java code). This code will also require the use of the mines jtk.

jython can be downloaded from http://www.jython.org/.
The mines jtk can be downloaded from https://github.com/dhale/jtk. Please read the corresponding readme.txt file.
Java SE JDK 7, which is the version used when writing the .java files, can be downloaded from http://www.oracle.com/technetwork/java/javase/downloads
########################################################################################

#########Organization of the folders####################################################
In the lfw folder is my workbench (bench folder) for research. Inside the bench folder,
I have the corresponding src folder (for our source code), the build folder 
contains the .class files corresponding to the .java files, and
the other folders and files are for building the .java files with gradle.

Inside the src folder, there is a tpt folder and a warpt folder. The tpt folder contains files for importing well logs, seismic, and horizons. You only need WellLog.java for well log correlation. The warpt folder contains our main research.

The wlw folder in the warpt folder contains files that were used to make images for the CWP sponsor's report and for slides presented at the 2014 CWP sponsor's meeting.
########################################################################################

#########Building our code##############################################################
To build the .java files for the warping with wavelets code, refer to “Building the Mines
JTK” section in the readme.txt file in the mines jtk.
For a Mac OS, the required folders and files to build the .java files are located in our bench folder.
If you are running a Mac OS, you will need to cd into the bench directory and type “sh gradlew” to build the .java files. If you
are using another operating system, please refer to “Building the Mines
JTK” section in the readme.txt file in the mines jtk.
########################################################################################

#########Research#######################################################################
Inside the warpt folder, you will find files pertaining to our research.

In each file, there should be a small description of the purpose of why the file exists
in the first few lines of the file.

The wlwDemo.py requires data from Teapot Dome...

If you have data that you want to run through the existing scripts, modify the methods that deal with extracting the data or feel free to create your own method.

If you have any questions, please do not hesitate to contact me at lwheeler@mines.edu.
########################################################################################
