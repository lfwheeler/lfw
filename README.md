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

Inside the src folder, there is a warpt folder which contains our main research. The wlw folder in the warpt folder contains files for current research and that were used to make images for the CWP sponsor's report and for slides presented at the 2014 CWP sponsor's meeting. The data folder contains contains data for five velocity logs (vlogs.txt), 11 porosity logs (plogs.txt), and 13 density logs (dlogs.txt) from Teapot Dome, provided by Rocky Mountain Oilfield Test Center, to be used with wlwDemo.py. A description of how the file was made is included in the first few lines of the file. 
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

A demo for using WellLogWarping.java is included in wlwDemo.py. Data to be used with this demo can be found in the data folder. 

If you have data that you want to run through the existing scripts, modify the methods that deal with extracting the data or feel free to create your own method.

If you have any questions, please do not hesitate to contact me at lwheeler@mines.edu.
########################################################################################
