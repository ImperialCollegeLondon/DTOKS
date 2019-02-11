\*! \mainpage DTOKSU
\section intro_sec Introduction
Author :	Luke Simons\n
Date :	06/02/2019\n
email:	ls5115@ic.ac.uk\n\n
Description :\n
This file describes the class structure and operation of DTOKSU code, \n
available at https://github.com/LukeSimons/DTOKS-U.\n
Tests exist in local directory '/Tests' with a readme file available at \n
'/Tests/README_TESTING'\n
\n
At Imperial College, London the code Dust in TOKamaks (DTOKS) was developed \n
in 2005 by James D. Martin, Michael Coppins and John Allen. Since, it has \n
been further developed by Minas Bacharis and others. DTOKS represents a \n
relatively independant and unique approach to computing dust tracks in \n
tokamaks and has hadsome implementation modelling dust in MAST, ITER, JET \n
and radiofrequency discharges.\n
\n
For more information on the physical models implemented, see the Thesis of \n
Dr. James Martin, and various papers:\n
J. D. Martin, M. Coppins, and G. F. Counsell, J. Nucl. Mater. 337–339, 114 \n
(2005).\n
Martin, J. D. (2006). Theory and Simulation of Dust in Tokamak Plasmas, \n
(July).\n
M. Bacharis, M. Coppins, and J. E. Allen, Phys. Rev. E - Stat. Nonlinear, \n
Soft Matter Phys. 82, 1 (2010).\n
M. Bacharis, M. Coppins, W. Fundamenski, and J. E. Allen, Plasma Phys. \n
Control. Fusion 54, 085010 (2012).\n


\section directories_sec Directory Structure in DTOKSU

Files Included With This Project,\n
build:\n
	FindNetCDF.cmake\n
\n
Config_Files:\n
	DTOKSU_Config.cfg\n
\n
Dependencies:\n
	install_netcdf4.sh\n
\n
src:\n
	Beryllium.cpp   Deuterium.cpp  ForceModel.cpp      HeatingModel.cpp  \n
	MathHeader.cpp  Molybdenum.cpp Tungsten.cpp        ChargingModel.cpp  \n
	DTOKSU.cpp      Functions.cpp  Iron.cpp             Matter.cpp      \n
	solveMOMLEM.cpp	Constants.cpp  DTOKSU_Manager.cpp  Graphite.cpp \n
	Lithium.cpp     Model.cpp      threevector.cpp\n
\n
include:\n
	Beryllium.h    Constants.h     DTOKSU.h    ForceModel.h GrainStructs.h  \n
	HeatingModel.h Lithium.h       Matter.h    Molybdenum.h solveMOMLEM.h \n
	Tungsten.h     ChargingModel.h Deuterium.h DTOKSU_Manager.h \n
	Functions.h    Graphite.h      Iron.h      MathHeader.h Model.h  \n
	PlasmaData.h  threevector.h\n
\n
PlasmaData/PlasmaGenerator:\n
	build.sh PlasmaGenerator.cpp\n
\n
Tests:\n
	IntegrationTests/ ModelTests/ UnitTests/ CMakeLists.txt README_TESTING\n
\n
\n
\secton How to Install/Setup DTOKSU
\n
DTOKS-U depends on several other packages for principally for providing the \n
configuration file handling (config4cpp) and for extracting data from NETCDF \n
files. For this reason, several other dependancies must first be installed. \n
The setup.sh and install_netcdf4.sh scripts will attempt to perform this \n
automatically though are liable to fail! \n
\n
By Default, DTOKS-U is built with netcdf capability enabled. To disable this,\n
you must avoid running the install_netcdf4.sh script, by commenting this line\n
out of the setup.sh script.\n
\n
To Install DTOKSU, follow steps 1 to 5 given below:\n
1) First, we need to download the files. To Git clone the repository\n
\n
	git clone https://github.com/LukeSimons/DTOKS-U.git\n
\n
\n
2) Run the setup script, this will in turn install the dependancies using the \n
install_netcdf4.sh script. IF YOU DON'T WANT THESE DEPENDENCIES, COMMENT THIS\n
LINE OUT OF THE SETUP SCRIPT BEFORE RUNNING.\n
\n
	. setup.sh\n
\n
\n
3) If there are no significant errors or conflicts after this, setup the \n
CMAKE build directory. By default, the NETCDF support is enabled and the test\n
directory is not build. So there are a few options available\n
\n
a) Configure DTOKS core files with config4cpp dependency\n
\n
	cd build\n
	cmake ../. -DBUILD_NETCDF=OFF\n
\n
b) Configure DTOKS core files with config4cpp and build test directory\n
\n
	cd build\n
	cmake ../. -DBUILD_TESTS=ON\n
\n
c) Configure DTOKS core files with config4cpp and NETCDF dependency\n
\n
	cd build\n
	cmake ../.\n
\n
\n
4) DTOKSU can now be built with CMAKE.\n
	make\n
\n
5a) The code can then be run in default mode:\n
\n
	cd ..\n
	./bin/dtoksu -h\n
	./bin/dtoksu\n
\n
\n
5b) Alternatively, if the test directory has been built following step 3b), \n
the tests can be ran\n
\n
	cd ..\n
	./bin/init_test -h\n
	./bin/init_test -m Gravity\n
	./bin/unit_test -h\n
	./bin/unit_test -m DeltaSec\n
\n
\n
\section operation_sec How to Run DTOKSU
The built executable and library files are kept in the /bin/ directory. There \n
are four executables named:\n
	dtoksu  int_test model_test  unit_test\n
\n
For details on how to use the three test executables, see the file:\n
Tests/README_TESTING\n
\n
The 'dtoksu' executable is used to run the code, with the following command\n
line options available:\n
	\n
./bin/dtoks -h <help message> -c <Config File> -t <Inital Temperature> \n
	-m <'w', 'g', 'b', 'd' or 'f'> 	-s <Initial Radius> \n
	-vr <Radial Velocity> -vt <Angular Velocity> -vz <longitudinal velocity>\n
	-rr <Radial position> -rt <Angular position> -rz <longitudinal position>\n
	-op <Output File Pre-fix> -om <MetaData filename>\n
\n
\n
\section classes_sec DTOKSU Class Structure and Design
DTOKSU follows an object oriented programing (oop) style with a few different \n
class heirarchies and data structures being used to manage the information. \n
These have been designed with simplicity and security in mind and with \n
maximal modular capability.\n
\n
The two principal abstract base classes are the Matter.h class and the \n
Model.h class. These classes and their children are operated by the wrapper\n
class DTOKSU to conduct the entire simulation. The class DTOKSU_Manager is \n
used to setup and correctly configure instances of DTOKSU.\n
\n
\subsection matter_subsec Matter class
The Matter class is used to contain all the information about the dust grain \n
but also to simulate all the physical models directly related to the \n
information contained. It owns an instance of both the 'GrainData' and \n
'ElementConsts' which are found in the GrainStructs.h file. This stores the\n
non-constant and constant information about the material. In addition to \n
this, Matter.h stores information on what models are being allowed to vary \n
with temperature. The matter class requires that classes which inherit from \n
it define the functions "update_radius()", "update_heatcapacity()" and \n
"update_vapourpressure()" which define the variation of radius, heat capacity\n
and vapour pressure with temperature.\n
\n
It contains functions able to alter the description of a dust grain including \n
it's motion, these are namely "update_temperature", "update_motion", \n
"update_charge" and "update_mass". Finally, a range of getter and one setter \n
methods exist for classes external to the heirarchy, mainly DTOKSU to access \n
this information. This class requires the definition of the "threevector.h" \n
class, the "Constants.h", additionaly "Functions.h" and obviously the data \n
structures in "GrainStructs.h"\n
\n
The class can be constructed simply by passing a pointer to the element \n
constants data structure which contains all the unmutable information about a \n
materials properties. The default values for it's variable properaties are \n
set by default.\n
\n
Inforamtion about the plasma and dust particle are held within the \n
PlasmaData.h and GrainStructs.h structures. Constants.h contains physical \n
constants and preprocessor directives for debugging Functions.h contains a \n
warning function, a rounding function and as well as methods used to \n
calculate physical models. The DTOKSU_Manager.h class is used to read in and\n
hold the information about a plasma background with spatial dependance over a\n
cartesian grid. The PlasmaGenerator.cpp file is used to generate grids of \n
data for different plasma backgrounds which can be fed to the PlasmaGrid \n
class.\n
\n
\subsection model_subsec Model.h class
The Model.h class is used as the base for solving the four principal \n
equations governing the behaviour of the dust. These are the force balance, \n
current balance and the equation for the heat flux. Instances of the Model \n
class hold a pointer to the Matter class which is established to be the \n
pointer to the same Matter object as initialised by DTOKSU_Manager. They \n
also share a pointer to the same PlasmaData object. The Model class has \n
information about the Plasma Grid and has member data specifying the \n
accuracy of the models, the time step, and total time elapsed as measured by\n
that model. \n
\n
The Model class owns two file streams, one of which is used to record the \n
plasma data at each time step and one of which records the model data.\n
\n
\subsection plasmadata_subsec PlasmaData/PlasmaBackground
First, the dust particle is located within the Plasma grid and assigned \n
integer coordinates giving it's position. The plasma velocity is defined by \n
the direction of the magnetic field and the mean plasma velocity provided the \n
densities sum is non-zero. The value of the potential on the grid points are \n
taken from input files for each machine type (MAST, ITER, Magnum-PSI, EAST or\n
DIII-D). The electric field in the radial direction is calculated from the \n
difference in potential between two points and similarly calculated for the \n
other three dimensions in cylindrical coordinates.\n
\n
-----------------------------------------------------------------------------\n
\n
-----------------------------------------------------------------------------\n
*/
