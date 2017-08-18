#include "ForceModel.h"
#include "MathHeader.h" // This is weird

int NeutralDragTest(char Element, bool DragSwitch){
	clock_t begin = clock();
	// ********************************************************** //
	// FIRST, define program default behaviour

	// Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
	char EmissivityModel = 'c'; 	// Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
	char ExpansionModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
													// and (z)ero expansion
	char HeatCapacityModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
	char BoilingModel = 'y'; 	// Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
													// and (t)homson

 	// Parameters describing the heating model
	double Size=1e-6; 		// m
	double Temp=280;		// K
	double Potential = 1;		// Normalised Potential
	Matter *Sample;			// Define the sample matter type

	// Set to true all Force models that are wanted
	bool Gravity = true;		// IS ON!
	bool Centrifugal = true;	// IS ON!
	bool Lorentz = true;		// IS ON!
	bool IonDrag = true;		// IS ON!
	bool NeutralDrag = DragSwitch;	// IS ON!

	PlasmaData *Pdata = new PlasmaData();
	Pdata->NeutralDensity	= 10e18; 	// m^-3, Neutral Density
	Pdata->IonDensity	= 10e18; 	// m^-3, Ion Density
	Pdata->ElectronDensity	= 10e18;	// m^-3, Electron Density
	Pdata->NeutralTemp	= 10*1.16e4;	// K, Neutral Temperature, convert from eV
	Pdata->IonTemp		= 10*1.16e4;	// K, Ion Temperature, convert from eV
	Pdata->ElectronTemp	= 10*1.16e4;	// K, Electron Temperature, convert from eV
	threevector PlasmaVel(0.0,0.0,0.5*sqrt(Kb*Pdata->NeutralTemp/Mp));
	Pdata->PlasmaVel	= PlasmaVel;

	std::array<bool,5> ForceModels  = {Gravity,Centrifugal,Lorentz,IonDrag,NeutralDrag};

	// Models and ConstModels are placed in an array in this order:
	std::array<char, 4> ConstModels =
		{ EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};

	if	(Element == 'W'){ 
		Sample = new Tungsten(Size,Temp,ConstModels);
	}else if (Element == 'B'){ 
		Sample = new Beryllium(Size,Temp,ConstModels);
	}else if (Element == 'F'){
		Sample = new Iron(Size,Temp,ConstModels);
	}else if (Element == 'G'){
		Sample = new Graphite(Size,Temp,ConstModels);
	}else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}
	Sample->set_potential(Potential);
	threevector vinit(2.0,3.0,5.0);
	threevector zeros(1.0,0.0,0.0);
	Sample->update_motion(zeros,vinit,0.0);
	// START NUMERICAL MODEL
	std::string filename;
	if( DragSwitch ) 	filename = "Data/NeutralDragOn.txt";
	else 			filename = "Data/NeutralDragOff.txt";
	ForceModel MyModel(filename,0.1,ForceModels,Sample,Pdata);
	double Mass = Sample->get_mass();
	MyModel.UpdateTimeStep();
	size_t imax(100000);

	for( size_t i(0); i < imax; i ++)
		MyModel.Force();

	threevector ModelVelocity = Sample->get_velocity();

	delete Pdata;
	delete Sample;
	// END NUMERICAL MODEL
	
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nErrorEstimateTest 1 completed in " << elapsd_secs << "s\n";

	return 1;
}
