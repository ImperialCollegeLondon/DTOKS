#include "ForceModel.h"

int ConstantForceTest(char Element){
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
	double Size=5e-8; 		// m
	double Temp=280;		// K
	double Potential = 1;		// Normalised Potential
	Matter *Sample;			// Define the sample matter type

	// Set to true all Force models that are wanted
	bool Gravity = true;		// IS ON!
	bool Centrifugal = false;	// IS OFF!
	bool Lorentz = true;		// IS ON!
	bool IonDrag = false;		// IS OFF!

	PlasmaData Pdata;
//	threevector PlasmaVelocity(10, 5, 3); // Taken from initial for DTOKS
//	Pdata.PlasmaVel = PlasmaVelocity;
	Pdata.ElectronTemp = 100*1.16e4;	// K, Electron Temperature, convert from eV
	threevector GravityForce(0, 0, -9.81);
	Pdata.Gravity = GravityForce;
	threevector Efield(1, -2, 3);
	Pdata.ElectricField = Efield;
	threevector Bfield(0.0, 0.0, 0.0);
	Pdata.MagneticField = Bfield;

	std::array<bool,4> ForceModels  = {Gravity,Centrifugal,Lorentz,IonDrag};

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
	ForceModel MyModel("out_ConstantForcingTest.txt",1.0,ForceModels,Sample,Pdata);

	double Mass = Sample->get_mass();
	MyModel.UpdateTimeStep();
	size_t imax(100);

	for( size_t i(0); i < imax; i ++)
		MyModel.Force();

	threevector ModelVelocity = Sample->get_velocity();

	double ConvertKelvsToeV(8.621738e-5);
	// NOTE: this is the charge to mass ratio for a negative grain only...
	double qtom = -3.0*epsilon0*Pdata.ElectronTemp*ConvertKelvsToeV*Sample->get_potential()
			/(pow(Sample->get_radius(),2)*Sample->get_density());
	threevector AnalyticVelocity = (Efield*qtom + Pdata.Gravity)*(imax*MyModel.get_timestep());
	
	double ReturnVal = 0;

	if( (ModelVelocity - AnalyticVelocity).mag3() == 0 ) 				ReturnVal = 1;
	else if( (ModelVelocity - AnalyticVelocity).mag3()/ModelVelocity.mag3() < 0.01 ) 	ReturnVal = 2;
	else										ReturnVal = -1;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModelVel = " << ModelVelocity << "m s^-1 : AnalyticVel = " << AnalyticVelocity << "m s^-1";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelVelocity.mag3()/AnalyticVelocity.mag3()) <<"%\n*****\n\n";

	delete Sample;

	return ReturnVal;
}
