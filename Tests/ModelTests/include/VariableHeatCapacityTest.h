#include "HeatingModel.h"

int VariableHeatCapacityTest(char Element, bool VaryHeatCapacity){
	clock_t begin = clock();
	// ********************************************************** //
	// FIRST, define program default behaviour

	// Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
	char EmissivityModel = 'c'; 	// Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
	char ExpansionModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
													// and (z)ero expansion

	char HeatCapacityModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
	if (VaryHeatCapacity) HeatCapacityModel = 'v';

	char BoilingModel = 'y'; 	// Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
													// and (t)homson
	char TimeStepType = 'f';	// Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
													// and (f)ixed steps
	std::string Name="constant";	// Describes heating model

 	// Parameters describing the heating model
	double Size=1e-6; 		// m
	double Temp=280;		// K
	double TimeStep=1e-9;		// s
	double Potential = 1;
	Matter *Sample;			// Define the sample matter type

	// Set to true all heating models that are wanted
	bool RadiativeCooling = true;
	bool EvaporativeCooling = false;
	bool NewtonCooling = false;		// This model is equivalent to Electron and Ion heat flux terms
	// Plasma heating terms
	bool NeutralHeatFlux = false;
	bool ElectronHeatFlux = true;
	bool IonHeatFlux = true;
	bool NeutralRecomb = true;
	// Electron Emission terms
	bool TEE = true;
	bool SEE = true;

	PlasmaData *Pdata = new PlasmaData;
	Pdata->NeutralDensity	= 1e18; 	// m^-3, Neutral Density
	Pdata->IonDensity	= 1e18; 	// m^-3, Ion Density
	Pdata->ElectronDensity	= 1e18;	// m^-3, Electron Density
	Pdata->NeutralTemp	= 10*1.16e4;	// K, Neutral Temperature, convert from eV
	Pdata->IonTemp		= 10*1.16e4;	// K, Ion Temperature, convert from eV
	Pdata->ElectronTemp	= 10*1.16e4;	// K, Electron Temperature, convert from eV
	Pdata->AmbientTemp = 0;

	// Models and ConstModels are placed in an array in this order:
	std::array<bool, HMN> Models = 
		{RadiativeCooling, EvaporativeCooling, NewtonCooling, IonHeatFlux, ElectronHeatFlux, NeutralHeatFlux, 
		NeutralRecomb, SEE, TEE, false };
	std::array<char, CM> ConstModels =
		{ EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel,'n'};

	if	(Element == 'W'){ 
		Sample = new Tungsten(Size,Temp,ConstModels);
	}else if (Element == 'B'){ 
		Sample = new Beryllium(Size,Temp,ConstModels);
	}else if (Element == 'F'){
		Sample = new Iron(Size,Temp,ConstModels);
	}else if (Element == 'G'){
		Sample = new Graphite(Size,Temp,ConstModels);
	}else if (Element == 'D'){
                Sample = new Deuterium(Size,Temp,ConstModels);
        }else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}
	threevector xinit(1.15,0.0,-1.99);// default injection right hand side
	threevector vinit(0.0,0.0,0.0);
	Sample->update_motion(xinit,vinit,0.0);

	Sample->set_potential(Potential);

	std::string filename;
	if( VaryHeatCapacity ) 	filename = "Data/HeatCapacityNonConst.txt";
	else 			filename = "Data/HeatCapacityConst.txt";

	HeatingModel MyModel(filename,1.0,Models,Sample,Pdata);
	MyModel.UpdateTimeStep();
	MyModel.Vapourise();
	
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

	std::cout << "\n\n*****\n\nErrorEstimateTest 5 completed in " << elapsd_secs << "s\n";

	return 0.0;
}
