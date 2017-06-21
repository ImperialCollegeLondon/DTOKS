#include "HeatingModel.h"

int ConstantHeatingTest(char Element){
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
	char TimeStepType = 'f';	// Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
													// and (f)ixed steps
	std::string Name="constant";	// Describes heating model

 	// Parameters describing the heating model
	double Power=1e-8;		// Kilo-Watts power in addition to heating model powers
	double Size=5e-8; 		// m
	double Temp=280;		// K
	double TimeStep=1e-10;		// s
	Matter *Sample;			// Define the sample matter type

	// Set to true all heating models that are wanted
	bool RadiativeCooling = false;
	bool EvaporativeCooling = false;
	bool NewtonCooling = false;		// This model is equivalent to Electron and Ion heat flux terms
	// Plasma heating terms
	bool NeutralHeatFlux = true;
	bool ElectronHeatFlux = true;
	bool IonHeatFlux = true;
	bool NeutralRecomb = true;
	// Electron Emission terms
	bool TEE = false;
	bool SEE = false;

	PlasmaData Pdata;

	bool PlasmaHeating = false; 		// If we want plasma heating terms turned off
	if( !PlasmaHeating ){
		NeutralHeatFlux = false;
		ElectronHeatFlux = false;
		IonHeatFlux = false;
		NeutralRecomb = false;
	}

	// Models and ConstModels are placed in an array in this order:
	std::array<bool, 9> Models = 
		{RadiativeCooling, EvaporativeCooling, NewtonCooling, IonHeatFlux, ElectronHeatFlux, NeutralHeatFlux, 
		NeutralRecomb, SEE, TEE };
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

	HeatingModel MyModel("out_ConstantHeatingTest.txt",1.0,Models,Sample,Pdata);

	double Mass = Sample->get_mass();
	MyModel.set_PowerIncident(Power);
	MyModel.UpdateTimeStep();

	MyModel.Vapourise();

	double ModelTime = MyModel.get_totaltime();
	
	double p1 = (Mass/Power)*((Sample->get_superboilingtemp()-Temp)*Sample->get_heatcapacity());
	double p2 = (Mass*Sample->get_latentfusion())/Power;
	double p3 = (Mass*Sample->get_latentvapour())/Power;
//	std::cout << "\np1 = " << p1 << "\np2 = " << p2 << "\np3 = " << p3;
	double AnalyticTime = p1 + p2 + p3;
	
	double ReturnVal = 0;

	if( ModelTime == AnalyticTime ) 			ReturnVal = 1;
	else if( fabs(1-ModelTime/AnalyticTime) < 0.01 ) 	ReturnVal = 2;
	else							ReturnVal = -1;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModelTime = " << ModelTime << "s : AnalyticTime = " << AnalyticTime << "s";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelTime/AnalyticTime) <<"%\n*****\n\n";

	delete Sample;

	return ReturnVal;
}
