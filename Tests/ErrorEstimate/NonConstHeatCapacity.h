#include "Constants.h"
#include <math.h>
#include <iostream>
#include "HeatingModel.h"

double HeatCapacityInt(double Temp){
	double Ret = 0.54212*Temp-2.42667e-6*pow(Temp,2)/2-90.2725*log10(Temp)-43449.3/Temp
					-1.59309e7/(2*pow(Temp,2))+1.43688e9/(3*pow(Temp,3));

	return Ret*4.184;
}

int NonConstHeatCapacity(char Element){
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
	double TimeStep=1e-12;		// s
	Matter *Sample;			// Define the sample matter type

	// Set to true all heating models that are wanted
	bool RadiativeCooling = false;
	bool EvaporativeCooling = false;
	bool NewtonCooling = false;		// This model is equivalent to Electron and Ion heat flux terms
	// Plasma heating terms
	bool NeutralHeatFlux = false;
	bool ElectronHeatFlux = false;
	bool IonHeatFlux = false;
	bool NeutralRecomb = false;
	// Electron Emission terms
	bool TEE = false;
	bool SEE = false;

	PlasmaData Pdata;
	Pdata.NeutralDensity = 3e19;		// m^-3, Neutral density
	Pdata.ElectronDensity = 8e17;	 	// m^-3, Electron density
	Pdata.Potential = 1;			// arb, assumed negative, potential normalised to dust temperature, (-e*phi)/(Kb*Td)
	Pdata.IonTemp = 100*1.16e5;	 	// K, Ion Temperature
	Pdata.ElectronTemp = 100*1.16e5;	// K, Electron Temperature, convert from eV
	Pdata.NeutralTemp = 100*1.16e5; 	// K, Neutral Temperature, convert from eV
	Pdata.AmbientTemp = 0;

	// Models and ConstModels are placed in an array in this order:
	std::array<bool, 9> Models = 
		{RadiativeCooling, EvaporativeCooling, NewtonCooling, IonHeatFlux, ElectronHeatFlux, NeutralHeatFlux, 
		NeutralRecomb, SEE, TEE };
	std::array<char, 4> ConstModels =
		{ EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};

	if	(Element == 'W'){ 
		Sample = new Tungsten(Size,Temp,ConstModels);
		TimeStep=1e-12;
	}else if (Element == 'B'){ 
		Sample = new Beryllium(Size,Temp,ConstModels);
		TimeStep=1e-11;
	}else if (Element == 'F'){
		Sample = new Iron(Size,Temp,ConstModels);
		TimeStep=1e-13;
	}else if (Element == 'G'){
		Sample = new Graphite(Size,Temp,ConstModels);
		TimeStep=1e-11;
	}else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}

	HeatingModel MyModel(Name,Element,Power,Models,TimeStep,Sample,Pdata);
//	MyModel.Vapourise(("out_ConstantHeatingTest_" + Element + ".txt").c_str(),ConstModels,TimeStepType);
	MyModel.Vapourise("out_ConstantHeatingTest.txt",ConstModels,TimeStepType);

	double ModelTime = MyModel.get_totaltime();
	double ModelEnergy = ModelTime*Power; // In kJ
	std::cout << "\nModel Time = " << ModelTime;
	std::cout << "\nModel Energy = " << ModelEnergy;

	// *********************** BEGIN ANALYTICAL MODEL ************************************ //

	int ReturnVal = 0;

	double mass = 1e-9;

	double Ti = 280;
	double Tf = 4000;
	double Cv0 = HeatCapacityInt(280);
	double Cv1 = HeatCapacityInt(4000);
	
	double AnalyticalEnergy1 = mass*(Cv0+Cv1);
//	std::cout << "\nTotal Energy = " << Energy << "kJ"; 

	double AnalyticalEnergy2 = mass*2.0*(Tf-Ti);
//	std::cout << "\nTotal Energy = " << Energy << "kJ"; 

	// ********************************************************** //
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nErrorEstimate 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModel energy = " << ModelEnergy << "kJ : Analytic Energy1 = " << AnalyticalEnergy1 << "kJ : Analytic energy 2 = " << AnalyticalEnergy2;
	std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelEnergy/AnalyticalEnergy1) <<"%\n*****\n\n";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelEnergy/AnalyticalEnergy2) <<"%\n*****\n\n";

	return ReturnVal;
}
