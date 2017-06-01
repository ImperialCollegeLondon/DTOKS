#include "DTOKSU.h"
#include <vector>

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "Options:\n"
	<< "\t-h,--help\t\tShow this help message\n\n"
	<< "\t-d,--deltatime DELTATIME\t\tSpecify time step\n\n"
	<< "\t-t,--temperature TEMPERATURE\tdouble variable defining the initial temperature\n"
	<< "\t-m,--material MATERIAL\t\tchar variable giving the dust grain element, possible values 'w', 'g', 'b' and 'f'\n"
	<< "\t\t\t\t\t(W): Tungsten, (G): Graphite, (B): Beryllium or (F): Iron\n\n"
	<< "\t-s,--size SIZE\t\t\tdouble variable giving the radius of the grain\n\n"
	<< "\t-b,--boiling BOILING\t\tchar variable with possible values 'y', 'n', 's' and 't'Corresponding to :\n"
	<< "\t\t\t\t\t(y)es, (n)o, (s)uper and (t)homson\n"
	<< "\t-e,--evaporation EVAPORATION\tboolean variable, true turns evaporation on, false turns evaporation off.\n\n"
	<< "\t-p,--plasmatemp PLASMATEMP\tdouble variable, determines temperature of plasma in ev.\n\n"
	<< "\t-w,--power POWER\tdouble variable, determines additional heating in kW dust is subject to.\n\n"
	<< std::endl;
}

template<typename T> int InputFunction(int &argc, char* argv[], int &i, std::stringstream &ss0, T &Temp){
	if (i + 1 < argc) { // Make sure we aren't at the end of argv!
		i+=1;
		ss0 << argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
		ss0 >> Temp;
		ss0.clear(); ss0.str("");
		return 0;
	}else{ // Uh-oh, there was no argument to the destination option.
		std::cerr << "\noption requires argument." << std::endl;
		return 1;
	}

}
/*
	std::cout << "\n\n************************************* BEGIN SETUP (1) ************************************* \n\n";

	DTOKSU MyDtoks;
	int errcode = MyDtoks.Run();

	std::cout << "\n\n************************************* BEGIN SETUP (2) ************************************* \n\n";
*/
int main(int argc, char* argv[]){


	// ***************************************** INITIALISE PLASMA ***************************************** //
	char machine='m';	// Initialise plasma grid
	plasmagrid *Background = new plasmagrid('h',machine,0.01);	

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
	char Element='W';		// Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
//	double Power=0;			// Kilo-Watts power in addition to heating model powers
	double Size=1e-6; 		// m
	double Temp=350;		// K
	double TimeStep=1e-5;		// s
//	std::shared_ptr<Matter> Sample;	// Define the sample matter type
	Matter * Sample;		// Define the sample matter type

	// ------------------- HEATING MODELS ------------------- //
	// Set to true all heating models that are wanted
	bool RadiativeCooling = false;
	bool EvaporativeCooling = false;
	bool NewtonCooling = false;		// This model is equivalent to Electron and Ion heat flux terms
	bool NeutralHeatFlux = true; 		// Plasma heating terms
	bool ElectronHeatFlux = true;
	bool IonHeatFlux = true;
	bool NeutralRecomb = false;
	bool TEE = false;			// Electron Emission terms
	bool SEE = false;
	bool PlasmaHeating = true; 		// If we want plasma heating terms turned off
	if( !PlasmaHeating ){
		NeutralHeatFlux = false;
		ElectronHeatFlux = false;
		IonHeatFlux = false;
		NeutralRecomb = false;
	}

	// ------------------- FORCING MODELS ------------------- //
        bool Gravity = false;
        bool Centrifugal = false;
        bool Lorentz = false;
        bool IonDrag = false;	

	// ------------------- CHARGING MODELS ------------------- //
        bool DTOKSOML = false;

	// Plasma Data
	PlasmaData Pdata;
	Pdata.NeutralDensity =  10e19;//10e17;	// m^-3, Neutral density
	Pdata.ElectronDensity = 10e19;//10e18; 	// m^-3, Electron density
	double NumOfev = 10;
	Pdata.IonTemp = NumOfev*1.16e5;	 	// K, Ion Temperature
	Pdata.ElectronTemp = NumOfev*1.16e5; 	// K, Electron Temperature, convert from eV
	Pdata.NeutralTemp = NumOfev*1.16e5; 	// K, Neutral Temperature, convert from eV


	std::vector <std::string> sources;
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 	|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--deltatime" 	|| arg == "-d" )	InputFunction(argc,argv,i,ss0,TimeStep);
		else if( arg == "--temperature" || arg == "-t" )	InputFunction(argc,argv,i,ss0,Temp);
		else if( arg == "--material" 	|| arg == "-m" ) 	InputFunction(argc,argv,i,ss0,Element);
		else if( arg == "--size" 	|| arg == "-s" )	InputFunction(argc,argv,i,ss0,Size);
		else if( arg == "--boiling" 	|| arg == "-b" )	InputFunction(argc,argv,i,ss0,BoilingModel);
		else if( arg == "--evaporation" || arg == "-e" )	InputFunction(argc,argv,i,ss0,EvaporativeCooling);
		else if( arg == "--plasmatemp" 	|| arg == "-p" )	InputFunction(argc,argv,i,ss0,NumOfev);
//		else if( arg == "--power" 	|| arg == "-w" )	InputFunction(argc,argv,i,ss0,Power);
                else{
			sources.push_back(argv[i]);
		}
	}
	// Models and ConstModels are placed in an array in this order:
	std::array<bool, 9> HeatModels = 
		{RadiativeCooling, EvaporativeCooling, NewtonCooling, IonHeatFlux, ElectronHeatFlux, NeutralHeatFlux, 
		NeutralRecomb, SEE, TEE };
	std::array<bool,4> ForceModels  = {Gravity,Centrifugal,Lorentz,IonDrag};
	std::array<bool,1> ChargeModels = {DTOKSOML};
	std::array<char,4> ConstModels  = {EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};

	std::array<double,3> AccuracyLevels = {1.0,1.0,1.0};

	if 	(Element == 'W') Sample = new Tungsten(Size,Temp,ConstModels);
	else if (Element == 'B') Sample = new Beryllium(Size,Temp,ConstModels);
	else if (Element == 'F') Sample = new Iron(Size,Temp,ConstModels);
	else if (Element == 'G') Sample = new Graphite(Size,Temp,ConstModels);
	else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}

//	DTOKSU MyDtoks2(TimeStep, AccuracyLevels, Sample, Pdata, HeatModels, ForceModels, ChargeModels);
	DTOKSU MyDtoks2(TimeStep, AccuracyLevels, Sample, Background, HeatModels, ForceModels, ChargeModels);
	
	int errcode2 = MyDtoks2.Run();

	std::cout << "\n\n * MAIN SCRIPT COMPLETE * \n\n";
	return 0;
}
