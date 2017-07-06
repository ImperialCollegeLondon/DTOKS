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

	clock_t begin = clock();
	std::array<char,4> ConstModels  = {'c','c','c','y'};
//	Matter *Sample1 = new Tungsten(1e-6,300,ConstModels);

	// ***************************************** INITIALISE PLASMA ***************************************** //
	char Plasma='h';
	char Machine='m';	// Initialise plasma grid
	double Spacing =0.01;	// metres

	PlasmaGrid Pgrid(Plasma,Machine,Spacing);

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
	char Element='B';		// Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
//	double Power=0;			// Kilo-Watts power in addition to heating model powers
	double Size=1e-6; 		// m
	double Temp=300;		// K
	double TimeStep=1e-5;		// s
//	std::shared_ptr<Matter> Sample;	// Define the sample matter type

	Matter * Sample;		// Define the sample matter type


	// ------------------- HEATING MODELS ------------------- //
	// Set to true all heating models that are wanted
	bool RadiativeCooling = 1;
	bool EvaporativeCooling = 1;
	bool NewtonCooling = 0;			// This model is equivalent to Electron and Ion heat flux terms
	bool NeutralHeatFlux = 1; 		// Plasma heating terms
	bool ElectronHeatFlux = 1;
	bool IonHeatFlux = 1;
	bool NeutralRecomb = 1;
	bool TEE = 1;				// Electron Emission terms
	bool SEE = 1;
	bool PlasmaHeating = 1; 		// If we want plasma heating terms turned off
	// NOTE: For Negative dust with RE=0 and Te = Ti, the Ion and Electron heat flux will be identical!
	if( !PlasmaHeating ){
		NeutralHeatFlux = false;
		ElectronHeatFlux = false;
		IonHeatFlux = false;
		NeutralRecomb = false;
	}

	// ------------------- FORCING MODELS ------------------- //
        bool Gravity = 1;
        bool Centrifugal = 1;
        bool Lorentz = 1;
        bool IonDrag = 1;	
        bool NeutralDrag = 1;	

	// ------------------- CHARGING MODELS ------------------- //
        bool DTOKSOML = true;

	// Plasma Data
	// NOTE: 16/06/17, current plasma data matches DTOKS initial plasma data for debugging purposes.

	PlasmaData *Pdata = new PlasmaData;
	Pdata->NeutralDensity =  1e18;//10e17;	// m^-3, Neutral density
	Pdata->IonDensity =  3.8e18;//10e17;	// m^-3, Ion density
	Pdata->ElectronDensity = 1e18;//10e18; 	// m^-3, Electron density
	double NumOfev = 1.5;
	Pdata->IonTemp = NumOfev*1.16e4;	// K, Ion Temperature
	Pdata->ElectronTemp = 1.188*1.16e4; 	// K, Electron Temperature, convert from eV
	Pdata->NeutralTemp = NumOfev*1.16e4; 	// K, Neutral Temperature, convert from eV
	threevector PlasmaVelocity(501.33, 7268.5, 914.947); // Taken from initial for DTOKS
	Pdata->PlasmaVel = PlasmaVelocity;
	threevector Efield(-13.673, 0, -27.925);
	Pdata->ElectricField = Efield;
	threevector Bfield(0.0226868, 0.328923, 0.0414043);
	Pdata->MagneticField = Bfield;

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
	std::array<bool,5> ForceModels  = {Gravity,Centrifugal,Lorentz,IonDrag,NeutralDrag};
	std::array<bool,1> ChargeModels = {DTOKSOML};
	ConstModels  = {EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};
	
	// Accuracy Levels correspond to Charging, Heating and Forcing respectively
	std::array<double,3> AccuracyLevels = {1.0,0.1,0.1};

	if 	(Element == 'W') Sample = new Tungsten(Size,Temp,ConstModels);
	else if (Element == 'B') Sample = new Beryllium(Size,Temp,ConstModels);
	else if (Element == 'F') Sample = new Iron(Size,Temp,ConstModels);
	else if (Element == 'G') Sample = new Graphite(Size,Temp,ConstModels);
	else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}

	threevector xinit(1.15,0.0,-1.99);// default injection right hand side
	threevector vinit(0.0,0.0,100.0);
	Sample->update_motion(xinit,vinit);

	std::cout << "\n\n * GENERATE DTOKS * \n\n";
//	DTOKSU MyDtoks1(TimeStep, AccuracyLevels, Sample, Pdata, HeatModels, ForceModels, ChargeModels);
	DTOKSU MyDtoks2(TimeStep, AccuracyLevels, Sample, Pgrid, HeatModels, ForceModels, ChargeModels);
	
	std::cout << "\n\n * RUN DTOKS * \n\n";
//	int errcode1 = MyDtoks1.Run();
	int errcode2 = MyDtoks2.Run();

	std::cout << "\n\n * MAIN SCRIPT COMPLETE * \n\n";
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nCOMPLETED IN : " << elapsd_secs << "s\n\n";

	return 0;
}
