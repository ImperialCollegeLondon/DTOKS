#include "DTOKSU.h"
#include "Breakup.h"

#include <vector>
#include <config4cpp/Configuration.h>

static void show_usage(std::string name){
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h,--help\t\tShow this help message\n\n"
	<< "\t-t,--temperature TEMPERATURE\tdouble variable defining the initial temperature\n"
	<< "\t-m,--material MATERIAL\t\tchar variable giving the dust grain element, possible values 'w', 'g', 'b' and 'f'\n"
	<< "\t\t\t\t\t(W): Tungsten, (G): Graphite, (B): Beryllium or (F): Iron\n\n"
	<< "\t-s,--size SIZE\t\t\tdouble variable giving the radius of the grain\n\n";
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

int main(int argc, char* argv[]){

	clock_t begin = clock();

	// ------------------- PARSE CONFIGURATION FILE ------------------- //
	std::string Filename = "Data/DTOKSU.txt";
	char Plasma='h';
	char Machine='p';	// Initialise plasma grid
	float xSpacing = 0.00234375;	// 0.15/64.0; x-dimensional spacing (Radial) metres
	float zSpacing = 0.04; 	// 1.0/25.0;  y-dimensional spacing (z) metres
	
	char Element='W';
	float size=0.5e-6;
	float Temp=300;

	char EmissivityModel = 'c';
	char ExpansionModel = 'c';
	char HeatCapacityModel = 'v'; 
	char BoilingModel = 't';

	// ------------------- GROUP MODELS ------------------- //
	std::array<bool, 9> HeatModels;
	std::array<bool,6> ForceModels;
	std::array<bool,3> ChargeModels;
	std::array<char,4> ConstModels;
	std::array<float,3> AccuracyLevels;
	PlasmaData *Pdata = new PlasmaData;

// ------------------- PROCESS USER-INPUT ------------------- //
	std::vector <std::string> sources;
	std::string Config_Filename = "DTOKSU_Config.cfg";
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if( arg == "--config"      || arg == "-c" )	InputFunction(argc,argv,i,ss0,Config_Filename);
		else{
			sources.push_back(argv[i]);
		}
	}
	threevector PlasmaVelocity(501.33, 7268.5, 914.947); // Taken from initial for DTOKS
	threevector Efield(-13.673, 0, -27.925);
	threevector Bfield(0.0226868, 0.328923, 0.0414043);
	Pdata->PlasmaVel 		= PlasmaVelocity;
	Pdata->ElectricField 	= Efield;
	Pdata->MagneticField 	= Bfield;
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	try {
        cfg->parse(Config_Filename.c_str());
        Filename = cfg->lookupString("", "Filename");
        Plasma   = cfg->lookupString("plasmagrid", "Plasma")[0];
        Machine  = cfg->lookupString("plasmagrid", "Machine")[0];
        xSpacing = cfg->lookupFloat("plasmagrid", "xSpacing");
        zSpacing = cfg->lookupFloat("plasmagrid", "zSpacing");
		Element  = cfg->lookupString("dust", "Element")[0];
        size     = cfg->lookupFloat("dust", "size");
        Temp     = cfg->lookupFloat("dust", "Temp");
        ConstModels = 
			{
				cfg->lookupString("variablemodels", "EmissivityModel")[0], cfg->lookupString("variablemodels", "ExpansionModel")[0], 
				cfg->lookupString("variablemodels", "HeatCapacityModel")[0], cfg->lookupString("variablemodels", "BoilingModel")[0]
			};
        HeatModels = 
			{
				cfg->lookupBoolean("heatingmodels","RadiativeCooling"), cfg->lookupBoolean("heatingmodels","EvaporativeCooling"), 
				cfg->lookupBoolean("heatingmodels","NewtonCooling"), cfg->lookupBoolean("heatingmodels","IonHeatFlux"), 
				cfg->lookupBoolean("heatingmodels","ElectronHeatFlux"), cfg->lookupBoolean("heatingmodels","NeutralHeatFlux"), 
				cfg->lookupBoolean("heatingmodels","NeutralRecomb"), cfg->lookupBoolean("heatingmodels","SEE"), 
				cfg->lookupBoolean("heatingmodels","TEE")
			};
		ForceModels =
			{
				cfg->lookupBoolean("forcemodels","Gravity"), cfg->lookupBoolean("forcemodels","Centrifugal"),
				cfg->lookupBoolean("forcemodels","Lorentz"), cfg->lookupBoolean("forcemodels","IonDrag"),
				cfg->lookupBoolean("forcemodels","HybridDrag"), cfg->lookupBoolean("forcemodels","NeutralDrag")
			};
		ChargeModels =
			{
				cfg->lookupBoolean("chargemodels","DTOKSOML"), cfg->lookupBoolean("chargemodels","SchottkyOML"),
				cfg->lookupBoolean("chargemodels","DTOKSWell")
			};
		
		AccuracyLevels = 
			{
				cfg->lookupFloat("accuracylevels", "force"), cfg->lookupFloat("accuracylevels", "heat"), 
				cfg->lookupFloat("accuracylevels", "charge")
			};
		Pdata->IonDensity 		= cfg->lookupFloat("plasmadata","IonDensity");
		Pdata->ElectronDensity 	= cfg->lookupFloat("plasmadata","ElectronDensity");
		Pdata->NeutralDensity 	= cfg->lookupFloat("plasmadata","NeutralDensity");
		Pdata->IonTemp 			= cfg->lookupFloat("plasmadata","IonTemp");
		Pdata->ElectronTemp 	= cfg->lookupFloat("plasmadata","ElectronTemp");
		Pdata->NeutralTemp 		= cfg->lookupFloat("plasmadata","NeutralTemp");

    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        cfg->destroy();
        return 1;
    }

//	PlasmaGrid Pgrid(Plasma,Machine,xSpacing,zSpacing);
    PlasmaGrid Pgrid('h','p',0.00234375,0.04);


// ------------------- PROCESS USER-INPUT ------------------- //
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 		|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--temperature" || arg == "-t" )	InputFunction(argc,argv,i,ss0,Temp);
		else if( arg == "--material" 	|| arg == "-m" ) 	InputFunction(argc,argv,i,ss0,Element);
		else if( arg == "--size" 		|| arg == "-s" )	InputFunction(argc,argv,i,ss0,size);
		else{
			sources.push_back(argv[i]);
		}
	}


    // ------------------- INITIALISE META_DATA FILE ------------------- //
    std::ofstream MetaDataFile;	// Data file for containing the run information
    time_t now = time(0);		// Get the time of simulation
	char * dt = ctime(&now);
    MetaDataFile.open(Filename);
	MetaDataFile << "## Run Data File ##\n";
	MetaDataFile << "#Date:\t" << dt;


	// ------------------- INITIALISE DUST ------------------- //
	Matter * Sample;		// Define the sample matter type
	if 	(Element == 'W')     Sample = new Tungsten(size,Temp,ConstModels);
	else if (Element == 'B') Sample = new Beryllium(size,Temp,ConstModels);
	else if (Element == 'F') Sample = new Iron(size,Temp,ConstModels);
	else if (Element == 'G') Sample = new Graphite(size,Temp,ConstModels);
	else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}

	// Coordinates are r, theta, z.
	// z is the vertical direction.
	// r and theta map out the toroidal plane.
	// Therefore a plot in r, z provides a poloidal cross section
//	CAUTION! WHEN THETA IS SET TO 0.0, WE CAN GET NEGATIVE RADII
	threevector xinit(0.145,0.01,0.05);	// Coordinates are r,theta,z
//	threevector xinit(1.15,0.0,-1.99);// default injection right hand side
	threevector vinit(0.0,0.0,0.0);
	double InitRotationalFreq(0.0);
	Sample->update_motion(xinit,vinit,InitRotationalFreq);


	// ------------------- PRINT METADATA ------------------- //
	MetaDataFile << "\n\n#DUST PARAMETERS" 
		<<"\nElem (arb)\tRadius (m)\tTemp (K)\txinit (m s^-1)\t\tvinit (m s^-1)\n"
		<<Element<<"\t\t"<<size<<"\t\t"<<Temp<<"\t\t"<<xinit<<"\t\t"<<vinit<<"\n"
		<<"\n\n#PLASMA PARAMETERS"
		<<"\nMachine\tgas\tgridx\tgridz\tgridtheta\tdlx\tdlz\n"
		<<Pgrid.get_machine()<<"\t"<<Pgrid.get_gas()<<"\t"<<Pgrid.get_gridx()<<"\t"<<Pgrid.get_gridz()<<"\t"
		<<Pgrid.get_gridtheta()<<"\t"<<Pgrid.get_dlx()<<"\t"<<Pgrid.get_dlz()<<"\n"
		<<"\nxmin (m)\txmax (m)\tzmin (m)\tzmax (m)\n"
		<<Pgrid.get_gridxmin()<<"\t\t"<<Pgrid.get_gridxmax()<<"\t\t"<<Pgrid.get_gridzmin()<<"\t\t"<<Pgrid.get_gridzmax()
		<<"\n\nNn (m^-3)\tNi (m^-3)\tNe (m^-3)\n"
		<<Pdata->NeutralDensity<<"\t\t"<<Pdata->IonDensity<<"\t\t"<<Pdata->ElectronDensity
		<<"\n\nTn (K)\t\tTi (K)\t\tTe (K)\n"
		<<Pdata->NeutralTemp<<"\t\t"<<Pdata->IonTemp<<"\t\t"<<Pdata->ElectronTemp<<"\t"
		<<"\n\nPvel (m s^-1)\t\tE (V m^-1)\t\tB (T)\n"
		<<PlasmaVelocity<<"\t"<<Efield<<"\t"<<Bfield<<"\n"
		<<"\n\n##MODEL SWITHES\n#HEATING MODELS\n"
		<<"RadiativeCooling:\t" << HeatModels[0] << "\n"<<"EvaporativeCooling:\t" << HeatModels[1] << "\n"
		<<"NewtonCooling:\t\t" << HeatModels[2] << "\n"<<"NeutralHeatFlux:\t" << HeatModels[3] << "\n"
		<<"ElectronHeatFlux:\t" << HeatModels[4] << "\n"
		<<"IonHeatFlux:\t\t" << HeatModels[5] << "\n"<<"NeutralRecomb:\t\t" << HeatModels[6] << "\n"
		<<"TEE:\t\t\t" << HeatModels[7] << "\n"<<"SEE:\t\t\t" << HeatModels[8] << "\n"
		<<"\n#FORCING MODELS\n"
		<<"Gravity:\t\t" << ForceModels[0] << "\n"<<"Centrifugal:\t\t" << ForceModels[1] << "\n"
		<<"Lorentz:\t\t" << ForceModels[2] << "\n"<<"IonDrag:\t\t" << ForceModels[3] << "\n"
		<<"HybridDrag:\t\t" << ForceModels[4] << "\n"<<"NeutralDrag:\t\t" << ForceModels[5] << "\n"
		<<"\n#CHARGING MODELS\n"
		<<"DTOKSOML:\t\t" << ChargeModels[0] << "\n"<<"SchottkyOML:\t\t" << ChargeModels[1] << "\n" << "DTOKSWell:\t\t" << ChargeModels[2];


	std::cout << "\n\n * GENERATE DTOKS * \n\n";
//	DTOKSU *MyDtoks1 = new DTOKSU(AccuracyLevels, Sample, Pdata, HeatModels, ForceModels, ChargeModels);
	DTOKSU *MyDtoks2 = new DTOKSU(AccuracyLevels, Sample, Pgrid, HeatModels, ForceModels, ChargeModels);

	std::cout << "\n\n * RUN DTOKS * \n\n";
	Breakup Break(MyDtoks2, Sample);
	Break.Run();
//	MyDtoks2->Run();

	std::cout << "\n\n * MAIN SCRIPT COMPLETE * \n\n";
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;	
	MetaDataFile << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	MetaDataFile.close();
	cfg->destroy();
//	Pgrid.datadump(); // Print the plasma grid data
	return 0;
}
