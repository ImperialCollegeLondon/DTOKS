//#include "DTOKSU.h"
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
	<< "\t-s,--size SIZE\t\t\tdouble variable giving the radius of the grain\n\n"
	<< "\t-vr,--rvel RVEL\t\t\tfloat variable defining radial velocity\n\n"
	<< "\t-vt,--thetavel THETAVEL\tfloat variable defining angular velocity\n\n"
	<< "\t-vz,--zvel ZVEL\t\t\tfloat variable defining lognitudinal velocity\n\n"
	<< "\t-rr,--rpos RPOS\t\t\tfloat variable defining radial position\n\n"
	<< "\t-rt,--thetapos THETAPOS\tfloat variable defining angular position\n\n"
	<< "\t-rz,--zpos ZPOS\t\t\tfloat variable defining longitudinal position\n\n"
	<< "\t-op,--output OUTPUT\t\tstring variable defining the filename prefix to write to\n\n"
	<< "\t-om,--metadata METADATA\tstring variable defining the MetaData filename to write to\n\n";
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
	std::string MetaDataFilename = "Data/DTOKSU.txt";
	std::string DataFilePrefix = "Data/DTOKSU";
//	std::string Config_Filename = "Config_Files/DTOKSU_Config_Magnum-PSI.cfg";
//	std::string Config_Filename = "Config_Files/DTOKSU_Config_MAST.cfg";
	std::string Config_Filename = "Config_Files/DTOKSU_Config_JET.cfg";
	
	std::vector <std::string> sources;
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 		|| arg == "-h" ){	show_usage( argv[0]); return 0; 		}
		else if( arg == "--config"      || arg == "-c" )	InputFunction(argc,argv,i,ss0,Config_Filename);
		else{
			sources.push_back(argv[i]);
		}
	}

	

	// Coordinates are r, theta, z.
	// z is the vertical direction.
	// r and theta map out the toroidal plane.
	// Therefore a plot in r, z provides a poloidal cross section
//	threevector xinit(0.147,0.01,0.1575);	// MAGNUM-PSI: Coordinates are r,theta,z
//	threevector xinit(1.15,0.0,-1.99);// default injection right hand side
//	threevector vinit(-1.4,0.0,0.0);		// MAGNUM-PSI
//	threevector vinit(0.0,0.0,100.0);

	std::string plasma_file = "Models/PlasmaData/MagnumPSI/Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9.nc";
	char Plasma='h';
	char Machine='p';	// Initialise plasma grid
	float xSpacing = 0.00234375;	// 0.15/64.0; x-dimensional spacing (Radial) metres
//	float zSpacing = 0.04; 	// 1.0/25.0;  y-dimensional spacing (z) metres, Magnum-PSI_Prelim_B1.41_L1.0.nc
	float zSpacing = 0.095; // 1.9/20.0;  y-dimensional spacing (z) metres, Magnum-PSI_Experiment_Homogeneous-B-Field_B0.1_L1.9.nc

	char Element='W';
	float size=0.5e-6;
	float Temp=300;
	float InitRotationalFreq(0.0);
	float rpos(0.147);
	float thetapos(0.01);
	float zpos(0.1575);
	float rvel(-1.4);
	float thetavel(0.0);
	float zvel(0.0);

	char EmissivityModel = 'v';
	char ExpansionModel = 'c';
	char HeatCapacityModel = 'v'; 
	char BoilingModel = 't';

	// ------------------- GROUP MODELS ------------------- //
	std::array<bool, 9> HeatModels;
	std::array<bool,6> ForceModels;
	std::array<bool,4> ChargeModels;
	std::array<char,4> ConstModels;
	std::array<float,3> AccuracyLevels;
	PlasmaData *Pdata = new PlasmaData;

	// ------------------- PROCESS USER-INPUT ------------------- //

	threevector PlasmaVelocity(501.33, 7268.5, 914.947); // Taken from initial for DTOKS
	threevector Efield(-13.673, 0, -27.925);
	threevector Bfield(0.0226868, 0.328923, 0.0414043);
	Pdata->PlasmaVel 		= PlasmaVelocity;
	Pdata->ElectricField 	= Efield;
	Pdata->MagneticField 	= Bfield;
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	std::cout << "\n* Reading configuration file: " << Config_Filename << " *";
	try {
        cfg->parse(Config_Filename.c_str());
        MetaDataFilename = cfg->lookupString("", "Filename");
        DataFilePrefix = cfg->lookupString("", "DataFilePrefix");
        plasma_file = cfg->lookupString("plasmagrid", "Filename");
        Plasma   = cfg->lookupString("plasmagrid", "Plasma")[0];
        Machine  = cfg->lookupString("plasmagrid", "Machine")[0];
        xSpacing = cfg->lookupFloat("plasmagrid", "xSpacing");
        zSpacing = cfg->lookupFloat("plasmagrid", "zSpacing");
		Element  = cfg->lookupString("dust", "Element")[0];
        size     = cfg->lookupFloat("dust", "size");
        Temp     = cfg->lookupFloat("dust", "Temp");
		rpos = cfg->lookupFloat("dust", "dynamics.rpos");
		thetapos = cfg->lookupFloat("dust", "dynamics.thetapos");
		zpos = cfg->lookupFloat("dust", "dynamics.zpos");
		rvel = cfg->lookupFloat("dust", "dynamics.rvel");
		thetavel = cfg->lookupFloat("dust", "dynamics.thetavel");
		zvel = cfg->lookupFloat("dust", "dynamics.zvel");
		InitRotationalFreq = cfg->lookupFloat("dust", "dynamics.InitRotationalFreq");
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
				cfg->lookupBoolean("chargemodels","DTOKSWell"), cfg->lookupBoolean("chargemodels","PHL")
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
    std::cout << "\n* Configuration file read successfully! *\n\n* Creating PlasmaGrid Object *";
    std::cout << "\n\t* Plasma:\t" << Plasma << "\n\t* Machine:\t" << Machine;
    std::cout << "\n\t* xSpacing:\t" << xSpacing << "\n\t* zSpacing:\t" << zSpacing;

    PlasmaGrid Pgrid(plasma_file,Plasma,Machine,xSpacing,zSpacing);
    std::cout << "\n* PlasmaGrid created successfully! *\n\n* Processing command line input *";

// ------------------- PROCESS USER-INPUT ------------------- //
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if( arg == "--temperature" || arg == "-t"  )	InputFunction(argc,argv,i,ss0,Temp);
		else if( arg == "--material" 	|| arg == "-m"  ) 	InputFunction(argc,argv,i,ss0,Element);
		else if( arg == "--size" 		|| arg == "-s"  )	InputFunction(argc,argv,i,ss0,size);
		else if( arg == "--rvel" 		|| arg == "-vr"  )	InputFunction(argc,argv,i,ss0,rvel);
		else if( arg == "--thetavel"    || arg == "-vt"  )	InputFunction(argc,argv,i,ss0,thetavel);
		else if( arg == "--zvel" 		|| arg == "-vz"  )	InputFunction(argc,argv,i,ss0,zvel);
		else if( arg == "--rpos" 		|| arg == "-rr" )	InputFunction(argc,argv,i,ss0,rpos);
		else if( arg == "--thetapos"    || arg == "-rt"  )	InputFunction(argc,argv,i,ss0,thetapos);
		else if( arg == "--zpos" 		|| arg == "-rz"  )	InputFunction(argc,argv,i,ss0,zpos);
		else if( arg == "--output" 		|| arg == "-op"  )	InputFunction(argc,argv,i,ss0,DataFilePrefix);
		else if( arg == "--MetaData" 		|| arg == "-om"  )	InputFunction(argc,argv,i,ss0,MetaDataFilename);
		else{
			sources.push_back(argv[i]);
		}
	}
	std::cout << "\n* Command line input processed successfully! *\n\n";

    // ------------------- INITIALISE META_DATA FILE ------------------- //
    std::cout << "* Creating MetaDataFile: " << MetaDataFilename << " *\n\n";
    std::ofstream MetaDataFile;	// Data file for containing the run information
    time_t now = time(0);		// Get the time of simulation
	char * dt = ctime(&now);
    MetaDataFile.open(MetaDataFilename);
	MetaDataFile << "## Run Data File ##\n";
	MetaDataFile << "#Date:\t" << dt;


	// ------------------- INITIALISE DUST ------------------- //
	std::cout << "* Creating Matter object *\n\t* Element:\t" << Element;
	Matter * Sample;		// Define the sample matter type
	if 	(Element == 'W')     Sample = new Tungsten(size,Temp,ConstModels);
	else if (Element == 'B') Sample = new Beryllium(size,Temp,ConstModels);
	else if (Element == 'F') Sample = new Iron(size,Temp,ConstModels);
	else if (Element == 'G') Sample = new Graphite(size,Temp,ConstModels);
	else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}
	threevector xinit(rpos,thetapos,zpos);
	threevector vinit(rvel,thetavel,zvel);
	Sample->update_motion(xinit,vinit,InitRotationalFreq);
	std::cout << "\n\t* xinit:\t" << xinit << "\n\t* vinit:\t" << vinit;
	std::cout << "\n\t* InitRotation:\t" << InitRotationalFreq << "\n";

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
		<<"DTOKSOML:\t\t" << ChargeModels[0] << "\n"<<"SchottkyOML:\t\t" << ChargeModels[1] << "\n" 
		<< "DTOKSWell:\t\t" << ChargeModels[2] << "\n" << "PHL:\t\t" << ChargeModels[3];


	std::cout << "\n\n * CREATING DTOKS * \n";
//	DTOKSU *MyDtoks1 = new DTOKSU(AccuracyLevels, Sample, Pdata, HeatModels, ForceModels, ChargeModels);
	DTOKSU *MyDtoks2 = new DTOKSU(AccuracyLevels, Sample, Pgrid, HeatModels, ForceModels, ChargeModels);

	MyDtoks2->OpenFiles(DataFilePrefix,0);
	std::cout << "\n * DTOKS SUCCESSFULLY INITIALISED * \n";

	std::cout << "\n * RUNNING DTOKS * \n";
	Breakup Break(MyDtoks2, Sample);
	Break.Run();
//	MyDtoks2->Run();
//	MyDtoks2->CloseFiles();

	std::cout << "\n\n * DTOKS COMPLETED SUCCESSFULLY * \n\n";
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;	
	MetaDataFile << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	MetaDataFile.close();
	cfg->destroy();
//	Pgrid.datadump(); // Print the plasma grid data
	return 0;
}
