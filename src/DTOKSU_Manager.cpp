//#define PAUSE
//#define DM_DEBUG
#include "DTOKSU_Manager.h"

// Default constructor
DTOKSU_Manager::DTOKSU_Manager(){
	DM_Debug("\n\nIn DTOKSU_Manager::DTOKSU_Manager()\n\n");
	Config_Status = -1;
};

// Parameterised constructor, from command line.
// Call the configure function with command line options to configure as well as construct.
DTOKSU_Manager::DTOKSU_Manager(int argc, char* argv[]){
	DM_Debug("\n\nIn DTOKSU_Manager::DTOKSU_Manager(int argc, char* argv[])\n\n");
	std::cout << "\n * CONFIGURING DTOKS * \n";
	Config_Status = -1;
	Config_Status = Configure(argc,argv);
};

// Parameterised constructor, from command line and config_file.
// Call the configure function with command line options to configure as well as construct.
DTOKSU_Manager::DTOKSU_Manager(int argc, char* argv[], std::string filename){
	DM_Debug("\n\nIn DTOKSU_Manager::DTOKSU_Manager(int argc, char* argv[], std::string filename = Config/DTOKSU_Config.cfg)\n\n");
	std::cout << "\n * CONFIGURING DTOKS * \n";
	Config_Status = -1;
	Config_Status = Configure(argc,argv,filename);
};


void DTOKSU_Manager::config_message()const{
	if( Config_Status == -1 ){ // Configuration hasn't been processed
		std::cout << "\n * DTOKS UNCONFIGURED * \n";
	}if( Config_Status == 0 ){ // Configuration was processed successfully
		std::cout << "\n * DTOKS SUCCESSFULLY CONFIGURED * \n";
	}else if( Config_Status == 1 ){	// Help message was displayed
		std::cout << "\n\n * HELP MESSAGE DISPLAYED * \n\n";
	}else if( Config_Status == 2 ){ // Failed to read configuration file
		std::cout << "\n\n * ERROR CODE 2! FAILURE PARSING CONFIGURATION FILE * \n\n";
	}else if( Config_Status == 3 ){ // Failed to configure plasma data
		std::cout << "\n\n * ERROR CODE 3! FAILURE CONFIGURING PLASMA DATA * \n\n";
	}else if( Config_Status == 4 ){ // Invalid element option selected
		std::cout << "\n\n * ERROR CODE 4! FAILURE CONFIGURING SAMPLE * \n\n";
	}else{
		std::cout << "\n\n * UNKNOWN CONFIGURATION STATUS * \n\n";
	}
}

// Remind user of the command line options
void DTOKSU_Manager::show_usage(std::string name)const{
	std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
	<< "\n\nOptions:\n"
	<< "\t-h, --help\t\tShow this help message\n\n"
	<< "\t-c  --config CONFIG\t\tSpecify the configuration file path\n\n"
	<< "\t-t, --temperature TEMPERATURE\tdouble variable defining the initial temperature\n\n"
	<< "\t-m, --material MATERIAL\t\tchar variable giving the dust grain element, possible values 'w', 'g', 'b' and 'f'\n"
	<< "\t\t\t\t\t(W): Tungsten, (G): Graphite, (B): Beryllium or (F): Iron\n\n"
	<< "\t-s, --size SIZE\t\t\tdouble variable giving the radius of the grain\n\n"
	<< "\t-vr,--rvel RVEL\t\t\tfloat variable defining radial velocity\n\n"
	<< "\t-vt,--thetavel THETAVEL\t\tfloat variable defining angular velocity\n\n"
	<< "\t-vz,--zvel ZVEL\t\t\tfloat variable defining lognitudinal velocity\n\n"
	<< "\t-rr,--rpos RPOS\t\t\tfloat variable defining radial position\n\n"
	<< "\t-rt,--thetapos THETAPOS\t\tfloat variable defining angular position\n\n"
	<< "\t-rz,--zpos ZPOS\t\t\tfloat variable defining longitudinal position\n\n"
	<< "\t-op,--output OUTPUT\t\tstring variable defining the filename prefix to write to\n\n"
	<< "\t-om,--metadata METADATA\t\tstring variable defining the MetaData filename to write to\n\n";
}

// Process input from user into variable Temp of type T
template<typename T> int DTOKSU_Manager::input_function(int &argc, char* argv[], int &i, std::stringstream &ss0, T &Temp)const{
	if (i + 1 < argc) { // Make sure we aren't at the end of argv!
		i+=1;
		ss0 << argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
		ss0 >> Temp;
		ss0.clear(); ss0.str("");
		return 0;
	}else{ // Uh-oh, there was no argument to the destination option.f
		std::cerr << "\noption requires argument." << std::endl;
		return 1;
	}
}

// Configure DTOKS ready to run with a particular plasma and matter sample
// which are established from configuration filename
int DTOKSU_Manager::Configure(std::string Config_Filename){
	char* argv[0] = {};
	return Configure(0,argv,Config_Filename);
}

// Configure DTOKS ready to run with a particular plasma and matter sample
// which are established from command line options and DTOKSU_Config.cfg configuration file
int DTOKSU_Manager::Configure(int argc, char* argv[], std::string Config_Filename){
	DM_Debug("\n\nIn DTOKSU_Manager::Configure(int argc, char* argv[], std::string Config_Filename=Config_Files/DTOKSU_Config.cfg)\n\n");

	// ------------------- PARSE CONFIGURATION FILE ------------------- //
	std::string MetaDataFilename = "Data/DTOKSU.txt";
	std::string DataFilePrefix = "Data/DTOKSU";
	std::string dir_name = "PlasmaData/Magnum-PSI/";

	// Check user input for specifying the configuration file.
	// We want other command line options to over-ride the configuration file.
	std::vector <std::string> sources;
	std::stringstream ss0;
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if     ( arg == "--help" 		|| arg == "-h" ){	
			Config_Status = 1; show_usage( argv[0] ); return Config_Status; 		
		}else if( arg == "--config"      || arg == "-c" ){
			input_function(argc,argv,i,ss0,Config_Filename);
		}else{
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

	// ------------------- DUST VARIABLE DEFAULTS ------------------- //
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
	bool ContinuousPlasma = false;

	// ------------------- GROUP MODELS ------------------- //
	std::array<bool,HMN> HeatModels;
	std::array<bool,FMN> ForceModels;
	std::array<bool,CMN> ChargeModels;
	std::array<char,4> ConstModels;
	std::array<float,DTOKSU::MN> AccuracyLevels;

	// ------------------- PROCESS CONFIGURATION FILE ------------------- //
//	threevector PlasmaVelocity(501.33, 7268.5, 914.947); // Taken from initial for DTOKS
//	threevector Efield(-13.673, 0, -27.925);
//	threevector Bfield(0.0226868, 0.328923, 0.0414043);
	config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
	std::cout << "\n* Reading configuration file: " << Config_Filename << " *";
	try {
        cfg->parse(Config_Filename.c_str());
        MetaDataFilename 		= cfg->lookupString("", "Filename");
        DataFilePrefix 			= cfg->lookupString("", "DataFilePrefix");
        ContinuousPlasma		= cfg->lookupBoolean("plasma", "ContinuousPlasma");
        dir_name 				= cfg->lookupString("plasma", "plasmagrid.Dirname");
        Pgrid.gas   			= cfg->lookupString("plasma", "plasmagrid.Plasma")[0];
        Pgrid.device  			= cfg->lookupString("plasma", "plasmagrid.Machine")[0];
        Pgrid.dlx 				= cfg->lookupFloat("plasma", "plasmagrid.xSpacing");
        Pgrid.dlz	 			= cfg->lookupFloat("plasma", "plasmagrid.zSpacing");
        Pdata.IonDensity 		= cfg->lookupFloat("plasma","plasmadata.IonDensity");
		Pdata.ElectronDensity 	= cfg->lookupFloat("plasma","plasmadata.ElectronDensity");
		Pdata.NeutralDensity 	= cfg->lookupFloat("plasma","plasmadata.NeutralDensity");
		Pdata.IonTemp 			= cfg->lookupFloat("plasma","plasmadata.IonTemp");
		Pdata.NeutralTemp 		= cfg->lookupFloat("plasma","plasmadata.NeutralTemp");
		Pdata.ElectronTemp 		= cfg->lookupFloat("plasma","plasmadata.ElectronTemp");
		Pdata.AmbientTemp 		= cfg->lookupFloat("plasma","plasmadata.AmbientTemp");
		Pdata.mi 				= cfg->lookupFloat("plasma","plasmadata.Mi");
		Element  				= cfg->lookupString("dust", "Element")[0];
        size     				= cfg->lookupFloat("dust", "size");
        Temp     				= cfg->lookupFloat("dust", "Temp");
		rpos 					= cfg->lookupFloat("dust", "dynamics.rpos");
		thetapos 				= cfg->lookupFloat("dust", "dynamics.thetapos");
		zpos 					= cfg->lookupFloat("dust", "dynamics.zpos");
		rvel 					= cfg->lookupFloat("dust", "dynamics.rvel");
		thetavel 				= cfg->lookupFloat("dust", "dynamics.thetavel");
		zvel 					= cfg->lookupFloat("dust", "dynamics.zvel");
		InitRotationalFreq 		= cfg->lookupFloat("dust", "dynamics.InitRotationalFreq");
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
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        cfg->destroy();
        Config_Status = 2;
        return Config_Status;
    }
    cfg->destroy();
    std::cout << "\n* Configuration file read successfully! *";
    std::cout << "\nContinuousPlasma = " << ContinuousPlasma;

	// ------------------- CONFIGURE PLASMAGRID ------------------- //
	if( !ContinuousPlasma ){
		std::cout << "\n\n* Full Machine Simulation! *\n\n* Creating PlasmaGrid_Data Structure *";
    	std::cout << "\n\t* Plasma:\t" << Pgrid.gas << "\n\t* Machine:\t" << Pgrid.device;
    	std::cout << "\n\t* xSpacing:\t" << Pgrid.dlx << "\n\t* zSpacing:\t" << Pgrid.dlz;
    	if( configure_plasmagrid(dir_name) != 0 ){ // Failed to configure plasma data
    		std::cerr << "\nFailed to configure plasma data!";
    		Config_Status = 3;
    		return Config_Status;
    	}
		std::cout << "\n* PlasmaGrid_Data Structure created successfully! *\n\n* Processing command line input *";
    }else{
    	std::cout << "\n\n* ContinuousPlasma! *";
    	std::cout << "\n* PlasmaData Structure created successfully! *\n\n* Processing command line input *";
    }



	// ------------------- PROCESS USER-INPUT ------------------- //
	for (int i = 1; i < argc; ++i){ // Read command line input
		std::string arg = argv[i];
		if( 	 arg == "--temperature" || arg == "-t"   )	input_function(argc,argv,i,ss0,Temp);
		else if( arg == "--material" 	|| arg == "-m"   ) 	input_function(argc,argv,i,ss0,Element);
		else if( arg == "--size" 		|| arg == "-s"   )	input_function(argc,argv,i,ss0,size);
		else if( arg == "--rvel" 		|| arg == "-vr"  )	input_function(argc,argv,i,ss0,rvel);
		else if( arg == "--thetavel"    || arg == "-vt"  )	input_function(argc,argv,i,ss0,thetavel);
		else if( arg == "--zvel" 		|| arg == "-vz"  )	input_function(argc,argv,i,ss0,zvel);
		else if( arg == "--rpos" 		|| arg == "-rr"  )	input_function(argc,argv,i,ss0,rpos);
		else if( arg == "--thetapos"    || arg == "-rt"  )	input_function(argc,argv,i,ss0,thetapos);
		else if( arg == "--zpos" 		|| arg == "-rz"  )	input_function(argc,argv,i,ss0,zpos);
		else if( arg == "--output" 		|| arg == "-op"  )	input_function(argc,argv,i,ss0,DataFilePrefix);
		else if( arg == "--MetaData" 	|| arg == "-om"  )	input_function(argc,argv,i,ss0,MetaDataFilename);
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
	threevector xinit(rpos,thetapos,zpos);
	threevector vinit(rvel,thetavel,zvel);
	std::cout << "* Creating Matter object *\n\t* Element:\t" << Element;
	if 	(Element == 'W')     Sample = new Tungsten(size,Temp,ConstModels,xinit,vinit);
	else if (Element == 'B') Sample = new Beryllium(size,Temp,ConstModels,xinit,vinit);
	else if (Element == 'F') Sample = new Iron(size,Temp,ConstModels,xinit,vinit);
	else if (Element == 'G') Sample = new Graphite(size,Temp,ConstModels,xinit,vinit);
	else{ 
		std::cerr << "\nInvalid Option entered for Element";
		Config_Status = 4;
		return Config_Status;
	}
	std::cout << "\n\t* xinit:\t" << Sample->get_position() << "\n\t* vinit:\t" << Sample->get_velocity();
	std::cout << "\n\t* InitRotation:\t" << Sample->get_rotationalfreq() << "\n";

	// ------------------- PRINT METADATA / CREATE DTOKSU ------------------- //
	MetaDataFile << "\n\n#DUST PARAMETERS" 
		<<"\nElem (arb)\tRadius (m)\tTemp (K)\txinit (m s^-1)\t\tvinit (m s^-1)\n"
		<<Sample->get_elem()<<"\t\t"<<Sample->get_radius()<<"\t\t"<<Sample->get_temperature()
		<<"\t\t"<<Sample->get_position()<<"\t\t"<<Sample->get_velocity()<<"\n";

	if( ContinuousPlasma ){
		MetaDataFile <<"\n\n#PLASMA PARAMETERS"
			<<"\nMachine\tgas\tgridx\tgridz\tgridtheta\tdlx\tdlz\n"
			<<Pgrid.device<<"\t"<<Pgrid.gas<<"\t"<<Pgrid.gridx<<"\t"<<Pgrid.gridz<<"\t"
			<<Pgrid.gridtheta<<"\t"<<Pgrid.dlx<<"\t"<<Pgrid.dlz<<"\n"
			<<"\nxmin (m)\txmax (m)\tzmin (m)\tzmax (m)\n"
			<<Pgrid.gridxmin<<"\t\t"<<Pgrid.gridxmax<<"\t\t"<<Pgrid.gridzmin<<"\t\t"<<Pgrid.gridzmax << "\n";
		Sim = new DTOKSU(AccuracyLevels, Sample, Pdata, HeatModels, ForceModels, ChargeModels);
	}else if( !ContinuousPlasma ){
		MetaDataFile <<"\n\nNn (m^-3)\tNi (m^-3)\tNe (m^-3)\n"
			<<Pdata.NeutralDensity<<"\t\t"<<Pdata.IonDensity<<"\t\t"<<Pdata.ElectronDensity
			<<"\n\nTn (K)\t\tTi (K)\t\tTe (K)\n"
			<<Pdata.NeutralTemp<<"\t\t"<<Pdata.IonTemp<<"\t\t"<<Pdata.ElectronTemp<<"\t"
			<<"\n\nPvel (m s^-1)\t\tE (V m^-1)\t\tB (T)\n"
			<<Pdata.PlasmaVel<<"\t"<<Pdata.ElectricField<<"\t"<<Pdata.MagneticField<<"\n";
		Sim = new DTOKSU(AccuracyLevels, Sample, Pgrid, HeatModels, ForceModels, ChargeModels);
	}
	MetaDataFile <<"\n\n##MODEL SWITHES\n#HEATING MODELS\n"
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
		<< "DTOKSWell:\t\t" << ChargeModels[2] << "\n" << "PHL:\t\t\t" << ChargeModels[3];
	MetaDataFile.close();

	Sim->OpenFiles(DataFilePrefix,0);
	Config_Status = 0;
	config_message();
	return Config_Status;
}

// Function to configure plasma grid ready for simulation. 
// Plasma data is read from plasma_dirname+filename where filename is a hard-coded string
int DTOKSU_Manager::configure_plasmagrid(std::string plasma_dirname){
	DM_Debug("\n\nIn DTOKSU_Manager::configure_plasmagrid(std::string plasma_dirname)\n\n");
	// Plasma parameters
	if(Pgrid.device=='m'){
		Pgrid.gridx = 121;
		Pgrid.gridz = 281;
		Pgrid.gridtheta=0;
		Pgrid.gridxmin = 0.2;
		Pgrid.gridzmin = -2.0;
		Pgrid.gridxmax = 1.4;
		Pgrid.gridzmax = 0.80;
		std::cout << "\n\n\tCalculation for MAST" << std::endl;
	}else if(Pgrid.device=='i'){
		Pgrid.gridx = 451;
		Pgrid.gridz = 951;
		Pgrid.gridtheta=0;
		Pgrid.gridxmin = 4.0;
		Pgrid.gridzmin = -4.7;
		Pgrid.gridxmax = 8.25;
		Pgrid.gridzmax = 4.80;
		std::cout << "\n\n\tCalculation for ITER" << std::endl;
	}else if(Pgrid.device=='j'){
		Pgrid.gridx = 251;
		Pgrid.gridz = 401;
		Pgrid.gridtheta=0;
		Pgrid.gridxmin = 1.5;
		Pgrid.gridzmin = -2.0;
		Pgrid.gridxmax = 4.0;
		Pgrid.gridzmax = 2.0;
		std::cout << "\n\n\tCalculation for JET" << std::endl;
	}else if(Pgrid.device=='d'){
		Pgrid.gridx = 180;
		Pgrid.gridz = 400;
		Pgrid.gridtheta=0;
		Pgrid.gridxmin = 0.2;
		Pgrid.gridzmin = -2.0;
		Pgrid.gridxmax = 2.0;
		Pgrid.gridzmax = 4.0;
		std::cout << "\n\n\tCalculation for Double Null MAST 17839files shot" << std::endl;
	}else if(Pgrid.device=='p'){
		Pgrid.gridx = 64;
		Pgrid.gridz = 20;
		Pgrid.gridtheta=64;
		Pgrid.gridxmin = 0.0;
		Pgrid.gridzmin = 0.0;
		Pgrid.gridxmax = 0.15;
		Pgrid.gridzmax = 1.9;
//		Pgrid.gridzmax = 1.0;
		std::cout << "\n\n\t* Calculation for Magnum-PSI *\n";
	}else{ 
		std::cout << "Invalid tokamak" << std::endl;
	}

	if( Pgrid.gas == 'h' ){ // For Hydrogen plasma, we have hydrogen ions of mass Mp
		Pgrid.mi = Mp;
	}else{
		std::cout << "Invalid gas, Pgrid.gas = " << Pgrid.gas << "!\n";
		return 1;
	}

	//impurity.open("output///impurity///impurity.vtk");
	int readstatus(-1);
	try{ // Read data
		readstatus = read_data(plasma_dirname);
		if(readstatus != 0)
			throw readstatus;
	}	
	catch( int e ){
		std::cerr << "Exception opening/reading/closing file\nError status: " << e;
		return readstatus;
	}
	return 0;
}

int DTOKSU_Manager::read_data(std::string plasma_dirname){
	P_Debug("\tIn DTOKSU_Manager::read_data(std::string plasma_dirname)\n\n");
	// 
	if(Pgrid.device=='p'){ // Note, grid flags will be empty 
		return read_MPSIdata(plasma_dirname);
	}else{
		std::ifstream scalars,threevectors,gridflagfile;
		if(Pgrid.device=='m'){
			scalars.open(plasma_dirname+"b2processed.dat");
			threevectors.open(plasma_dirname+"b2processed2.dat");
			gridflagfile.open(plasma_dirname+"locate.dat");
		}else if(Pgrid.device=='i'){
			scalars.open(plasma_dirname+"b2processed.dat");
			threevectors.open(plasma_dirname+"b2processed2.dat");
			gridflagfile.open(plasma_dirname+"locate.dat");
		}else if(Pgrid.device=='j'){
			scalars.open(plasma_dirname+"Interpolated_Data85974PostBolomValid_PlasmaData.txt");
			threevectors.open(plasma_dirname+"Interpolated_Data85974PostBolomValid_BFieldData.txt");
			gridflagfile.open(plasma_dirname+"locate.dat");
		}else if(Pgrid.device=='d'){
			scalars.open(plasma_dirname+"b2processed.dat");
			threevectors.open(plasma_dirname+"b2processed2.dat");
			gridflagfile.open(plasma_dirname+"locate.dat");
		}else{
			std::cerr << "\nInvalid device char!";
			return 1;
		}

		// Open files to read
		assert( scalars.is_open() );
		assert( threevectors.is_open() );
		assert( gridflagfile.is_open() );

		// Preallocate size of vectors
		Pgrid.Te 		= std::vector<std::vector<double>>(Pgrid.gridx,std::vector<double>(Pgrid.gridz));
		Pgrid.Ti 		= Pgrid.Te;
		Pgrid.na0		= Pgrid.Te;
		Pgrid.na1		= Pgrid.Te;
		Pgrid.po 		= Pgrid.Te;
		Pgrid.ua0		= Pgrid.Te;
		Pgrid.ua1		= Pgrid.Te;
		Pgrid.bx 		= Pgrid.Te;
		Pgrid.by 		= Pgrid.Te;
		Pgrid.bz 		= Pgrid.Te;
		Pgrid.x  		= Pgrid.Te;
		Pgrid.z  		= Pgrid.Te;
		Pgrid.gridflag	= std::vector<std::vector<int>>(Pgrid.gridx,std::vector<int>(Pgrid.gridz));

		// Throw away variables which read in data which is unimportnat.
		char dummy_char;
		double dummy_dub;
		// Ignore first line of file
		for(unsigned int i=0; i<=19; i++){
			scalars >> dummy_char;
			threevectors >> dummy_char;
		}
		// Now actually loop over the grid and feed in the data into the vectors
		for(unsigned int k=0; k<=Pgrid.gridz-1; k++){
			for(unsigned int i=0; i<=Pgrid.gridx-1; i++){
				scalars >> Pgrid.x[i][k] >> Pgrid.z[i][k] >> Pgrid.Te[i][k] >> Pgrid.Ti[i][k] 
						>> Pgrid.na0[i][k] >> Pgrid.na1[i][k] >> Pgrid.po[i][k] >> Pgrid.ua0[i][k]
						>> Pgrid.ua1[i][k];
				threevectors >> dummy_dub >> dummy_dub >> Pgrid.bx[i][k] >> Pgrid.bz[i][k] >> Pgrid.by[i][k];
				gridflagfile >> dummy_dub >> dummy_dub >> Pgrid.gridflag[i][k];
				// For some reason, very small non-zero values are being assigned to the 'zero' values being read in.
				// This should correct for this...

				if( Pgrid.device != 'j' )								{	Pgrid.bz[i][k]  = -Pgrid.bz[i][k]; }
				if( Pgrid.bx[i][k] != 0 && Pgrid.bx[i][k] < 1e-100 ) 	{	Pgrid.bx[i][k]  = 0.0; } 
				if( Pgrid.by[i][k] != 0 && Pgrid.by[i][k] < 1e-100 ) 	{	Pgrid.by[i][k]  = 0.0; } 
				if( Pgrid.bz[i][k] != 0 && Pgrid.bz[i][k] < 1e-100 ) 	{	Pgrid.bz[i][k]  = 0.0; }
				if( Pgrid.Te[i][k] != 0 && Pgrid.Te[i][k] < 1e-100 ) 	{	Pgrid.Te[i][k]  = 0.0; } 
				if( Pgrid.Ti[i][k] != 0 && Pgrid.Ti[i][k] < 1e-100 ) 	{	Pgrid.Ti[i][k]  = 0.0; } 
				if( Pgrid.na0[i][k]!= 0 && Pgrid.na0[i][k] < 1e-100 )	{	Pgrid.na0[i][k] = 0.0; }
				if( Pgrid.na1[i][k]!= 0 && Pgrid.na1[i][k] < 1e-100 )	{	Pgrid.na1[i][k] = 0.0; }
				if( Pgrid.po[i][k] != 0 && Pgrid.po[i][k] < 1e-100 ) 	{	Pgrid.po[i][k]  = 0.0; } 
				if( Pgrid.ua0[i][k]!= 0 && Pgrid.ua0[i][k] < 1e-100 )	{	Pgrid.ua0[i][k] = 0.0; }
				if( Pgrid.ua1[i][k]!= 0 && Pgrid.ua1[i][k] < 1e-100 )	{ 	Pgrid.ua1[i][k] = 0.0; }
			}
		}
		scalars.close();
		threevectors.close();
		gridflagfile.close();
	}
	return 0; // return success!
}
// *************************************** READING FUNCTIONS *************************************** //

// for Magnum-PSI, we need to read a NET-cdf file which is special
// This function does all the necessary effort of extracting this information.
int DTOKSU_Manager::read_MPSIdata(std::string plasma_dirname){
	P_Debug("\tDTOKSU_Manager::read_MPSIdata(std::string plasma_dirname)\n\n");
	
	const int NC_ERR = 2;
	std::string filename = "Magnum-PSI_Experiment_Homogeneous-B-Field_B0.2_L1.9.nc";
	std::cout << "\n\t\t* Reading data from: " << plasma_dirname + filename << " *";

	float electron_dens_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
	float electron_temp_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
	float electron_Vele_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
	float ion_dens_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
	float ion_Veli_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
	float Potential_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
	float Bxy_mat[Pgrid.gridx][Pgrid.gridz];

	// Change the error behavior of the netCDF C++ API by creating an
	// NcError object. Until it is destroyed, this NcError object will
	// ensure that the netCDF C++ API silently returns error codes on
	// any failure, and leaves any other error handling to the calling
	// program. In the case of this example, we just exit with an
	// NC_ERR error code.
	NcError err(NcError::silent_nonfatal);

	
	NcFile dataFile((plasma_dirname+filename).c_str(), NcFile::ReadOnly);
	if(!dataFile.is_valid())
		return NC_ERR;
	
	if (dataFile.num_dims() != 3 || dataFile.num_vars() != 8 ||
		dataFile.num_atts() != 0 || dataFile.rec_dim() != 0)
		return NC_ERR;
	// We get back a pointer to each NcVar we request. Get the
	// latitude and longitude coordinate variables.

	NcVar *Ne, *e_Temp, *Vel_e, *Ni, *Vel_i, *Potential; //, *B_xy;
	if (!(Ne = dataFile.get_var("Ne")))
		return NC_ERR;
	if (!(e_Temp = dataFile.get_var("Te")))
		return NC_ERR;
	if (!(Vel_e = dataFile.get_var("Vel_e")))
		return NC_ERR;
	if (!(Ni = dataFile.get_var("Ni")))
		return NC_ERR;
	if (!(Vel_i = dataFile.get_var("Vel_i")))
		return NC_ERR;
	if (!(Potential = dataFile.get_var("Potential_Phi")))
		return NC_ERR;
//	if (!(B_xy = dataFile.get_var("B_xy")))
//		return NC_ERR;

	if (!Ne->get(&electron_dens_mat[0][0][0], Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
		return NC_ERR;
	if (!e_Temp->get(&electron_temp_mat[0][0][0], Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
		return NC_ERR;
	if (!Vel_e->get(&electron_Vele_mat[0][0][0], Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
		return NC_ERR;
	if (!Ni->get(&ion_dens_mat[0][0][0], Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
		return NC_ERR;
	if (!Vel_i->get(&ion_Veli_mat[0][0][0], Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
		return NC_ERR;
	if (!Potential->get(&Potential_mat[0][0][0], Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
		return NC_ERR;
//	if (!B_xy->get(&Bxy_mat[0][0], Pgrid.gridx, Pgrid.gridz))
//		return NC_ERR;
	
	for(unsigned int i=0; i< Pgrid.gridx; i++){
		for(unsigned int k=0; k< Pgrid.gridz; k++){
			Pgrid.Ti[i][k] = 1.0*echarge;
			Pgrid.Te[i][k] = electron_temp_mat[i][k][0];
			Pgrid.na0[i][k] = ion_dens_mat[i][k][0];
			Pgrid.na1[i][k] = electron_dens_mat[i][k][0];
			Pgrid.po[i][k] = Potential_mat[i][k][0];
			Pgrid.ua0[i][k] = ion_Veli_mat[i][k][0];
			Pgrid.ua1[i][k] = electron_Vele_mat[i][k][0];
			Pgrid.bx[i][k] = 0.0;
			Pgrid.by[i][k] = 0.0;
			Pgrid.bz[i][k] = 0.4;
//			Pgrid.bz[i][k] = Bxy_mat[i][k];
		}
	}
	return 0;
}

// Run DTOKS Normally a single time
int DTOKSU_Manager::Run(){
	DM_Debug("\n\nIn DTOKSU_Manager::Run()\n\n");
	if( Config_Status != 0 ){
		std::cerr << "\nDTOKSU Is not configured! Please configure first.";
		config_message();
		return 1;
	}

	std::cout << "\n * RUNNING DTOKS * \n";
	clock_t begin = clock();	// Measure start time

	// Actually running DTOKS
	int RunStatus = Sim->Run();

	clock_t end = clock();		// Measure end time
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;	
	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	std::cout << "\n\n * DTOKS COMPLETED SUCCESSFULLY * \n\n";

	return RunStatus;
}

// Run DTOKS many times with breakup turned on
int DTOKSU_Manager::Breakup(){
	DM_Debug("\n\nIn DTOKSU_Manager::Breakup()\n\n");

	if( Config_Status != 0 ){
		std::cerr << "\nDTOKSU Is not configured! Please configure first.";
		config_message();
		return 1;
	}
	
	
	std::cout << "\n * RUNNING DTOKS * \n";

	threevector Zeroes(0.0,0.0,0.0);

	clock_t begin = clock();

	unsigned int p(1);
	unsigned int i(1);
	
	double seed=std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937 randnumber(seed);
	std::uniform_real_distribution<double> rad(0.0, 1.0); // Uniformly Randomly Distributed Variable between 0.0 and 1.0

	// Vector of Grain Data to record the final state of each branch
	std::vector<GrainData> GDvector;

	// Vectors to retain the time of each simulation branch so that the global time is preserved
	std::vector<double> HMTime;
	std::vector<double> FMTime;
	std::vector<double> CMTime;
	// Loop over the number of split paths.
	for( unsigned int j(0);	j < i; j ++){
		DM_Debug("\nSimulating NEGATIVE Branch "); DM_Debug(i); DM_Debug(" : "); DM_Debug(j);
		DM_Debug("\nStart Pos = "); DM_Debug(Sample->get_position()); 
		DM_Debug("\nVelocity = "); DM_Debug(Sample->get_velocity());
		DM_Debug("\nMass = "); DM_Debug(Sample->get_mass()); 
		DM_Debug("\nTemperature = "); DM_Debug(Sample->get_temperature()); DM_Debug("\n");

		// When breakup occurs and a path forks, track it. If it breaks up, track the subsequent particle
		// Repeat until the end condition is no-longer breakup, i.e while return of DTOKSU object isn't 3.
		while( Sim->Run() == 3 ){ // DTOKSU_Manager has occured...

			// Close data files and open new ones, with new names based off index
			Sim->CloseFiles();
			Sim->OpenFiles("Data/breakup",p);
			
			// Reset breakup so that it's recorded with breakup turned off
			Sample->reset_breakup();

			// Reset the end point data with the same position, no rotation and heading off in negative direction
			// Rotation occurs in random direction perpendicular to magnetic field and velocity as per theory
			double VelocityMag = 2*PI*(Sample->get_radius())*Sample->get_rotationalfreq();
			threevector Unit(2.0*rad(randnumber)-1.0,2.0*rad(randnumber)-1.0,2.0*rad(randnumber)-1.0);
			threevector VelocityUnitVec = (Unit.getunit()^Sim->get_bfielddir()).getunit();
			threevector dvMinus = VelocityMag*VelocityUnitVec;
			Sample->update_motion(Zeroes,dvMinus,-Sample->get_rotationalfreq());

			// Record dust end conditions for later when we simulate other half of fork
			GDvector.push_back(Sample->get_graindata());
			HMTime.push_back(Sim->get_HMTime());
			FMTime.push_back(Sim->get_FMTime());
			CMTime.push_back(Sim->get_CMTime());
	
			// Change the dust velocity, the mass has already been halved in Matter. 
			// Add the velocity twice over as we took it away in one direction
			threevector dvPlus = -2.0*dvMinus;
			Sample->update_motion(Zeroes,dvPlus,0.0);

			DM_Debug("\nSimulating POSITIVE Branch "); DM_Debug(i); DM_Debug(" : "); DM_Debug(j);
			DM_Debug("\nStart Pos = "); DM_Debug(Sample->get_position());
			DM_Debug("\nVelocity = "); DM_Debug(Sample->get_velocity());
			DM_Debug("\nMass = "); DM_Debug(Sample->get_mass()); 
			DM_Debug("\nTemperature = "); DM_Debug(Sample->get_temperature()); DM_Debug("\n");
			// increment counters of number of forks and number of positive forks.
			// p is used for recording index of file
			i = i + 1;
			p = p + 1;
		} // Run the simulation again if breakup occured!
		
		DM_Debug("\n***** START OF : DUST DIDN'T BREAKUP *****\n!");
		DM_Debug(i); DM_Debug(" : "); DM_Debug(j);
		// Close data files and open new ones
		Sim->CloseFiles();
		if( GDvector.size() > 0 ){ // If we have at least one breakup event, i.e stored end point data
			// Re-initialise simulation with new Velocity... 
			// But we've already done this now when we saved the data previously...
			// So we're good to go, just copy over the previous end conditions
			Sample->update_graindata(GDvector[j]);


			// Reset Model Times, the only reason for this is to make the plotting work correctly...
			// Ensure that the time of the models is global and not local to the track
			Sim->ResetModelTime(HMTime[j]-Sim->get_HMTime(),FMTime[j]-Sim->get_FMTime(),CMTime[j]-Sim->get_CMTime());

			// Open files 
			Sim->OpenFiles("Data/breakup",p);
			p = p + 1;
			DM_Debug("\nSimulating POSITIVE Branch "); DM_Debug(i); DM_Debug(" : "); DM_Debug(j);
                        DM_Debug("\nStart Pos = "); DM_Debug(Sample->get_position()); 
                        DM_Debug("\nVelocity = "); DM_Debug(Sample->get_velocity());
                        DM_Debug("\nMass = "); DM_Debug(Sample->get_mass());
                        DM_Debug("\nTemperature = "); DM_Debug(Sample->get_temperature()); DM_Debug("\n");
			Pause();
		}
		DM_Debug("\n***** END OF : DUST DIDN'T BREAKUP *****\n!");
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;	
	std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
	std::cout << "\n\n * DTOKS COMPLETED SUCCESSFULLY * \n\n";
	//	Pgrid.datadump(); // Print the plasma grid data
	return 0;
}
