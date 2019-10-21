/** @file DTOKSU_Manager.cpp
 *  @brief Implementation of class DTOKSU_Manager
 *  
 *  Implement the member functions of the DTOKSU_Manager class
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 */

#include "DTOKSU_Manager.h"

DTOKSU_Manager::DTOKSU_Manager(){
    DM_Debug("In DTOKSU_Manager::DTOKSU_Manager()\n\n");
    Config_Status = -1;
};

DTOKSU_Manager::DTOKSU_Manager(int argc, char* argv[]){
    DM_Debug("In DTOKSU_Manager::DTOKSU_Manager(int argc, "
        << "char* argv[])\n\n");
    std::cout << "\n * CONFIGURING DTOKS * \n";
    Config_Status = -1;

    //!< Call the configure function with command line options to configure as 
    //!< well as construct.
    Config_Status = Configure(argc,argv);
};

DTOKSU_Manager::DTOKSU_Manager(int argc, char* argv[], std::string filename){
    DM_Debug("In DTOKSU_Manager::DTOKSU_Manager(int argc, char* argv[], "
        << "std::string filename = Config/DTOKSU_Config.cfg)\n\n");
    std::cout << "\n * CONFIGURING DTOKS * \n";
    Config_Status = -1;

    //!< Call the configure function with command line options to configure as 
    //!< well as construct.
    Config_Status = Configure(argc,argv,filename);
};


void DTOKSU_Manager::config_message()const{
    DM_Debug("  In DTOKSU_Manager::config_message()\n\n");
    
    if( Config_Status == -3 ){       //!< Configuration processed successfully
        std::cout << "\n * DTOKS SUCCESSFULLY CONFIGURED WITHOUT BREAKUP * \n";
    }else if( Config_Status == -2 ){ //!< Configuration processed successfully
        std::cout << "\n * DTOKS SUCCESSFULLY CONFIGURED WITH BREAKUP * \n";
    }else if( Config_Status == -1 ){ //!< Configuration hasn't been processed
        std::cout << "\n * DTOKS UNCONFIGURED * \n";
    }else if( Config_Status == 0 ){  //!< Configuration processed successfully
        std::cout << "\n * DTOKS SUCCESSFULLY CONFIGURED * \n";
    }else if( Config_Status == 1 ){  //!< Help message was displayed
        std::cout << "\n\n * HELP MESSAGE DISPLAYED * \n\n";
    }else if( Config_Status == 2 ){  //!< Failed to read configuration file
        std::cout << "\n\n * ERROR CODE 2! FAILURE PARSING CONFIGURATION FILE"
            << " * \n\n";
    }else if( Config_Status == 3 ){  //!< Failed to configure plasma data
        std::cout << "\n\n * ERROR CODE 3! FAILURE CONFIGURING PLASMA DATA *"
            << " \n\n";
    }else if( Config_Status == 4 ){  //!< Invalid element option selected
        std::cout << "\n\n * ERROR CODE 4! FAILURE CONFIGURING SAMPLE * \n\n";
    }else if( Config_Status == 5 ){  //!< Invalid element option selected
        std::cout << "\n\n * ERROR CODE 5! FAILURE CONFIGURING BREAKUP * \n\n";
    }else if( Config_Status == 6 ){  //!< Failed to configure boundary data
        std::cout << "\n\n * ERROR CODE 6! FAILURE CONFIGURING BOUNDARY DATA *"
            << " \n\n";
    }else if( Config_Status == 7 ){  //!< Failed to configure current terms data
        std::cout << "\n\n * ERROR CODE 7! FAILURE CONFIGURING CURRENT TERMS *"
            << " \n\n";
    }else{
        std::cout << "\n\n * UNKNOWN CONFIGURATION STATUS * \n\n";
    }
}

bool DTOKSU_Manager::check_pdata_range()const{
    DM_Debug("  In DTOKSU_Manager::check_pdata_range()\n\n");
    bool ReturnVal = true;
    if( Pdata.IonDensity > Overflows::Density 
        || Pdata.IonDensity < Underflows::Density ){
        std::cout << "\n* Error, Configured Ion Density outside valid range! *";
        std::cout << "\n* Density must be in range: " << Underflows::Density 
            << "<= Density <=" << Overflows::Density;  
        ReturnVal = false;
    }
    if( Pdata.ElectronDensity > Overflows::Density 
        || Pdata.ElectronDensity < Underflows::Density ){
        std::cout << "\n* Error, Configured Electron Density outside valid "
            << "range! *";
        std::cout << "\n* Density must be in range: " << Underflows::Density 
            << "<= Density <=" << Overflows::Density;  
        ReturnVal = false;
    }
    if( Pdata.NeutralDensity > Overflows::Density 
        || Pdata.NeutralDensity < Underflows::Density ){
        std::cout << "\n* Error, Configured Neutral Density outside valid ";
        std::cout << "range! *\n* Density must be in range: " 
            << Underflows::Density << "<= Density <=" << Overflows::Density;  
        ReturnVal = false;
    }
    if( Pdata.IonTemp > Overflows::Temperature 
        || Pdata.IonTemp < Underflows::Temperature ){
        std::cout << "\n* Error, Configured Ion Temperature outside valid ";
        std::cout << "range! *\n* Temperature must be in range: " 
            << Underflows::Temperature << "<= Temperature <=" 
            << Overflows::Temperature;  
        ReturnVal = false;
    }
    if( Pdata.ElectronTemp > Overflows::Temperature 
        || Pdata.ElectronTemp < Underflows::Temperature ){
        std::cout << "\n* Error, Configured Electron Temperature outside ";
        std::cout << "valid range! *\n* Temperature must be in range: " 
            << Underflows::Temperature << "<= Temperature <=" 
            << Overflows::Temperature;  
        ReturnVal = false;
    }
    if( Pdata.NeutralTemp > Overflows::Temperature 
        || Pdata.NeutralTemp < Underflows::Temperature ){
        std::cout << "\n* Error, Configured Neutral Temperature outside valid ";
        std::cout << "range! *\n* Temperature must be in range: " 
            << Underflows::Temperature << "<= Temperature <=" 
            << Overflows::Temperature;  
        ReturnVal = false;
    }
    if( Pdata.AmbientTemp > Overflows::Temperature 
        || Pdata.AmbientTemp < Underflows::Temperature ){
        std::cout << "\n* Error, Configured Ambient Temperature outside valid";
        std::cout << " range! *\n* Temperature must be in range: " 
            << Underflows::Temperature << "<= Temperature <=" 
            << Overflows::Temperature;  
        ReturnVal = false;
    }
    if( Pdata.PlasmaVel.mag3() > Overflows::PlasmaVel 
        || Pdata.PlasmaVel.mag3() < Underflows::PlasmaVel ){
        std::cout << "\n* Error, Configured Plasma velocity outside valid ";
        std::cout << "range! *\n* Plasma Velocity Magnitude must be in range: " 
            << Underflows::PlasmaVel << "<= PlasmaVel <=" 
            << Overflows::PlasmaVel;  
        ReturnVal = false;
    }
    if( Pdata.ElectricField.mag3() > Overflows::Field 
        || Pdata.ElectricField.mag3() < Underflows::Field ){
        std::cout << "\n* Error, Configured Electric Field outside valid ";
        std::cout << "range! *\n* Field Magnitude must be in range: " 
            << Underflows::Field << "<= Field <=" << Overflows::Field;    
        ReturnVal = false;
    }
    if( Pdata.MagneticField.mag3() > Overflows::Field 
        || Pdata.MagneticField.mag3() < Underflows::Field ){
        std::cout << "\n* Error, Configured Magnetic Field Field outside ";
        std::cout << "valid range! *\n* Field Magnitude must be in range: " 
            << Underflows::Field << "<= Field <=" << Overflows::Field;    
        ReturnVal = false;
    }
    if( Pdata.Gravity.mag3() > Overflows::Field 
        || Pdata.Gravity.mag3() < Underflows::Field ){
        std::cout << "\n* Error, Configured Gravity Field outside valid ";
        std::cout << "range! *\n* Field Magnitude must be in range: " 
            << Underflows::Field << "<= Field <=" << Overflows::Field;    
        ReturnVal = false;
    }
    return ReturnVal;
}

void DTOKSU_Manager::show_usage(std::string name)const{
    DM_Debug("  In DTOKSU_Manager::show_usage(std::string name)\n\n");
    std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
    << "\n\nOptions:\n"
    << "\t-h, --help\t\tShow this help message\n\n"
    << "\t-c, --config CONFIG\t\tSpecify the configuration file path\n\n"
    << "\t-i, --imax IMAX\t\t\t(int), Specify the number of grains\n\n"
    << "\t-se,--seed SEED\t\t\tSpecify the seed for simulations\n\n"
    << "\t-vd,--veldist VELDIST\t\tSpecify velocity dist:\n"
    << "\t\t\t\t\tuniform, lognormal, none\n\n"
    << "\t-sd,--sizedist SIZEDIST\t\tSpecify size dist:\n"
    << "\t\t\t\t\tuniform, exponential, lognormal, none\n\n"
    << "\t-vpo,--velparone VELPARONE\tSpecify the 1st par of velocity dist\n\n"
    << "\t-vpt,--velpartwo VELPARTWO\tSpecify the 2nd par of velocity dist\n\n"
    << "\t-sdo,--sizeparone SIZEPARONE\tSpecify the 1st par of size dist\n\n"
    << "\t-sdt,--sizepartwo SIZEPARTWO\tSpecify the 2nd par of size dist\n\n"
    << "\t-t, --temperature TEMPERATURE\tdouble the initial temperature\n\n"
    << "\t-m, --material MATERIAL\t\tchar variable giving the dust grain "
    << "element, possible values 'w', 'g', 'b', 'd', 'l', 'm' and 'f'\n"
    << "\t\t\t\t\t(W): Tungsten, (G): Graphite, (B): Beryllium, (D): Deuterium "
    << "(L): Lithium, (M): Molybdenum or (F): Iron\n\n"
    << "\t-s, --size SIZE\t\t\tdouble the radius of the grain\n\n"
    << "\t-vr,--rvel RVEL\t\t\tfloat radial velocity\n\n"
    << "\t-vt,--thetavel THETAVEL\t\tfloat angular velocity\n\n"
    << "\t-vz,--zvel ZVEL\t\t\tfloat lognitudinal velocity\n\n"
    << "\t-rr,--rpos RPOS\t\t\tfloat radial position\n\n"
    << "\t-rt,--thetapos THETAPOS\t\tfloat angular position\n\n"
    << "\t-rz,--zpos ZPOS\t\t\tfloat longitudinal position\n\n"
    << "\t-op,--output OUTPUT\t\tstring the filename prefix to write to\n\n"
    << "\t-om,--metadata METADATA\t\tstring the MetaData filename to write\n\n";
}

template<typename T> int DTOKSU_Manager::input_function(int &argc, char* argv[],
int &i, std::stringstream &ss0, T &Temp)const{
    DM_Debug("  In DTOKSU_Manager::input_function(int &argc, char* argv[],"
        << " int &i, std::stringstream &ss0, T &Temp)const\n\n");
    if (i + 1 < argc) { //!< Make sure we aren't at the end of argv!
        i+=1; //!< Increment 'i' so we don't get the next argv[i].
        ss0 << argv[i];
        ss0 >> Temp;
        ss0.clear(); ss0.str("");
        return 0;
    }else{ //!< Uh-oh, there was no argument to the destination option.
        std::cerr << "\noption requires argument." << std::endl;
        return 1;
    }
}

//!< Configure DTOKS ready to run with a particular plasma and matter sample
//!< which are established from configuration filename
int DTOKSU_Manager::Configure(std::string Config_Filename){
    DM_Debug("  In DTOKSU_Manager::Configure(std::string Config_Filename)"
        << "\n\n");
    char* argv[0] = {};
    return Configure(0,argv,Config_Filename);
}

//!< Configure DTOKS ready to run with a particular plasma and matter sample
//!< which are established from command line options and DTOKSU_Config.cfg 
//!< configuration file
int DTOKSU_Manager::Configure(int argc, char* argv[], 
std::string Config_Filename){
    DM_Debug("  In DTOKSU_Manager::Configure(int argc, char* argv[], "
        << "std::string Config_Filename=Config_Files/DTOKSU_Config.cfg)\n\n");

    //!< Default data file prefix and plasma data directory
    std::string MetaDataFilename = "Data/DTOKSU.txt";
    std::string DataFilePrefix = "Data/DTOKSU";
    std::string PlasmaData_dir = "PlasmaData/";
    std::string WallData_dir = "PlasmaData/";
    std::string CoreData_dir = "PlasmaData/";


    // ------------------- PARSE COMMAND LINE INPUT ------------------- //
    // Check user input for specifying the configuration file.
    // We want other command line options to over-ride the configuration file.
    std::vector <std::string> sources;
    std::stringstream ss0;
    for (int i = 1; i < argc; ++i){ // Read command line input
        std::string arg = argv[i];
//	std::cout << "\narg[" << i << "] = " << arg;
        if     ( arg == "--help"        || arg == "-h" ){   
            Config_Status = 1; show_usage( argv[0] ); return Config_Status;         
        }else if( arg == "--config"      || arg == "-c" ){
            input_function(argc,argv,i,ss0,Config_Filename);
        }else{
            sources.push_back(argv[i]);
        }
    }


    // ------------------- PARALLELISATION VARIABLES ------------------- //
    unsigned long long imax = 1;
    unsigned long long Seed = 0;
    double SizeParOne = 0;
    double SizeParTwo = 0;
    double VelParOne = 0;
    double VelParTwo = 0;

    std::string SizeDistributionType = "none";
    std::string VelDistributionType = "none";
    int ValueTwo(0);


    // ------------------- DUST VARIABLE DEFAULTS ------------------- //
    char Element='W';
    char IonSpecies='h';
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
    std::vector<CurrentTerm*> CurrentTerms;
    std::vector<ForceTerm*> ForceTerms;
    std::vector<HeatTerm*> HeatTerms;
    std::array<char,CM> ConstModels;
    std::array<float,DTOKSU::MN> AccuracyLevels;

    config4cpp::StringVector CfgStringVec;


    // ------------------- PROCESS CONFIGURATION FILE ------------------- //
//  threevector PlasmaVelocity(501.33, 7268.5, 914.947); 
//  threevector Efield(-13.673, 0, -27.925);
//  threevector Bfield(0.0226868, 0.328923, 0.0414043);
    config4cpp::Configuration *  cfg = config4cpp::Configuration::create();
    std::cout << "\n* Reading configuration file: " << Config_Filename << " *";
    try {
        cfg->parse(Config_Filename.c_str());
        MetaDataFilename = cfg->lookupString("", "Filename");
        DataFilePrefix   = cfg->lookupString("", "DataFilePrefix");
        ContinuousPlasma = cfg->lookupBoolean("plasma", "ContinuousPlasma");
        IonSpecies       = cfg->lookupString("plasma", "Plasma")[0];
        Pdata.Z          = cfg->lookupFloat("plasma", "MeanIonization");
        if( !ContinuousPlasma ){
            PlasmaData_dir= cfg->lookupString("plasma","plasmagrid.Plasmadir");
            WallData_dir  = cfg->lookupString("plasma","plasmagrid.Walldir");
            CoreData_dir  = cfg->lookupString("plasma","plasmagrid.Coredir");
            Pgrid.device  = cfg->lookupString("plasma","plasmagrid.Machine")[0];
            Pgrid.dlx     = cfg->lookupFloat ("plasma","plasmagrid.xSpacing");
            Pgrid.dlz     = cfg->lookupFloat ("plasma","plasmagrid.zSpacing");
        }
        Pdata.IonDensity      
            = cfg->lookupFloat("plasma","plasmadata.IonDensity");
        Pdata.ElectronDensity 
            = cfg->lookupFloat("plasma","plasmadata.ElectronDensity");
        Pdata.NeutralDensity  
            = cfg->lookupFloat("plasma","plasmadata.NeutralDensity");
        Pdata.IonTemp         
            = cfg->lookupFloat("plasma","plasmadata.IonTemp");
        Pdata.NeutralTemp     
            = cfg->lookupFloat("plasma","plasmadata.NeutralTemp");
        Pdata.ElectronTemp    
            = cfg->lookupFloat("plasma","plasmadata.ElectronTemp");
        Pdata.AmbientTemp     
            = cfg->lookupFloat("plasma","plasmadata.AmbientTemp");

        cfg->lookupList("plasma","plasmadata.PlasmaVelocity",CfgStringVec);
        if( CfgStringVec.length() != 3 ){
            std::cout << "\nError, reading configuration file!\n"
                << "plasmadata.PlasmaVelocity length != 3";
        }else{
            Pdata.PlasmaVel = threevector(std::stof(CfgStringVec[0]),
                std::stof(CfgStringVec[1]),std::stof(CfgStringVec[2]));
        }
        cfg->lookupList("plasma","plasmadata.Gravity",CfgStringVec);
        if( CfgStringVec.length() != 3 ){
            std::cout << "\nError, reading configuration file!\n"
                << "plasmadata.Gravity length != 3";
        }else{
            Pdata.Gravity = threevector(std::stof(CfgStringVec[0]),
                std::stof(CfgStringVec[1]),std::stof(CfgStringVec[2]));
        }
        cfg->lookupList("plasma","plasmadata.Efield",CfgStringVec);
        if( CfgStringVec.length() != 3 ){
            std::cout << "\nError, reading configuration file!\n"
                << "plasmadata.Efield length != 3";
        }else{          
            Pdata.ElectricField = threevector(std::stof(CfgStringVec[0]),
                std::stof(CfgStringVec[1]),std::stof(CfgStringVec[2]));
        }
        cfg->lookupList("plasma","plasmadata.Bfield",CfgStringVec);
        if( CfgStringVec.length() != 3 ){
            std::cout << "\nError, reading configuration file!\n"
                << "plasmadata.Bfield length != 3";
        }else{          
            Pdata.MagneticField = threevector(std::stof(CfgStringVec[0]),
                std::stof(CfgStringVec[1]),std::stof(CfgStringVec[2]));
        }
        Element  = cfg->lookupString("dust", "Element")[0];
        size     = cfg->lookupFloat("dust", "size");
        Temp     = cfg->lookupFloat("dust", "Temp");
        rpos     = cfg->lookupFloat("dust", "dynamics.rpos");
        thetapos = cfg->lookupFloat("dust", "dynamics.thetapos");
        zpos     = cfg->lookupFloat("dust", "dynamics.zpos");
        rvel     = cfg->lookupFloat("dust", "dynamics.rvel");
        thetavel = cfg->lookupFloat("dust", "dynamics.thetavel");
        zvel     = cfg->lookupFloat("dust", "dynamics.zvel");
        InitRotationalFreq 
            = cfg->lookupFloat("dust", "dynamics.InitRotationalFreq");
        ConstModels = 
            {
                cfg->lookupString("variablemodels", "EmissivityModel")[0], 
                cfg->lookupString("variablemodels", "ExpansionModel")[0], 
                cfg->lookupString("variablemodels", "HeatCapacityModel")[0], 
                cfg->lookupString("variablemodels", "BoilingModel")[0],
                cfg->lookupString("variablemodels", "BreakupModel")[0]
            };
        HeatModels = 
            {
                cfg->lookupBoolean("heatingmodels","RadiativeCooling"),
                cfg->lookupBoolean("heatingmodels","EvaporativeCooling"),
                cfg->lookupBoolean("heatingmodels","NewtonCooling"),
                cfg->lookupBoolean("heatingmodels","NeutralHeatFlux"),
                cfg->lookupBoolean("heatingmodels","SOMLIonHeatFlux"),
                cfg->lookupBoolean("heatingmodels","SMOMLIonHeatFlux"),
                cfg->lookupBoolean("heatingmodels","DTOKSIonHeatFlux"),
                cfg->lookupBoolean("heatingmodels","DUSTTIonHeatFlux"),
                cfg->lookupBoolean("heatingmodels","OMLElectronHeatFlux"),
                cfg->lookupBoolean("heatingmodels","PHLElectronHeatFlux"),
                cfg->lookupBoolean("heatingmodels","DTOKSElectronHeatFlux"),
                cfg->lookupBoolean("heatingmodels","SOMLNeutralRecombination"),
                cfg->lookupBoolean("heatingmodels","SMOMLNeutralRecombination"),
                cfg->lookupBoolean("heatingmodels","DTOKSNeutralRecombination"),
                cfg->lookupBoolean("heatingmodels","SEE"),
                cfg->lookupBoolean("heatingmodels","DTOKSSEE"),
                cfg->lookupBoolean("heatingmodels","TEE"),
                cfg->lookupBoolean("heatingmodels","DTOKSTEE")
            };
        ForceModels =
            {
                cfg->lookupBoolean("forcemodels","Gravity"),
                cfg->lookupBoolean("forcemodels","Lorentz"), 
                cfg->lookupBoolean("forcemodels","SOMLIonDrag"), 
                cfg->lookupBoolean("forcemodels","SMOMLIonDrag"), 
                cfg->lookupBoolean("forcemodels","DTOKSIonDrag"),
                cfg->lookupBoolean("forcemodels","DUSTTIonDrag"), 
                cfg->lookupBoolean("forcemodels","HybridDrag"), 
                cfg->lookupBoolean("forcemodels","LloydDrag"), 
                cfg->lookupBoolean("forcemodels","NeutralDrag"), 
                cfg->lookupBoolean("forcemodels","RocketForce")
            };
        ChargeModels =
            {
                cfg->lookupBoolean("chargemodels","OMLe"), 
                cfg->lookupBoolean("chargemodels","PHLe"), 
                cfg->lookupBoolean("chargemodels","DTOKSe"),
                cfg->lookupBoolean("chargemodels","THSe"), 
                cfg->lookupBoolean("chargemodels","OMLi"),  
                cfg->lookupBoolean("chargemodels","MOMLi"), 
                cfg->lookupBoolean("chargemodels","SOMLi"), 
                cfg->lookupBoolean("chargemodels","SMOMLi"),
                cfg->lookupBoolean("chargemodels","THSi"), 
                cfg->lookupBoolean("chargemodels","DTOKSi"), 
                cfg->lookupBoolean("chargemodels","TEE"),
                cfg->lookupBoolean("chargemodels","TEESchottky"), 
                cfg->lookupBoolean("chargemodels","SEE"), 
                cfg->lookupBoolean("chargemodels","CW"),
                cfg->lookupBoolean("chargemodels","MOMLWEM")
            };
        
        AccuracyLevels = 
            {
                cfg->lookupFloat("accuracylevels", "charge"), 
                cfg->lookupFloat("accuracylevels", "heat"),
                cfg->lookupFloat("accuracylevels", "force")
                
            };
    } catch(const config4cpp::ConfigurationException & ex) {
        std::cerr << ex.c_str() << std::endl;
        cfg->destroy();
        Config_Status = 2;
        return Config_Status;
    }

    //!< Tell user if the simulation is self-consistent.
    /*if( HeatModels[4] && HeatModels[5] && ForceModels[3] && ChargeModels[3] ){
        std::cout << "\n* SOML Self-Consistent Simulation! *";
    }else if( HeatModels[6] && HeatModels[7] && ForceModels[4] 
        && ChargeModels[4] ){
        std::cout << "\n* SMOML Self-Consistent Simulation! *";
    }else if( HeatModels[8] && HeatModels[9] && HeatModels[10] 
        && ChargeModels[8] ){
        std::cout << "\n* PHL Self-Consistent Simulation! *";
    }else if( HeatModels[12] && HeatModels[13] && HeatModels[14] 
        && HeatModels[15] && HeatModels[16] && ForceModels[5] 
        && (ChargeModels[7] || ChargeModels[8]) ){
        std::cout << "\n* DTOKS Simulation! Non-self Consistent fluxes! *";
    }else{
        std::cout << "\n* Non-self Consistent fluxes! *";
    }*/

    if( ChargeModels[0] )      CurrentTerms.push_back(new Term::OMLe());
    else if( ChargeModels[1] ) CurrentTerms.push_back(new Term::PHLe());
    else if( ChargeModels[2] ) CurrentTerms.push_back(new Term::THSe());
    else if( ChargeModels[3] ) CurrentTerms.push_back(new Term::DTOKSe());
    else{
        std::cout << "\n* No electron Current Model Chosen! *";
        Config_Status = 7;
        return Config_Status;
    }
    if( ChargeModels[4] )      CurrentTerms.push_back(new Term::OMLi());
    else if( ChargeModels[5] ) CurrentTerms.push_back(new Term::MOMLi());
    else if( ChargeModels[6] ) CurrentTerms.push_back(new Term::SOMLi());
    else if( ChargeModels[7] ) CurrentTerms.push_back(new Term::SMOMLi());
    else if( ChargeModels[8] ) CurrentTerms.push_back(new Term::THSi());
    else if( ChargeModels[9] ) CurrentTerms.push_back(new Term::DTOKSi());
    else{
        std::cout << "\n* No ion Current Model Chosen! *";
        Config_Status = 7;
        return Config_Status;
    }
    if( ChargeModels[10] )      CurrentTerms.push_back(new Term::TEEcharge());
    else if( ChargeModels[11] ) CurrentTerms.push_back(new Term::TEESchottky());
    if( ChargeModels[12] )      CurrentTerms.push_back(new Term::SEEcharge());
    if( ChargeModels[13] && !ChargeModels[14] ){
        CurrentTerms.clear();
        CurrentTerms.push_back(new Term::CW());
    }else if( !ChargeModels[13] && ChargeModels[14] ){
        CurrentTerms.clear();
        CurrentTerms.push_back(new Term::MOMLWEM());
    }else if( ChargeModels[13] && ChargeModels[14] ){
        std::cout << "\n* CW and MOMLWEM Not compatible models! *";
        Config_Status = 7;
        return Config_Status;
    }

    if( ForceModels[0] ) ForceTerms.push_back(new Term::Gravity());
    if( ForceModels[1] ) ForceTerms.push_back(new Term::LorentzForce());
    if( ForceModels[2] ) ForceTerms.push_back(new Term::SOMLIonDrag());
    else if( ForceModels[3] ) ForceTerms.push_back(new Term::SMOMLIonDrag());
    else if( ForceModels[4] ) ForceTerms.push_back(new Term::DTOKSIonDrag());
    else if( ForceModels[5] ) ForceTerms.push_back(new Term::DUSTTIonDrag());
    else if( ForceModels[6] ) ForceTerms.push_back(new Term::HybridIonDrag());
    if( ForceModels[7] ) ForceTerms.push_back(new Term::LloydIonDrag());
    if( ForceModels[8] ) ForceTerms.push_back(new Term::NeutralDrag());
    if( ForceModels[9] ) ForceTerms.push_back(new Term::RocketForce());
    
    if( HeatModels[0] )  HeatTerms.push_back(new Term::EmissivityModel());
    if( HeatModels[1] )  HeatTerms.push_back(new Term::EvaporationModel());
    if( HeatModels[2] )  HeatTerms.push_back(new Term::NewtonCooling());
    if( HeatModels[3] )  HeatTerms.push_back(new Term::NeutralHeatFlux());
    if( HeatModels[4] ) HeatTerms.push_back(new Term::OMLElectronHeatFlux());
    else if( HeatModels[5] ) HeatTerms.push_back(new Term::PHLElectronHeatFlux());
    else if( HeatModels[6] ) HeatTerms.push_back(new Term::DTOKSElectronHeatFlux());
    if( HeatModels[7] )  HeatTerms.push_back(new Term::SOMLIonHeatFlux());
    else if( HeatModels[8] )  HeatTerms.push_back(new Term::SMOMLIonHeatFlux());
    else if( HeatModels[9] ) HeatTerms.push_back(new Term::DTOKSIonHeatFlux());
    else if( HeatModels[10] ) HeatTerms.push_back(new Term::DUSTTIonHeatFlux());
    if( HeatModels[11] )  HeatTerms.push_back(new Term::SOMLNeutralRecombination());
    else if( HeatModels[12] )  HeatTerms.push_back(new Term::SMOMLNeutralRecombination());
    else if( HeatModels[13] ) HeatTerms.push_back(new Term::DTOKSNeutralRecombination());
    if( HeatModels[14] )  HeatTerms.push_back(new Term::SEE());
    else if( HeatModels[15] ) HeatTerms.push_back(new Term::DTOKSSEE());
    if( HeatModels[16] )  HeatTerms.push_back(new Term::TEE());
    else if( HeatModels[17] ) HeatTerms.push_back(new Term::DTOKSTEE());
    
    cfg->destroy();
    if( !check_pdata_range() ){
        Config_Status = 3;
        return Config_Status;
    }

    std::cout << "\n* Configuration file read successfully! *";
    Pause();

    //!< For Hydrogen plasma, we have hydrogen ions of mass Mp
    if( IonSpecies == 'h' ){ 
        Pgrid.mi = Mp;
        Pdata.mi = Mp;
        Pdata.A = 1.0;
        if( Pdata.Z > 1.0 ){
            std::cout << "\nPdata.Z < 1.0 for Hydrogen! Assuming Pdata.Z = 1.0";
        }
        Pdata.Z = 1.0;
    }else{
        std::cout << "Invalid IonSpecies, IonSpecies = " << IonSpecies << "!\n";
        return 1;
    }


    // ------------------- CONFIGURE PLASMAGRID ------------------- //
    Boundary_Data WallBound, CoreBound;
    if( !ContinuousPlasma ){
        std::cout << "\n\n* Full Machine Simulation! *\n\n"
            << "* Creating PlasmaGrid_Data Structure *";
        std::cout << "\n\t* Plasma:\t" << IonSpecies << "\n\t* Machine:\t" 
            << Pgrid.device;
        std::cout << "\n\t* xSpacing:\t" << Pgrid.dlx << "\n\t* zSpacing:\t" 
            << Pgrid.dlz;
        //!< Failed to configure plasma data
        if( configure_plasmagrid(PlasmaData_dir) != 0 ){ 
            std::cerr << "\nFailed to configure plasma data!";
            Config_Status = 3;
            return Config_Status;
        }
        std::cout << "\n* PlasmaGrid_Data Structure created successfully! *";
        if( WallData_dir != ""){
            //!< Failed to configure plasma data
            if( configure_boundary(WallData_dir,"WallData.txt",WallBound) != 0 )
            { 
                std::cerr << "\nFailed to configure plasma data!";
                Config_Status = 6;
                return Config_Status;
            }
            std::cout << "\n* Wall_Data Structure created successfully! *\n";
        }else{
            std::cout << "\n* No Wall Data Used! *\n";
        }
        if( CoreData_dir != "" ){
            if(Pgrid.device == 'p'){ // In case of magnum-PSI
                std::cout << "\n* Core boundary not valid for Magnum-PSI!"
                    << " Ignoring filepath:" << CoreData_dir << " *\n";
            }else{
                std::cout << "\n* Trajectories ending at Core Boundary! *\n";
                if( configure_boundary(CoreData_dir,"CoreData.txt",CoreBound) 
                    != 0 ){ //!< Failed to configure plasma data
                    std::cerr << "\nFailed to configure plasma data!";
                    Config_Status = 6;
                    return Config_Status;
                }
                std::cout 
                    << "\n* Core_Boundary Structure created successfully! *\n";
            }
        }else if(Pgrid.device != 'p'){ //!< In case of magnum-PSI
            std::cout << "\n* Trajectories pass through Core Boundary with "
                << "interpolated plasma data!";
        }
        
    }else{
        std::cout << "\n\n* ContinuousPlasma! *";
        std::cout << "\n* PlasmaData Structure created successfully! *\n";
    }
    Pause();


    // ------------------- PARSE COMMAND LINE INPUT ------------------- //
    std::cout << "\n* Processing command line input *";

    //!< Check user input for specifying the configuration file.
    //!< We want other command line options to over-ride the configuration file.
    for (int i = 1; i < argc; ++i){ //!< Read command line input
        std::string arg = argv[i];
        if( arg == "--imax"
            || arg == "-i"   ) input_function(argc,argv,i,ss0,imax);
        else if( arg == "--seed" 
            || arg == "-se" )  input_function(argc,argv,i,ss0,Seed);
        else if( arg == "--veldist" 
            || arg == "-vd" )  input_function(argc,argv,i,ss0,VelDistributionType);
        else if( arg == "--sizedist" 
            || arg == "-sd" )  input_function(argc,argv,i,ss0,SizeDistributionType);
        else if( arg == "--velparone" 
            || arg == "-vpo" )  input_function(argc,argv,i,ss0,VelParOne);
        else if( arg == "--velpartwo" 
            || arg == "-vpt" )  input_function(argc,argv,i,ss0,VelParTwo);
        else if( arg == "--sizeparone" 
            || arg == "-spo" )  input_function(argc,argv,i,ss0,SizeParOne);
        else if( arg == "--sizepartwo" 
            || arg == "-spt" )  input_function(argc,argv,i,ss0,SizeParTwo);
        else if( arg == "--temperature" 
            || arg == "-t"   ) input_function(argc,argv,i,ss0,Temp);
        else if( arg == "--material"    
            || arg == "-m"   ) input_function(argc,argv,i,ss0,Element);
        else if( arg == "--size"        
            || arg == "-s"   ) input_function(argc,argv,i,ss0,size);
        else if( arg == "--rvel"        
            || arg == "-vr"  ) input_function(argc,argv,i,ss0,rvel);
        else if( arg == "--thetavel"    
            || arg == "-vt"  ) input_function(argc,argv,i,ss0,thetavel);
        else if( arg == "--zvel"        
            || arg == "-vz"  ) input_function(argc,argv,i,ss0,zvel);
        else if( arg == "--rpos"        
            || arg == "-rr"  ) input_function(argc,argv,i,ss0,rpos);
        else if( arg == "--thetapos"    
            || arg == "-rt"  ) input_function(argc,argv,i,ss0,thetapos);
        else if( arg == "--zpos"        
            || arg == "-rz"  ) input_function(argc,argv,i,ss0,zpos);
        else if( arg == "--output"      
            || arg == "-op"  ) input_function(argc,argv,i,ss0,DataFilePrefix);
        else if( arg == "--MetaData"    
            || arg == "-om"  ) input_function(argc,argv,i,ss0,MetaDataFilename);
        else{
            sources.push_back(argv[i]);
        }
    }
    std::cout << "\n* Command line input processed successfully! *\n\n";
    Pause();


    // ------------------- BEGIN PARALLELISATION ------------------- //


    int run_status = 0;
    std::vector<std::mt19937> randnumbers;
    for(int p = 0; p < omp_get_max_threads(); p ++){
        randnumbers.push_back(std::mt19937(Seed+p));
    }

    #pragma omp parallel for private(size)
    for( int i = 0; i < imax; i ++ ){
        //std::cout << "\n" << omp_get_thread_num() << "/" << omp_get_num_threads();

        int Config_Status_Local(0);
        if( SizeDistributionType == "lognormal" ){
            std::lognormal_distribution<double> size_distribution(log(SizeParOne),SizeParTwo);
            size = size_distribution(randnumbers[omp_get_thread_num()]);
        }else if( SizeDistributionType == "uniform" ){
            std::uniform_real_distribution<double> size_distribution(SizeParOne,SizeParTwo);
            size = size_distribution(randnumbers[omp_get_thread_num()]);
        }else if( SizeDistributionType == "exponential" ){
            std::exponential_distribution<double> size_distribution(SizeParTwo);
            size = SizeParOne*size_distribution(randnumbers[omp_get_thread_num()]);
        }else if( SizeDistributionType == "none" ){
            
        }else{
            Config_Status_Local = 4;
        }


        double VelocityMag(0.0);   
        if( VelDistributionType == "lognormal" ){
            std::lognormal_distribution<double> velocity_distribution(log(VelParOne),VelParTwo);
            VelocityMag = velocity_distribution(randnumbers[omp_get_thread_num()]);
        }else if( VelDistributionType == "uniform" ){
            std::uniform_real_distribution<double> velocity_distribution(VelParOne,VelParTwo);
            VelocityMag = velocity_distribution(randnumbers[omp_get_thread_num()]);
        }else if( VelDistributionType == "none" ){

        }else{
            Config_Status_Local = 4;
        }
    //    std::uniform_real_distribution<double> Velocity(0,20);
    //    std::uniform_real_distribution<double> Positions(0,0.1); // UPPER OUTER DIVERTOR
    //    std::uniform_real_distribution<double> Positions(0,1.0); // UPPER INNER DIVERTOR
    //    std::uniform_real_distribution<double> Positions(0,0.25); // LOWER OUTER DIVERTOR
        std::uniform_real_distribution<double> Positions(0,0.03); // LOWER INNER DIVERTOR
    
        double Angle = custom_distribution(randnumbers[omp_get_thread_num()]);
        double DistanceAlongDivertor = Positions(randnumbers[omp_get_thread_num()]);
    
    //    rpos = 1.70705+DistanceAlongDivertor*0.242535625; // UPPER OUTER DIVERTOR
    //    zpos = 1.15-DistanceAlongDivertor*0.9701425001;   // UPPER OUTER DIVERTOR
    //    rpos = 1.4440497335701599;                              // UPPER INNER DIVERTOR
    //    zpos = 0.873889876-DistanceAlongDivertor*0.04618117229; // UPPER INNER DIVERTOR
    //    rpos = 1.76377+DistanceAlongDivertor*0.2594386678;  // LOWER OUTER DIVERTOR
    //    zpos = -1.15+DistanceAlongDivertor*0.9657595858; // LOWER OUTER DIVERTOR
        rpos = 1.423;  // LOWER INNER DIVERTOR
        zpos = -0.7282415630550609-DistanceAlongDivertor; // LOWER INNER DIVERTOR
    
    //    rvel = -VelocityMag*std::cos(Angle-atan(1.0/4.0)); // UPPER OUTER DIVERTOR
    //    zvel = VelocityMag*std::sin(Angle-atan(1.0/4.0));  // UPPER OUTER DIVERTOR
    //    rvel = VelocityMag*std::cos(Angle); // UPPER INNER DIVERTOR
    //    zvel = VelocityMag*std::sin(Angle); // UPPER INNER DIVERTOR
    //    rvel = -VelocityMag*std::cos(Angle+atan(1.0/4.0)); // LOWER OUTER DIVERTOR
    //    zvel = VelocityMag*std::sin(Angle+atan(1.0/4.0)); // LOWER OUTER DIVERTOR
        rvel = VelocityMag*std::cos(Angle); // LOWER INNER DIVERTOR
        zvel = VelocityMag*std::sin(Angle); // LOWER INNER DIVERTOR
//        #pragma omp critical
//        {
//        std::cout << "\n" << size << "\t" << VelocityMag << "\t" << rpos << "\t" << zpos 
//            << "\t" << rvel << "\t" << zvel << "\t" << Angle;
//        }
        
        // Change filenames
        std::string FinalMetaDataFilename = MetaDataFilename + "_" + std::to_string(i) + ".txt";
//        std::string FinalDataFilePrefix = DataFilePrefix + "_" + std::to_string(VelocityMag) + "ms_" + std::to_string(size) + "um";
        std::string FinalDataFilePrefix = DataFilePrefix;
    
        // ------------------- INITIALISE META_DATA FILE ------------------- //
    
        #pragma omp critical
        {
//        std::cout << "* Creating MetaDataFile: " << FinalMetaDataFilename 
//            << "_" << i << ".txt *\n\n";
        std::cout << "* Creating MetaDataFile: " << FinalMetaDataFilename << " *\n\n";
        }
        //!< Data file for containing the run information
        std::ofstream MetaDataFile; 
        time_t now = time(0); //!< Get the time of simulation
        char * dt = ctime(&now);
        MetaDataFile.open(FinalMetaDataFilename);
        MetaDataFile << "## Run Data File ##\n";
        MetaDataFile << "#Date:\t" << dt;
    
    
        // ------------------- INITIALISE DUST ------------------- //
        threevector xinit(rpos,thetapos,zpos);
        threevector vinit(rvel,thetavel,zvel);
        #pragma omp critical
        {
        std::cout << "* Creating Matter object *\n\t* Element:\t" << Element;
        }
        Matter *SampleLocal;
        if  (Element == 'W')     
            SampleLocal = new Tungsten(size,Temp,ConstModels,xinit,vinit);
        else if (Element == 'B') 
            SampleLocal = new Beryllium(size,Temp,ConstModels,xinit,vinit);
        else if (Element == 'F') 
            SampleLocal = new Iron(size,Temp,ConstModels,xinit,vinit);
        else if (Element == 'G') 
            SampleLocal = new Graphite(size,Temp,ConstModels,xinit,vinit);
        else if (Element == 'D') 
            SampleLocal = new Deuterium(size,Temp,ConstModels,xinit,vinit);
        else if (Element == 'M') 
            SampleLocal = new Molybdenum(size,Temp,ConstModels,xinit,vinit);
        else if (Element == 'L') 
            SampleLocal = new Lithium(size,Temp,ConstModels,xinit,vinit);
        else{ 
            std::cerr << "\nInvalid Option entered for Element";
            Config_Status_Local = 4;
            //return Config_Status;
        }
        #pragma omp critical
        {
        std::cout << "\n\t* xinit:\t" << SampleLocal->get_position() << "\n\t* vinit:\t"
            << SampleLocal->get_velocity();
        std::cout << "\n\t* InitRotation:\t" << SampleLocal->get_rotationalfreq() 
            << "\n";
        std::cout << "\n* Matter Object created successfully! *\n\n";
        }
        Pause();
    
    
        // ------------------- PRINT METADATA / CREATE DTOKSU ------------------- //
        MetaDataFile << "\n\n#DUST PARAMETERS" <<"\nElem (arb)\tRadius (m)\t"
            <<"Temp (K)\txinit (m s^-1)\t\tvinit (m s^-1)\n"<<SampleLocal->get_elem()
            <<"\t\t"<<SampleLocal->get_radius()<<"\t\t"<<SampleLocal->get_temperature()
            <<"\t\t"<<SampleLocal->get_position()<<"\t\t"<<SampleLocal->get_velocity()<<"\n";
        DTOKSU* SimLocal;
        PlasmaGrid_Data PgridLocal = Pgrid;
        PlasmaData PdataLocal = Pdata;
        if( !ContinuousPlasma ){
            MetaDataFile <<"\n\n#PLASMA GRID PARAMETERS"
                <<"\nMachine\tIonSpecies\tgridx\tgridz\tgridtheta\tdlx\tdlz\n"
                <<PgridLocal.device<<"\t"<<IonSpecies<<"\t"<<PgridLocal.gridx<<"\t"
                <<PgridLocal.gridz<<"\t"<<PgridLocal.gridtheta<<"\t"<<PgridLocal.dlx<<"\t"
                <<PgridLocal.dlz<<"\n"<<"\nxmin (m)\txmax (m)\tzmin (m)\tzmax (m)\n"
                <<PgridLocal.gridxmin<<"\t\t"<<PgridLocal.gridxmax<<"\t\t"<<PgridLocal.gridzmin
                <<"\t\t"<<PgridLocal.gridzmax << "\n";
            SimLocal = new DTOKSU(AccuracyLevels, SampleLocal, PgridLocal, PdataLocal, WallBound, 
                CoreBound, HeatTerms, ForceTerms, CurrentTerms);
        }else{
            MetaDataFile <<"\n\n#PLASMA DATA PARAMETERS"
                <<"\n\nNn (m^-3)\tNi (m^-3)\tNe (m^-3)\n"
                <<PdataLocal.NeutralDensity<<"\t\t"<<PdataLocal.IonDensity<<"\t\t"
                <<PdataLocal.ElectronDensity
                <<"\n\nTn (K)\t\tTi (K)\t\tTe (K)\n"
                <<PdataLocal.NeutralTemp<<"\t\t"<<PdataLocal.IonTemp<<"\t\t"
                <<PdataLocal.ElectronTemp<<"\t"
                <<"\n\nPvel (m s^-1)\t\tE (V m^-1)\t\tB (T)\n"
                <<PdataLocal.PlasmaVel<<"\t"<<PdataLocal.ElectricField<<"\t"
                <<PdataLocal.MagneticField<<"\n";
            SimLocal = new DTOKSU(AccuracyLevels, SampleLocal, PdataLocal, HeatTerms, ForceTerms,
                CurrentTerms);
        }
        MetaDataFile <<"\n\n##MODEL SWITHES\n#HEATING MODELS\n"
            << "RadiativeCooling:\t" << HeatModels[0] 
            << "\nEvaporativeCooling:\t" << HeatModels[1]
            << "\nNewtonCooling:\t\t" << HeatModels[2] 
            << "\nNeutralHeatFlux:\t" << HeatModels[3]
            << "\nSOMLIonHeatFlux:\t" << HeatModels[4] 
            << "\nSMOMLIonHeatFlux:\t" << HeatModels[5] 
            << "\nDTOKSIonHeatFlux:\t" << HeatModels[6] 
            << "\nDUSTTIonHeatFlux:\t" << HeatModels[7]
            << "\nOMLElectronHeatFlux:\t" << HeatModels[8]
            << "\nPHLElectronHeatFlux:\t" << HeatModels[9]
            << "\nDTOKSElectronHeatFlux:\t" << HeatModels[10] 
            << "\nSOMLNeutralRecomb:\t" << HeatModels[11]
            << "\nSMOMLNeutralRecomb:\t" << HeatModels[12]
            << "\nDTOKSNeutralRecomb:\t" << HeatModels[13]
            << "\nPHLSEE:\t\t\t" << HeatModels[14] 
            << "\nDTOKSSEE:\t\t" << HeatModels[15] 
            << "\nPHLTEE:\t\t\t" << HeatModels[16]
            << "\nDTOKSTEE:\t\t" << HeatModels[17]
            << "\n\n#FORCING MODELS\n"
            << "Gravity:\t\t" << ForceModels[0]
            << "\nLorentz:\t\t" << ForceModels[1] 
            << "\nSOMLIonDrag:\t\t" << ForceModels[2]
            << "\nSMOMLIonDrag:\t\t" << ForceModels[3] 
            << "\nDTOKSIonDrag:\t\t" << ForceModels[4]
            << "\nDUSTTIonDrag:\t\t" << ForceModels[5] 
            << "\nHybridDrag:\t\t" << ForceModels[6]
            << "\nLloydDrag:\t\t" << ForceModels[7]
            << "\nNeutralDrag:\t\t" << ForceModels[8]
            << "\nRocketForce:\t\t" << ForceModels[9]
            << "\n\n#CHARGING MODELS\n"
            << "OMLe:\t\t\t" << ChargeModels[0]
            << "\nPHLe:\t\t\t" << ChargeModels[1]
            << "\nDTOKSe:\t\t\t" << ChargeModels[2]
            << "\nTHSe:\t\t\t" << ChargeModels[3]
            << "\nOMLi:\t\t\t" << ChargeModels[4]
            << "\nMOMLi:\t\t\t" << ChargeModels[5]
            << "\nSOMLi:\t\t\t" << ChargeModels[6]
            << "\nSMOMLi:\t\t\t" << ChargeModels[7]
            << "\nTHSi:\t\t\t" << ChargeModels[8]
            << "\nDTOKSi:\t\t\t" << ChargeModels[9]
            << "\nTEE:\t\t\t" << ChargeModels[10] 
            << "\nTEESchottky:\t\t" << ChargeModels[11]
            << "\nSEE:\t\t\t" << ChargeModels[12]
            << "\nCW:\t\t\t" << ChargeModels[13]
            << "\nMOMLWEM:\t\t" << ChargeModels[14] << "\n";
        MetaDataFile.close();
    
        SimLocal->OpenFiles(FinalDataFilePrefix,i);
        if( ConstModels[4] == 'n' || ConstModels[4] == 'e' ){
            Config_Status_Local = -3;
        }else if( ConstModels[4] == 'r' || ConstModels[4] == 'b' ){
            Config_Status_Local = -2;
        }else{
            Config_Status_Local = 5;
        }
        run_status = Run_Local(SimLocal,Config_Status_Local);
        Pause();
    }
    return run_status;
    //config_message();
    //return Config_Status;
}

//!< Function to configure plasma grid ready for simulation. 
//!< Plasma data is read from plasma_dirname+filename where filename is a 
//!< hard-coded string
int DTOKSU_Manager::configure_plasmagrid(std::string plasma_dirname){
    DM_Debug("  In DTOKSU_Manager::configure_plasmagrid(std::string "
        << "plasma_dirname)\n\n");
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
        Pgrid.gridx = 141;
        Pgrid.gridz = 301;
        Pgrid.gridtheta=0;
        Pgrid.gridxmin = 1.0;
        Pgrid.gridzmin = -1.5;
        Pgrid.gridxmax = 2.4;
        Pgrid.gridzmax = 1.5;
        std::cout << "\n\n\tCalculation for DIII-D" << std::endl;
    }else if(Pgrid.device=='n'){
        Pgrid.gridx = 180;
        Pgrid.gridz = 400;
        Pgrid.gridtheta=0;
        Pgrid.gridxmin = 0.2;
        Pgrid.gridzmin = -2.0;
        Pgrid.gridxmax = 2.0;
        Pgrid.gridzmax = 4.0;
        std::cout << "\n\n\tCalculation for Double Null MAST 17839files shot" 
            << std::endl;
    }else if(Pgrid.device=='p'){
        Pgrid.gridx = 64;
        Pgrid.gridz = 20;
        Pgrid.gridtheta=64;
        Pgrid.gridxmin = 0.0;
        Pgrid.gridzmin = 0.0;
        Pgrid.gridxmax = 0.15;
        Pgrid.gridzmax = 1.9;
//      Pgrid.gridzmax = 1.0;
        std::cout << "\n\n\t* Calculation for Magnum-PSI *\n";
    }else if(Pgrid.device=='e'){
        Pgrid.gridx = 141;
        Pgrid.gridz = 301;
        Pgrid.gridtheta=0;
        Pgrid.gridxmin = 1.0;
        Pgrid.gridzmin = -1.5;
        Pgrid.gridxmax = 2.4;
        Pgrid.gridzmax = 1.5;
        std::cout << "\n\n\tCalculation for EAST" << std::endl;
    }else if(Pgrid.device=='t'){
        Pgrid.gridx = 120;
        Pgrid.gridz = 280;
        Pgrid.gridtheta=0;
        Pgrid.gridxmin = 0.2;
        Pgrid.gridzmin = -2.0;
        Pgrid.gridxmax = 1.4;
        Pgrid.gridzmax = 0.8;
        std::cout << "\n\n\tCalculation for TEST plasma" << std::endl;
    }else{ 
        std::cout << "Invalid tokamak" << std::endl;
    }

    //impurity.open("output///impurity///impurity.vtk");
    int readstatus(-1);
    try{ //!< Read data
        readstatus = read_data(plasma_dirname);
        if( readstatus == 1)
            std::cerr << "\n* Warning! Some PlasmaGrid values exceed Overflows "
                << "or Underflows! *\n\n";
        else if(readstatus == 2)
            throw PlasmaFileReadFailure();

    }catch( PlasmaFileReadFailure &e ){
        std::cout << e.what();
        std::cerr << "\nRead error status: " << readstatus;
        return readstatus;
    }
    return 0;
}

//!< Function to configure boundary ready for simulation. 
//!< Wall data and core data are read from files in 'dirname/filename'
//!< into BoundaryData&BD. Data files must be polygons (contain more than 2 
//!< points)
int DTOKSU_Manager::configure_boundary(std::string dirname,
std::string filename, Boundary_Data& BD){
    DM_Debug("  In DTOKSU_Manager::configure_boundary(std::string dirname, "
        << "std::string filename, Boundary_Data& BD)\n\n")
    std::ifstream BoundaryGrid_File;
    double R_temp(0.0), Z_temp(0.0);
    char Dummy;
    BoundaryGrid_File.open(dirname+filename);
    assert( BoundaryGrid_File.is_open() );

//    std::cout << "\n" << dirname+filename;
    while( BoundaryGrid_File >> R_temp >> Dummy >>  Z_temp ){
        BD.Grid_Pos.push_back( std::make_pair(R_temp,Z_temp) );
//      std::cout << "\n" << R_temp << "\t" << Dummy << "\t" << Z_temp;
        if( R_temp <= 0.0 ){
            std::cerr << "\nError reading boundary data in DTOKSU_Manager::"
                << "configure_boundary! R_temp < 0";
            return 1;
        }
    }
    assert( BD.Grid_Pos.size() > 2 );
    BoundaryGrid_File.close();
    return 0;
}


int DTOKSU_Manager::read_data(std::string plasma_dirname){
    P_Debug("\tIn DTOKSU_Manager::read_data(std::string plasma_dirname)\n\n");
    // Preallocate size of vectors
    Pgrid.Te = std::vector<std::vector<double>>
        (Pgrid.gridx,std::vector<double>(Pgrid.gridz));
    Pgrid.Ti  = Pgrid.Te;
    Pgrid.Tn  = Pgrid.Te;
    Pgrid.Ta  = Pgrid.Te;
    Pgrid.na0 = Pgrid.Te;
    Pgrid.na1 = Pgrid.Te;
    Pgrid.na2 = Pgrid.Te;
    Pgrid.po  = Pgrid.Te;
    Pgrid.ua0 = Pgrid.Te;
    Pgrid.ua1 = Pgrid.Te;
    Pgrid.bx  = Pgrid.Te;
    Pgrid.by  = Pgrid.Te;
    Pgrid.bz  = Pgrid.Te;
    Pgrid.x   = Pgrid.Te;
    Pgrid.z   = Pgrid.Te;
    Pgrid.dm  = Pgrid.Te;
    Pgrid.gridflag  = std::vector<std::vector<int>>
        (Pgrid.gridx,std::vector<int>(Pgrid.gridz));
    int ReStat = 0;
    if(Pgrid.device=='p'){ //!< Note, grid flags will be empty 
        #ifdef NETCDF_SWITCH
        return read_MPSIdata(plasma_dirname);
        #else
        std::cerr << "\nNETCDF SUPPORT REQUIRED FOR MPSI DATA!";
        std::cerr << "\nRECOMPILE WITH NETCDF AND DEFINE NETCDF_SWITCH!\n\n";
        ReStat = 2;
        return 2;
        #endif
    }else{
        std::ifstream scalars,threevectors,gridflagfile;
	std::cout << "\n* Reading plasma data from " << plasma_dirname << "b2processed.dat *\n\n";
        scalars.open(plasma_dirname+"b2processed.dat");
        threevectors.open(plasma_dirname+"b2processed2.dat");
        gridflagfile.open(plasma_dirname+"locate.dat");

        //!< Open files to read
        assert( scalars.is_open() );
        assert( threevectors.is_open() );
        assert( gridflagfile.is_open() );

        //!< Throw away variables which read in data which is unimportant.
        std::string dummy_line;
        double dummy_dub;
        //!< Ignore first line of file
        getline(scalars,dummy_line);
        getline(threevectors,dummy_line);
        //!< Now loop over the grid and feed in the data into the vectors
        double convertJtoK = 7.242971666667e22;
        double converteVtoK = 11604.5250061657;
        for(unsigned int k=0; k<=Pgrid.gridz-1; k++){
            for(unsigned int i=0; i<=Pgrid.gridx-1; i++){
                //!< This is the read-in format without neutrals
//              scalars >> Pgrid.x[i][k] >> Pgrid.z[i][k] >> Pgrid.Te[i][k]
//                    >> Pgrid.Ti[i][k] >> Pgrid.na0[i][k] >> Pgrid.na1[i][k] 
//                    >> Pgrid.po[i][k] >> Pgrid.ua0[i][k] >> Pgrid.ua1[i][k];
                //!< This is the read-in format with neutrals
                scalars >> Pgrid.x[i][k] >> Pgrid.z[i][k] >> Pgrid.Ti[i][k]
                    >> Pgrid.Te[i][k] >> Pgrid.Tn[i][k] >> Pgrid.na0[i][k] 
                    >> Pgrid.na1[i][k] >> Pgrid.na2[i][k] >> Pgrid.po[i][k] 
                    >> Pgrid.ua0[i][k] >> Pgrid.ua1[i][k] >> dummy_dub;
//                std::cout << "\n" << Pgrid.x[i][k] << "\t" << Pgrid.z[i][k] << "\t" << Pgrid.Te[i][k]
//                    << "\t" << Pgrid.Ti[i][k] << "\t" << Pgrid.na0[i][k] << "\t" << Pgrid.na1[i][k]
//                    << "\t" << Pgrid.po[i][k] << "\t" << Pgrid.ua0[i][k] << "\t" << Pgrid.ua1[i][k];

                threevectors >> dummy_dub >> dummy_dub >> Pgrid.bx[i][k] >>
                    Pgrid.bz[i][k] >> Pgrid.by[i][k];
                gridflagfile >> dummy_dub >> dummy_dub >> Pgrid.gridflag[i][k];

                //!< For some reason, very small non-zero values are being 
                //!< assigned to the 'zero' values being read in.
                //!< This should correct for this...
                if( Pgrid.device == 'j' || Pgrid.device == 'i' || Pgrid.device == 'd' ){ 
	                Pgrid.Te[i][k] = convertJtoK*Pgrid.Te[i][k];
        	        Pgrid.Ti[i][k] = convertJtoK*Pgrid.Ti[i][k];
                	Pgrid.Tn[i][k] = convertJtoK*Pgrid.Tn[i][k];
                }else{
                        Pgrid.Te[i][k] = converteVtoK*Pgrid.Te[i][k];
                        Pgrid.Ti[i][k] = converteVtoK*Pgrid.Ti[i][k];
                        Pgrid.Tn[i][k] = converteVtoK*Pgrid.Tn[i][k];
		}
                if( fabs(Pgrid.bx[i][k]) > Overflows::Field 
                    || fabs(Pgrid.bx[i][k]) < Underflows::Field ){ ReStat = 1; } 
                if( fabs(Pgrid.by[i][k]) > Overflows::Field 
                    || fabs(Pgrid.by[i][k]) < Underflows::Field ){ ReStat = 1; } 
                if( fabs(Pgrid.bz[i][k]) > Overflows::Field 
                    || fabs(Pgrid.bz[i][k]) < Underflows::Field ){ ReStat = 1; }
                if( Pgrid.Te[i][k] > Overflows::Temperature 
                    || Pgrid.Te[i][k] < Underflows::Temperature ){ ReStat = 1; } 
                if( Pgrid.Ti[i][k] > Overflows::Temperature 
                    || Pgrid.Ti[i][k] < Underflows::Temperature ){ ReStat = 1; } 
                if( Pgrid.na0[i][k] > Overflows::Density 
                    || Pgrid.na0[i][k] < Underflows::Density )   { ReStat = 1; }
                if( Pgrid.na1[i][k] > Overflows::Density 
                    || Pgrid.na1[i][k] < Underflows::Density )   { ReStat = 1; }
                if( fabs(Pgrid.ua0[i][k]) > Overflows::PlasmaVel 
                    || fabs(Pgrid.ua0[i][k]) < Underflows::PlasmaVel ) 
                    { ReStat = 1; }
                if( fabs(Pgrid.ua1[i][k]) > Overflows::PlasmaVel 
                    || fabs(Pgrid.ua1[i][k]) < Underflows::PlasmaVel ) 
                    { ReStat = 1; }
                Pgrid.Ta[i][k] = Pdata.AmbientTemp;
                Pgrid.Tn[i][k] = Pdata.NeutralTemp;
                Pgrid.na2[i][k] = Pdata.NeutralDensity;
                Pgrid.dm[i][k] = 0.0;
            }
        }
        scalars.close();
        threevectors.close();
        gridflagfile.close();
    }
    return ReStat; //!< return success!
}
// *************************** READING FUNCTIONS *************************** //

//!< for Magnum-PSI, we need to read a NET-cdf file which is special
//!< This function does all the necessary effort of extracting this information.
#ifdef NETCDF_SWITCH
int DTOKSU_Manager::read_MPSIdata(std::string plasma_dirname){
    P_Debug("\tDTOKSU_Manager::read_MPSIdata(std::string plasma_dirname)\n\n");
    
    const int NC_ERR = 2;
    std::string filename 
    = "Magnum-PSI_Experiment_Homogeneous-B-Field_B0.4_L1.9.nc";
    std::cout << "\t\t* Reading data from: " << plasma_dirname + filename;
    std::cout << " *\n\t\t* gridx: " << Pgrid.gridx << "\t * gridz: " 
        << Pgrid.gridz << "\t * gridtheta: " << Pgrid.gridtheta << " *\n";

    float electron_dens_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
    float electron_temp_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
    float electron_Vele_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
    float ion_dens_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
    float ion_Veli_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
    float Potential_mat[Pgrid.gridx][Pgrid.gridz][Pgrid.gridtheta];
    float Bxy_mat[Pgrid.gridx][Pgrid.gridz];

    //< Change the error behavior of the netCDF C++ API by creating an
    //< NcError object. Until it is destroyed, this NcError object will
    //< ensure that the netCDF C++ API silently returns error codes on
    //< any failure, and leaves any other error handling to the calling
    //< program. In the case of this example, we just exit with an
    //< NC_ERR error code.
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
//  if (!(B_xy = dataFile.get_var("B_xy")))
//      return NC_ERR;

    if (!Ne->get(&electron_dens_mat[0][0][0], 
        Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
        return NC_ERR;
    if (!e_Temp->get(&electron_temp_mat[0][0][0], 
        Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
        return NC_ERR;
    if (!Vel_e->get(&electron_Vele_mat[0][0][0], 
        Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
        return NC_ERR;
    if (!Ni->get(&ion_dens_mat[0][0][0], 
        Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
        return NC_ERR;
    if (!Vel_i->get(&ion_Veli_mat[0][0][0], 
        Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
        return NC_ERR;
    if (!Potential->get(&Potential_mat[0][0][0], 
        Pgrid.gridx, Pgrid.gridz, Pgrid.gridtheta))
        return NC_ERR;
//  if (!B_xy->get(&Bxy_mat[0][0], Pgrid.gridx, Pgrid.gridz))
//      return NC_ERR;
    double ConvertJtoK(7.24297166e22); //!< Conversion factor from J to K
    for(unsigned int i=0; i< Pgrid.gridx; i++){
        for(unsigned int k=0; k< Pgrid.gridz; k++){

//            std::cout << "\n[" << i << "]" << "[" << k << "]";  std::cin.get();

            Pgrid.Ti[i][k] = fabs(electron_temp_mat[i][k][0])*ConvertJtoK;// Can't set Ti=fixed as this messes up the Dust potential! 2.5*11598.5895;
            Pgrid.Te[i][k] = fabs(electron_temp_mat[i][k][0])*ConvertJtoK;
            Pgrid.Tn[i][k] = Pdata.NeutralTemp;
            Pgrid.Ta[i][k] = Pdata.AmbientTemp;
            Pgrid.na0[i][k] = ion_dens_mat[i][k][0];
            Pgrid.na1[i][k] = electron_dens_mat[i][k][0];
            Pgrid.na2[i][k] = Pdata.NeutralDensity;
            Pgrid.po[i][k] = Potential_mat[i][k][0];
            Pgrid.ua0[i][k] = ion_Veli_mat[i][k][0];
            Pgrid.ua1[i][k] = electron_Vele_mat[i][k][0];
            Pgrid.bx[i][k] = 0.0;
            Pgrid.by[i][k] = 0.0;
            Pgrid.bz[i][k] = 0.4;
            Pgrid.x[i][k] = Pgrid.gridxmin+i*Pgrid.dlx;
            Pgrid.z[i][k] = Pgrid.gridzmin+k*Pgrid.dlz;
            Pgrid.gridflag  = std::vector<std::vector<int>>
                (Pgrid.gridx,std::vector<int>(Pgrid.gridz));
            Pgrid.dm[i][k] = 0.0;
        }
    }
    return 0;
}
#endif

int DTOKSU_Manager::ChargeTest(double accuracy
,std::vector<CurrentTerm*> CurrentTerms){
    DM_Debug("  In DTOKSU_Manager::ChargeTest(double accuracy"
        << ",std::vector<CurrentTerm*> CurrentTerms)\n\n");
    if( Config_Status != -2 && Config_Status != -3 ){
        std::cerr << "\nDTOKSU Is not configured! Please configure first.";
        config_message();
        return 1;
    }
    ChargingModel MyChargeModel("Data/DTOKS_ChargeTest.txt",accuracy,
        CurrentTerms,Sample,Pdata);
    bool cm_InGrid = MyChargeModel.update_plasmadata();

    //!< Need to manually update the first time as first step is not necessarily
    //!< heating
    Sample->update();       
    double ChargeTime   = MyChargeModel.UpdateTimeStep();
    MyChargeModel.Charge(ChargeTime);

}

int DTOKSU_Manager::ForceTest(double accuracy,std::vector<ForceTerm*> ForceTerms){
    DM_Debug("  In DTOKSU_Manager::ForceTest(double accuracy,"
        << " std::vector<ForceTerm*> ForceTerms)\n\n");
    if( Config_Status != -2 && Config_Status != -3 ){
        std::cerr << "\nDTOKSU Is not configured! Please configure first.";
        config_message();
        return 1;
    }
    ForceModel MyForceModel("Data/DTOKS_ForceTest.txt",accuracy,ForceTerms,
        Sample,Pdata);
    bool fm_InGrid = MyForceModel.update_plasmadata();

    //!< Manually update the first time as first step is not necessarily heating
    Sample->update();       
    double ForceTime    = MyForceModel.UpdateTimeStep();
    MyForceModel.Force(ForceTime);
}

int DTOKSU_Manager::HeatTest(double accuracy,std::vector<HeatTerm*> HeatTerms){
    DM_Debug("  In DTOKSU_Manager::HeatTest(double accuracy,"
        << " std::vector<HeatTerm*> HeatTerms)\n\n");
    if( Config_Status != -2 && Config_Status != -3 ){
        std::cerr << "\nDTOKSU Is not configured! Please configure first.";
        config_message();
        return 1;
    }
    HeatingModel MyHeatModel("Data/DTOKS_HeatTest.txt",accuracy,HeatTerms,
        Sample,Pdata);
    bool hm_InGrid = MyHeatModel.update_plasmadata();

    //!< Manually update the first time as first step is not necessarily heating
    Sample->update();       
    double HeatTime     = MyHeatModel.UpdateTimeStep();
    MyHeatModel.Heat(HeatTime);
}

// Run DTOKS Normally a single time
int DTOKSU_Manager::Run(){
    DM_Debug("  In DTOKSU_Manager::Run()\n\n");
    if( Config_Status != -2 && Config_Status != -3 ){
        std::cerr << "\nDTOKSU Is not configured! Please configure first.";
        config_message();
        return 1;
    }

    clock_t begin = clock();    // Measure start time

    // Actually running DTOKS
    int RunStatus(-1);
    if( Config_Status == -3 ){
        std::cout << "\n * RUNNING DTOKS * \n";
        RunStatus = Sim->Run();
    }else if( Config_Status == -2 )
        Breakup();
    else{
        std::cerr << "\nBreakup Is not configured! Please configure correctly.";
        config_message();
        return 1;
    }
    Sim->ImpurityPrint();

    clock_t end = clock();      // Measure end time
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;  
    std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
    std::cout << "\n\n * DTOKS COMPLETED SUCCESSFULLY * \n\n";

    return RunStatus;
}

// Run a DTOKS thread a single time
int DTOKSU_Manager::Run_Local(DTOKSU* SimLocal, int Config_Status_Local){
    DM_Debug("  In DTOKSU_Manager::Run_Local(int Config_Status_Local)\n\n");
    if( Config_Status_Local != -2 && Config_Status_Local != -3 ){
        std::cerr << "\nDTOKSU Is not configured! Please configure first.";
        return 1;
    }

    clock_t begin = clock();    // Measure start time

    // Actually running DTOKS
    int RunStatus(-1);
    if( Config_Status_Local == -3 ){
        std::cout << "\n * RUNNING DTOKS * \n";
        RunStatus = SimLocal->Run();
    }else if( Config_Status_Local == -2 ){
        std::cerr << "\n* Parallelised Breakup not currently Implemented! Running normally. *";
        RunStatus = SimLocal->Run();
    }else{
        std::cerr << "\nBreakup Is not configured! Please configure correctly.";
        return 1;
    }
    SimLocal->ImpurityPrint();

    clock_t end = clock();      // Measure end time
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;  
    std::cout << "\n\n*****\n\nCompleted in " << elapsd_secs << "s\n";
    std::cout << "\n\n * DTOKS COMPLETED SUCCESSFULLY * \n\n";

    return RunStatus;
}


// Run DTOKS many times with breakup turned on
void DTOKSU_Manager::Breakup(){
    DM_Debug("  In DTOKSU_Manager::Breakup()\n\n");

    if( Config_Status != -2 ){
        std::cerr << "\nDTOKSU Is not configured! Please configure first.";
        config_message();
        return;
    }
    
    
    std::cout << "\n * RUNNING DTOKS * \n";

    threevector Zeroes(0.0,0.0,0.0);

    clock_t begin = clock();

    unsigned int p(1);
    unsigned int i(1);
    
    double seed
        =std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 randnumber(seed);
    //!< Uniformly Randomly Distributed Variable between 0.0 and 1.0
    std::uniform_real_distribution<double> rad(0.0, 1.0); 

    //!< Vector of Grain Data to record the final state of each branch
    std::vector<GrainData> GDvector;

    //!< Vectors to retain the time of each simulation branch so that the global
    //!< time is preserved
    std::vector<double> HMTime;
    std::vector<double> FMTime;
    std::vector<double> CMTime;
    // Loop over the number of split paths.
    for( unsigned int j(0); j < i; j ++){
        DM_Debug("\tSimulating NEGATIVE Branch "); DM_Debug(i); DM_Debug(" : ");
        DM_Debug(j);
        DM_Debug("\n\tStart Pos = "); DM_Debug(Sample->get_position()); 
        DM_Debug("\n\tVelocity = "); DM_Debug(Sample->get_velocity());
        DM_Debug("\n\tMass = "); DM_Debug(Sample->get_mass()); 
        DM_Debug("\n\tTemperature = "); DM_Debug(Sample->get_temperature());
        DM_Debug("\n\t");

        //!< When breakup occurs and a path forks, track it. If it breaks up, 
        //!< track the subsequent particle Repeat until the end condition is no-
        //!< longer breakup, i.e while return of DTOKSU object isn't 3.
        while( Sim->Run() == 3 ){ // DTOKSU_Manager has occured...

            //!< Close data files and open new ones, with names based off index
            Sim->CloseFiles();
            Sim->OpenFiles("Data/breakup",p);
            
            //!< Reset breakup so that it's recorded with breakup turned off
            Sample->reset_breakup();

            //!< Reset the end point data with the same position, no rotation 
            //!< and heading off in negative direction
            //!< Rotation occurs in random direction perpendicular to magnetic
            //!< field and velocity as per theory
            double VelocityMag = 2*PI*(Sample->get_radius())*
                Sample->get_rotationalfreq();
            threevector Unit(2.0*rad(randnumber)-1.0,2.0*rad(randnumber)-
                1.0,2.0*rad(randnumber)-1.0);
            threevector VelocityUnitVec 
                = (Unit.getunit()^Sim->get_bfielddir()).getunit();
            double RandomlyDistributeRotationAndLinMom = rad(randnumber);
            threevector dvMinus = RandomlyDistributeRotationAndLinMom*
                VelocityMag*VelocityUnitVec; // This is a fudge
            Sample->update_motion(Zeroes,dvMinus,
                -RandomlyDistributeRotationAndLinMom*
                Sample->get_rotationalfreq()/2.0);

            //!< Record dust end conditions for when other half is simulated
            GDvector.push_back(Sample->get_graindata());
            HMTime.push_back(Sim->get_HMTime());
            FMTime.push_back(Sim->get_FMTime());
            CMTime.push_back(Sim->get_CMTime());
    
            //!< Change dust velocity, mass has already been halved in Matter. 
            //!< Add the velocity twice over as we took it away in one direction
            threevector dvPlus = -2.0*dvMinus;
            Sample->update_motion(Zeroes,dvPlus,0.0);

            DM_Debug("\nSimulating POSITIVE Branch "); 
            DM_Debug(i); DM_Debug(" : "); DM_Debug(j);
            DM_Debug("\nStart Pos = "); DM_Debug(Sample->get_position());
            DM_Debug("\nVelocity = "); DM_Debug(Sample->get_velocity());
            DM_Debug("\nMass = "); DM_Debug(Sample->get_mass()); 
            DM_Debug("\nTemperature = "); 
            DM_Debug(Sample->get_temperature()); DM_Debug("\n");
            //!< increment counters of number of forks and positive forks.
            //!< p is used for recording index of file
            i = i + 1;
            p = p + 1;
        } //!< Run the simulation again if breakup occured!
        
        DM_Debug("\n***** START OF : DUST DIDN'T BREAKUP *****\n!");
        DM_Debug(i); DM_Debug(" : "); DM_Debug(j);
        //!< Close data files and open new ones
        Sim->CloseFiles();
        if( GDvector.size() > 0 ){ // If we have at least one breakup event, 
            //!< Re-initialise simulation with new Velocity... 
            //!< i.e stored end point data But we've already done this now when
            //!< we saved the data previously...
            //!< So we're good to go, just copy over the previous end conditions
            Sample->set_graindata(GDvector[j]);


            //!< Reset Model Times, this is to make the plotting work correctly.
            //!< Ensure that time of the models is global and not local to track
            Sim->ResetModelTime(
                HMTime[j]-Sim->get_HMTime(),
                FMTime[j]-Sim->get_FMTime(),
                CMTime[j]-Sim->get_CMTime());

            //!< Open files 
            Sim->OpenFiles("Data/breakup",p);
            p = p + 1;
            DM_Debug("\nSimulating POSITIVE Branch "); 
            DM_Debug(i); DM_Debug(" : "); DM_Debug(j);
            DM_Debug("\nStart Pos = "); DM_Debug(Sample->get_position()); 
            DM_Debug("\nVelocity = "); DM_Debug(Sample->get_velocity());
            DM_Debug("\nMass = "); DM_Debug(Sample->get_mass());
            DM_Debug("\nTemperature = "); 
            DM_Debug(Sample->get_temperature()); DM_Debug("\n");
            Pause();
        }
        DM_Debug("\n***** END OF : DUST DIDN'T BREAKUP *****\n!");
    }

    //  Pgrid.datadump(); // Print the plasma grid data
}
