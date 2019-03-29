#include "HeatingModel.h"

int VariableHeatCapacityTest(char Element, bool VaryHeatCapacity){
    clock_t begin = clock();
    // ********************************************************** //
    // FIRST, define program default behaviour

    // Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
    char EmissivityModel = 'c';     // Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
    char ExpansionModel = 'c';  // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
                                                    // and (z)ero expansion

    char HeatCapacityModel = 'c';   // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
    if (VaryHeatCapacity) HeatCapacityModel = 'v';

    char BoilingModel = 'y';    // Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
                                                    // and (t)homson
    char TimeStepType = 'f';    // Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
                                                    // and (f)ixed steps
    std::string Name="constant";    // Describes heating model

    // Parameters describing the heating model
    double Radius=1e-6;       // m
    double Temp=280;        // K
    double TimeStep=1e-9;       // s
    double Potential = 1;
    Matter *Sample;         // Define the sample matter type

    // Set to true all heating models that are wanted
    bool RadiativeCooling = true;
    bool EvaporativeCooling = false;
    bool NewtonCooling = false;     // This model is equivalent to Electron and Ion heat flux terms

    // Plasma heating terms
    bool NeutralHeatFlux = true;
    bool SOMLIonHeatFlux = true;
    bool SOMLNeutralRecombination = true;
    bool SMOMLIonHeatFlux = true;
    bool SMOMLNeutralRecombination = true;

    // Electron Emission terms
    bool TEE = true;
    bool SEE = true;

    bool PHLElectronHeatFlux = true;
    bool OMLElectronHeatFlux = true;
    bool DTOKSTEE = true;
    bool DTOKSSEE = true;
    bool DTOKSIonHeatFlux = true;
    bool DTOKSNeutralRecomb = true;
    bool DTOKSElectronHeatFlux = true;
    bool DUSTTIonHeatFlux = true;

    PlasmaData *Pdata = new PlasmaData;
    Pdata->NeutralDensity   = 1e18;     // m^-3, Neutral Density
    Pdata->IonDensity   = 1e18;     // m^-3, Ion Density
    Pdata->ElectronDensity  = 1e18; // m^-3, Electron Density
    Pdata->NeutralTemp  = 10*1.16e4;    // K, Neutral Temperature, convert from eV
    Pdata->IonTemp      = 10*1.16e4;    // K, Ion Temperature, convert from eV
    Pdata->ElectronTemp = 10*1.16e4;    // K, Electron Temperature, convert from eV
    Pdata->mi           = Mp;        // kg, Ion Mass
    Pdata->A            = 1.0;       // arb, Atomic Number
    Pdata->AmbientTemp = 0;

    // Models and ConstModels are placed in an array in this order:
    std::array<bool, 18> HeatModels = 
        {RadiativeCooling, EvaporativeCooling, NewtonCooling, NeutralHeatFlux,
            OMLElectronHeatFlux, PHLElectronHeatFlux, DTOKSElectronHeatFlux,
            SOMLIonHeatFlux,  SMOMLIonHeatFlux, DTOKSIonHeatFlux, DUSTTIonHeatFlux,
            SOMLNeutralRecombination, SMOMLNeutralRecombination, DTOKSNeutralRecomb,
            SEE, DTOKSSEE, TEE, DTOKSTEE};

    std::vector<HeatTerm*> HeatTerms;
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
    
    std::array<char, CM> ConstModels =
        { EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel,'n'};

    if  (Element == 'W'){ 
        Sample = new Tungsten(Radius,Temp,ConstModels);
    }else if (Element == 'B'){ 
        Sample = new Beryllium(Radius,Temp,ConstModels);
    }else if (Element == 'F'){
        Sample = new Iron(Radius,Temp,ConstModels);
    }else if (Element == 'G'){
        Sample = new Graphite(Radius,Temp,ConstModels);
    }else if (Element == 'M'){
        Sample = new Molybdenum(Radius,Temp,ConstModels);
    }else if (Element == 'L'){
        Sample = new Lithium(Radius,Temp,ConstModels);
    }else if (Element == 'D'){
                Sample = new Deuterium(Radius,Temp,ConstModels);
        }else{ 
        std::cerr << "\nInvalid Option entered";
        return -1;
    }
    threevector xinit(1.15,0.0,-1.99);// default injection right hand side
    threevector vinit(0.0,0.0,0.0);
    Sample->update_motion(xinit,vinit,0.0);

    Sample->set_potential(Potential);

    std::string filename;
    if( VaryHeatCapacity )  filename = "Data/HeatCapacityNonConst.txt";
    else                    filename = "Data/HeatCapacityConst.txt";

    HeatingModel MyModel(filename,1.0,HeatTerms,Sample,Pdata);
    MyModel.UpdateTimeStep();
    MyModel.Vapourise();
    
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

    std::cout << "\n\n*****\n\nErrorEstimateTest 5 completed in " << elapsd_secs << "s\n";

    return 0.0;
}
