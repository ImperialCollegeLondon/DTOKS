#include "HeatingModel.h"

int EvaporativeCoolingTest(char Element, bool HeatSwitch){
    clock_t begin = clock();
    // ********************************************************** //
    // FIRST, define program default behaviour

    // Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
    char EmissivityModel = 'c';     // Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
    char ExpansionModel = 'c';  // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
                                                    // and (z)ero expansion
    char HeatCapacityModel = 'c';   // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
    char BoilingModel = 'y';    // Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
                                                    // and (t)homson
    char TimeStepType = 'f';    // Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
                                                    // and (f)ixed steps
    std::string Name="constant";    // Describes heating model

    // Parameters describing the heating model
    double Radius=1e-6;     // m
    double Temp=280;        // K
    double Potential = 1;
    Matter *Sample;         // Define the sample matter type

    // Set to true all heating models that are wanted
    bool RadiativeCooling = true;
    bool EvaporativeCooling = HeatSwitch;
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
    Pdata->AmbientTemp = 0;

    // Models and ConstModels are placed in an array in this order:
    std::array<bool, HMN> Models = 
        {RadiativeCooling, EvaporativeCooling, NewtonCooling, NeutralHeatFlux,
            SOMLIonHeatFlux, SOMLNeutralRecombination, SMOMLIonHeatFlux, 
            SMOMLNeutralRecombination, TEE, SEE, PHLElectronHeatFlux,
            OMLElectronHeatFlux, DTOKSTEE, DTOKSSEE, DTOKSIonHeatFlux,
            DTOKSNeutralRecomb, DTOKSElectronHeatFlux, DUSTTIonHeatFlux };
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
    }else if (Element == 'D'){
        Sample = new Deuterium(Radius,Temp,ConstModels);
    }else if (Element == 'M'){
        Sample = new Molybdenum(Radius,Temp,ConstModels);
    }else if (Element == 'L'){
        Sample = new Lithium(Radius,Temp,ConstModels);
    }else{ 
        std::cerr << "\nInvalid Option entered";
        return -1;
    }
    threevector xinit(1.15,0.0,-1.99);// default injection right hand side
    threevector vinit(0.0,0.0,0.0);
    Sample->update_motion(xinit,vinit,0.0);

    Sample->set_potential(Potential);

    std::string filename;
    if( HeatSwitch )    filename = "Data/EvaporativeCoolingOn.txt";
    else            filename = "Data/EvaporativeCoolingOff.txt";

    HeatingModel MyModel(filename,1.0,Models,Sample,Pdata);
    MyModel.UpdateTimeStep();
    MyModel.Vapourise();
    
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;

    std::cout << "\n\n*****\n\nErrorEstimateTest 3 completed in " << elapsd_secs << "s\n";

    return 0.0;
}
