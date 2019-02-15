#include "ForceModel.h"
#include "MathHeader.h" // This is weird

int NeutralDragTest(char Element, bool DragSwitch){
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

    // Parameters describing the heating model
    double Radius=1e-6;       // m
    double Temp=280;        // K
    double Potential = 1;       // Normalised Potential
    Matter *Sample;         // Define the sample matter type

    // Set to true all Force models that are wanted
    bool Gravity = true;        // IS ON!
    bool Centrifugal = true;    // IS ON!
    bool Lorentz = true;        // IS ON!
    bool SOMLIonDrag = true;    // IS ON!
    bool SMOMLIonDrag = true;   // IS ON!
    bool DTOKSIonDrag = true;   // IS ON!
    bool DUSTTIonDrag = true;   // IS ON!
    bool HybridDrag = true;   // IS ON!
    bool NeutralDrag = DragSwitch;  // IS ON!
    bool RocketForce = true;   // IS OFF!

    PlasmaData *Pdata = new PlasmaData();
    Pdata->NeutralDensity   = 1e18;     // m^-3, Neutral Density
    Pdata->IonDensity   = 1e18;     // m^-3, Ion Density
    Pdata->ElectronDensity  = 1e18;     // m^-3, Electron Density
    Pdata->NeutralTemp  = 10*1.16e4;    // K, Neutral Temperature, convert from eV
    Pdata->IonTemp      = 10*1.16e4;    // K, Ion Temperature, convert from eV
    Pdata->ElectronTemp = 10*1.16e4;    // K, Electron Temperature, convert from eV
    threevector PlasmaVel(0.0,0.0,0.5*sqrt(Kb*Pdata->NeutralTemp/Mp));
    threevector gravity(0.0,0.0,-9.81);
    Pdata->PlasmaVel    = PlasmaVel;    // ms^-1, Plasma Velocity
    Pdata->Gravity      = gravity;      // ms^-2, acceleration due to gravity
    Pdata->mi           = Mp;           // kg, Ion mass is proton mass 

    std::array<bool,FMN> ForceModels  = {Gravity,Centrifugal,Lorentz,SOMLIonDrag,SMOMLIonDrag,
        DTOKSIonDrag,DUSTTIonDrag,HybridDrag,NeutralDrag,RocketForce};

    // Models and ConstModels are placed in an array in this order:
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
    Sample->set_potential(Potential);
    threevector vinit(2.0,3.0,5.0);
    threevector zeros(1.0,0.0,0.0);
    Sample->update_motion(zeros,vinit,0.0);
    
    // START NUMERICAL MODEL
    std::string filename;
    if( DragSwitch )    filename = "Data/NeutralDragOn.txt";
    else                filename = "Data/NeutralDragOff.txt";
    ForceModel MyModel(filename,0.1,ForceModels,Sample,Pdata);
    double Mass = Sample->get_mass();
    MyModel.UpdateTimeStep();
    size_t imax(100000);

    for( size_t i(0); i < imax; i ++)
        MyModel.Force();

    threevector ModelVelocity = Sample->get_velocity();

    delete Pdata;
    delete Sample;
    // END NUMERICAL MODEL
    
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nErrorEstimateTest 1 completed in " << elapsd_secs << "s\n";

    return 1;
}
