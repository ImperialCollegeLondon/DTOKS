/** @file Matter.cpp
 *  @brief Implementation of class for matter object for use in physics models
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 *  @bug update_dim needs validating
 *  @bug update() the set heat capacity option results in a hard coded value.
 */

#include "Matter.h"

//!< Default values for St data structure describe Tungsten. see GrainStructs.h
struct GrainData MatterDefaults = { 
    false,          //!< Is Liquid,         fine
    false,          //!< Is Gas,            fine
    false,          //!< Breakup,           fine
    1e-6,           //!< Unheated Radius,   should be updated
    8.21e-14,       //!< Mass,              should be updated
    1e-6,           //!< Radius,            fine
    1.26e-11,       //!< Surface Area       should be updated
    4.189e-18,      //!< Volume             should be updated
    19600,          //!< Density            should be updated
    5555,           //!< SuperBoilingTemp   should be updated   (was 3000)
    270,            //!< Temperature        fine
    1e5,            //!< Vapour Pressure    should be updated
    1.0,            //!< Emissivity         should be updated
    1.0,            //!< Linear Expansion   should be updated
    0.5,            //!< Heat Capacity      should be updated
    0,              //!< DeltaSec           should be updated
    0,              //!< DeltaTherm         should be updated
    0,              //!< Potential          should be updated
    false,          //!< Is Positive        should be updated
    {0.0,0.0,0.0},  //!< Dust position      should be updated
    {0,0,0},        //!< Dust Velocity      case dependant
    0,              //!< Rotational Freq    fine
    0,              //!< FusionEnergy       fine
    0               //!< VapourEnergy       fine
};

//!< Make sure dimensions of material are self-consistent after constructor
Matter::Matter(const ElementConsts &elementconsts):Ec(elementconsts),
        St(MatterDefaults){
    M_Debug("\n\nIn Matter::Matter(...)\n\n");
    ConstModels = {'c','c','c','y'};
    update_dim();
    PreBoilMass = St.Mass;
};

Matter::Matter(double rad, const ElementConsts &elementconsts):Ec(elementconsts)
        ,St(MatterDefaults){
    M_Debug("\n\nIn Matter::Matter(...)\n\n");
    ConstModels = {'c','c','c','y'};
    St.Radius = rad;            // m
    St.UnheatedRadius = St.Radius;      // m
    update_dim();
    PreBoilMass = St.Mass;
    assert(St.Radius > 0 && St.UnheatedRadius > 0);
};

Matter::Matter(double rad, double temp, const ElementConsts &elementconsts):
        Ec(elementconsts),St(MatterDefaults){
    M_Debug("\n\nIn Matter::Matter(...)\n\n");
    ConstModels = {'c','c','c','y'};
    St.Radius = rad;                        // m
    St.UnheatedRadius = St.Radius;          // m
    St.Volume = (4*PI*pow(St.Radius,3))/3;
    St.SurfaceArea = 4*PI*pow(St.Radius,2);
    St.Temperature = temp;                  // K
    PreBoilMass = St.Mass;
    //!< if the size is not variable, set the initial mass according to density
    if( ConstModels[1] != 'v' && ConstModels[1] != 'V' )    
        St.Mass = Ec.RTDensity*St.Volume;
    else                            
        update();

    //!< Ascertain the initial state dependant on input temperature
    if( St.Temperature >= Ec.BoilingTemp ){
        St.VapourEnergy = Ec.LatentVapour*St.Mass;
        St.FusionEnergy = Ec.LatentFusion*St.Mass;
        St.Liquid = false;
        St.Gas = true;
    }else if( St.Temperature >= Ec.MeltingTemp ){
        St.FusionEnergy = Ec.LatentFusion*St.Mass;
        St.Liquid = true;
        St.Gas = false;
    }else{
        St.Liquid = false;
        St.Gas = false;
    }
    assert(St.Radius > 0 && St.UnheatedRadius > 0 && St.Temperature > 0 
        && St.Temperature < St.SuperBoilingTemp );
    
    M_Debug("\nSt.Radius = " << St.Radius << "\nSt.Mass = " << St.Mass 
        << "\nSt.Temperature = " << St.Temperature);
    M_Debug("\nEc.LatentVapour = " << Ec.LatentVapour << "\nEc.AtomicMass = " 
        << Ec.AtomicMass);
};

Matter::Matter(double rad, double temp, const ElementConsts &elementconsts, 
        std::array<char,CM> &constmodels):Ec(elementconsts),St(MatterDefaults){
    M_Debug("\n\nIn Matter::Matter(...)\n\n");
    ConstModels = constmodels;
    St.Radius = rad;                        // m
    St.UnheatedRadius = St.Radius;          // m
    St.Volume = (4*PI*pow(St.Radius,3))/3;
    St.SurfaceArea = 4*PI*pow(St.Radius,2);
    St.Temperature = temp;                  // K
    PreBoilMass = St.Mass;
    //!< if the size is not variable, set the initial mass according to density
    if( ConstModels[1] != 'v' && ConstModels[1] != 'V' )    
        St.Mass = Ec.RTDensity*St.Volume;
    else
        update();

    //!< Ascertain the initial state dependant on input temperature
    if( St.Temperature >= Ec.BoilingTemp ){
        St.VapourEnergy = Ec.LatentVapour*St.Mass;
        St.FusionEnergy = Ec.LatentFusion*St.Mass;
        St.Liquid = false;
        St.Gas = true;
    }else if( St.Temperature >= Ec.MeltingTemp ){
        St.FusionEnergy = Ec.LatentFusion*St.Mass;
        St.Liquid = true;
        St.Gas = false;
    }else{
        St.Liquid = false;
        St.Gas = false;
    }
    assert(St.Radius > 0 && St.UnheatedRadius > 0 && St.Temperature > 0 
        && St.Temperature < St.SuperBoilingTemp );

    M_Debug("\nSt.Radius = " << St.Radius << "\nSt.Mass = " << St.Mass 
        << "\nSt.Temperature = " << St.Temperature);
    M_Debug("\nEc.LatentVapour = " << Ec.LatentVapour << "\nEc.AtomicMass = " 
        << Ec.AtomicMass);
};


void Matter::update_dim(){
    M_Debug("\tIn Matter::update_dim():\n\n");

    //!< Determine radius and density
    if(ConstModels[1] == 'c' || ConstModels[1] == 'C'){
        St.Density = Ec.RTDensity;              //!< Density is RT density
        //!< Set radius from Mass
        St.Radius = pow((3*St.Mass)/(4*PI*Ec.RTDensity),1./3.); 
        //!< @bug THIS WHOLE PART NEEDS CHECKING...
    }else if(ConstModels[1] == 'v' || ConstModels[1] == 'V'){ 
        update_radius();
        //!< @bug SHOULD THIS BE : pow(2*St.LinearExpansion,3) ???
        St.Density = Ec.RTDensity / pow(St.LinearExpansion,3); 
        M2_Debug("\nUnheatedRadius was = " << St.UnheatedRadius);
        //!< Change the radius according to the amount of mass lost
        // This was Stefans logic, something to do with Differentials...
        St.UnheatedRadius = St.UnheatedRadius*
            (pow((3*St.Mass)/(4*PI*St.Density),1./3.)/St.Radius);
        M2_Debug("\nUnheatedRadius now = " << St.UnheatedRadius);

    }else if(ConstModels[1] == 's' || ConstModels[1] == 's'){
        St.Density = Ec.RTDensity; //!< Fix the density
    }else{
        std::cout << "\nError! In Matter::update_dim(char ConstModels[1])\n"
                                << "Invalid input for Expansion Model.";
        assert( strchr("vVcCsS",ConstModels[1]) );
    }
    St.Volume = (4*PI*pow(St.Radius,3))/3;
    St.SurfaceArea = 4*PI*pow(St.Radius,2);
    St.Mass = St.Density*St.Volume; 

    // Sanity Check, this may be an issue
    //  assert( abs((St.Density - St.Mass/St.Volume)/St.Density) < 0.000001 ); 
    
    M2_Debug("\nSt.Mass = " << St.Mass << "\nSt.Volume = " << St.Volume 
        << "\nSt.Density = " << St.Density << "\nSt.Mass/St.Volume = " 
        << St.Mass/St.Volume);

    if( St.Gas == false ) assert( St.Mass > MinMass );
}

//!< For variable emissivity, see: https://github.com/cfinch/Mie_scattering
void Matter::update_emissivity(){
    M_Debug("\tIn Matter::update_emissivity()\n\n");
    //!< Assume constant emissivity
    if( ConstModels[0] == 'c' || ConstModels[0] == 'C' ){ 
        St.Emissivity = St.Emissivity;
    //!< Get emissivity from file
    }else if( ConstModels[0] == 'f' || ConstModels[0] == 'F' ){ 
        if( St.Radius < 2e-8 || St.Radius > 1e-4 ){
//          std::cout << "\nError! In Matter::update_emissivity().\n"
//              << "Radius = " << St.Radius 
//              << " outside limit of Emissivty model.\n";
            static bool runOnce1 = true;
            WarnOnce(runOnce1,"Radius outside limit of Emissivty model.");
//          assert(St.Radius > 2e-8);
//          assert(St.Radius < 1e-4);
        }else{
            assert(St.Temperature > 275); //!< Check temperature is in table
            assert(St.Temperature <= Ec.BoilingTemp);
            std::string DirName, FileName("/Temp_");
            std::string TempString, txtstring(".txt");
            //!< Set DirName for correct elements Emissivity Data folder
            if(Ec.Elem == 'w' || Ec.Elem == 'W' ) 
                DirName = "EmissivityData/EmissivityDataTungsten";
            else if(Ec.Elem == 'f' || Ec.Elem == 'F' ) 
                DirName = "EmissivityData/EmissivityDataIron";
            else if(Ec.Elem == 'g' || Ec.Elem == 'G' ) 
                DirName = "EmissivityData/EmissivityDataGraphite";
            else if(Ec.Elem == 'b' || Ec.Elem == 'B' ) 
                DirName = "EmissivityData/EmissivityDataBeryllium";
            std::ifstream TempFile;
            std::stringstream ss;

            ss << round(St.Temperature);
            ss >> TempString;

            //!< Open the file with emissivity data
            TempFile.open((DirName+FileName+TempString+txtstring).c_str());
            bool Found = false;

            char buffer;
            double TestRad;
            //!< Extract data from file
            while( TempFile >> TestRad >> buffer >> St.Emissivity 
                && Found != true ){
                if(TestRad == round_to_digits(St.Radius,2)){ 
                    Found = true;
                    M_Debug( "\nround(St.Radius,9) = " 
                        << round_to_digits(St.Radius,2) << "\nTestRad = " 
                        << TestRad );
                    //Pause();
                }
            }
            ss.clear(); ss.str("");
            TempFile.close(); //!< close the file with emissivity data
        }
    }else{ //!< Invalid input for Emissivity model
        std::cout << "\nError! In Matter::update_emissivity()\n"
                << "Invalid input for Emissivity Model.";
        assert( strchr("cCfF",ConstModels[0]) );
    }
    //!< Ensure emissivity is correctly set
    if(St.Emissivity > 1.0 || St.Emissivity < 0.0 ){
        static bool runOnce2 = true;
        WarnOnce(runOnce2,"Emissivity > 1! Emissivity being forced equal to 1");
        St.Emissivity = 1.0;
    }
}

void Matter::update_boilingtemp(){
    M_Debug("\tIn Matter::update_boilingtemp()\n\n");
    // Conver LatentVapour from kJ to J for the Gas constant
    //  St.CapillaryPressure = exp((St.LatentVapour*St.AtomicMass)/
    //    (1000*R)*(1/St.Temperature() - 1/(St.MeltingTemp)));

    double CapillaryPressure = Ec.SurfaceTension*2/(St.Radius*101325);
    M2_Debug("\n\nCapillaryPressure = " << CapillaryPressure << "atmospheres.");
    if( ConstModels[3] == 's' || ConstModels[3] == 'S' ){ //!< For super boiling
        //!< We can't do Graphite, it has no vapour pressure
        assert(Ec.Elem != 'g' && Ec.Elem != 'G'); 
        St.SuperBoilingTemp = 1/(1/Ec.BoilingTemp 
            - R*log(CapillaryPressure)/(1000*Ec.LatentVapour*Ec.AtomicMass));
        if(Ec.BoilingTemp > St.SuperBoilingTemp){
            St.SuperBoilingTemp = Ec.BoilingTemp;
        }
    //!< For Thomson Boiling
    }else if( ConstModels[3] == 't' || ConstModels[3] == 'T' ){ 
        //!< We can't do Graphite, it has no vapour pressure
        assert(Ec.Elem != 'g' && Ec.Elem != 'G'); 
        //!< Derived from Kelvin and Clausius-Clapyeron Equation, 
        //!< makes a difference of ~5% for Grains of order 1e-8m.
        St.SuperBoilingTemp = Ec.BoilingTemp*exp(-(2*Ec.SurfaceTension)/
            (St.Radius*Ec.LatentVapour*1000*St.Density));
        //!< Kelvin Equation, Makes a difference of ~ 10% for Grains of size 
        //!< 1e-8m. This is an increasingly small effect for larger grains
        St.VapourPressure = St.VapourPressure*
            exp((2*Ec.SurfaceTension*Ec.AtomicMass)
            /(St.Density*St.Radius*R*St.Temperature));
    }else if( ConstModels[3] == 'y' || ConstModels[3] == 'Y' ) 
        St.SuperBoilingTemp = Ec.BoilingTemp;
    else if( ConstModels[3] == 'n' || ConstModels[3] == 'N' ) 
        St.SuperBoilingTemp = std::numeric_limits<double>::max();
    else{
        std::cout << "\nError! In Matter::update_boilingtemp()\n"
            << "Invalid input '" << ConstModels[3] << "' for Emissivity Model.";
        assert( strchr("yYnNsStT",ConstModels[3]) );
    }
    M2_Debug("\nR*log(CapillaryPressure)/(1000*Ec.LatentVapour*Ec.AtomicMass) "
        << "= " << R*log(CapillaryPressure)/(1000*Ec.LatentVapour*Ec.AtomicMass) 
        << "\nSt.SuperBoilingTemp = " << St.SuperBoilingTemp << "\n");
}

void Matter::update_state(double EnergyIn){
    M_Debug("\tIn Matter::update_state(double EnergyIn)\n\n");

    //!< HeatCapacity in Units of kJ/(kg K)
    St.Temperature += EnergyIn/(St.Mass*St.HeatCapacity); 
    double Epsilon(1e-15);  //!< This factor avoids numerical errors.

    //!< Account for latent heat of melted solid that needs to be re-frozen
    if( St.FusionEnergy > 0 && St.Temperature < Ec.MeltingTemp 
        && EnergyIn < 0 ){
        St.FusionEnergy += EnergyIn;
        St.Temperature = Ec.MeltingTemp;
        //!< If it has just re-frozen, St.Fusion Energy is negative, account for 
        //!< lost energy
        if( St.FusionEnergy < Epsilon ){ 
            St.Temperature = Ec.MeltingTemp+(St.FusionEnergy)
            /(St.Mass*St.HeatCapacity);
            St.FusionEnergy = 0;
            St.Liquid = false;
        }
    }
    //!< Case of Melting or Liquid
    if( St.Temperature < St.SuperBoilingTemp 
        && St.Temperature >= Ec.MeltingTemp ){ 
        //!< Factor of 1e-10 is to avoid numerical errors in comparison

        if( St.FusionEnergy < Ec.LatentFusion*St.Mass 
            && St.Temperature < Ec.MeltingTemp + 1.0 ){ //!< Must be melting
            //!< Add energy to Latent heat and set Temperature to Melting 
            //!< Temperature
            St.FusionEnergy += EnergyIn;
            St.Temperature = Ec.MeltingTemp;
            M2_Debug("\n*Liquid is melting*, T = " << St.Temperature
                << "\nFusionEnergy is < LatentFusion*Mass : " << St.FusionEnergy 
                << " < " << Ec.LatentFusion*St.Mass << "\nSt.Mass = " 
                << St.Mass);
            //!< if it melts fully
            if( St.FusionEnergy >= Ec.LatentFusion*St.Mass ){ 
                St.Liquid = true; St.Gas = false;
                St.Temperature = Ec.MeltingTemp + 
                        (St.FusionEnergy-Ec.LatentFusion*St.Mass)
                        /St.HeatCapacity;
            }
            M2_Debug( "\nEnergyIn = " << EnergyIn << "\nSt.FusionEnergy = " 
                << St.FusionEnergy << "\nEc.LatentFusion = " << Ec.LatentFusion 
                << "\nEc.LatentFusion*Mass = " << Ec.LatentFusion*St.Mass 
                << "\nMass = " << St.Mass);
        //!< Else it is in the liquid phase!
        }else{ St.Liquid = true; St.Gas = false; }  
    }else if( St.Temperature >= St.SuperBoilingTemp ){ //!< Boiling or Gas
        if( St.VapourEnergy == 0 ) PreBoilMass = St.Mass;
        //!< Case of boiling
        if( St.VapourEnergy < Ec.LatentVapour*PreBoilMass ){ 
            //!< Add energy to Latent heat and set Temperature to Melting 
            //!< Temperature
            M2_Debug( "\nEnergyIn = " << EnergyIn << "\nSt.VapourEnergy = " 
                << St.VapourEnergy << "\nEc.LatentVapour = " << Ec.LatentVapour 
                << "\nEc.LatentVapour*PreBoilMass = " 
                << Ec.LatentVapour*PreBoilMass << "\nMass = " << St.Mass 
                << "\nPreBoilMass = "  << PreBoilMass);

            update_mass(EnergyIn/Ec.LatentVapour);
            St.VapourEnergy += EnergyIn; 
            St.Temperature = St.SuperBoilingTemp;

            M2_Debug("\n*Liquid is boiling*, T = " << St.Temperature 
                << " K\nVapourEnergy is < St.LatentVapour*Mass : " 
                << St.VapourEnergy << " < " << Ec.LatentVapour*St.Mass);
            //!< Check if it has boiled
            if( St.VapourEnergy > (Ec.LatentVapour*Ec.RTDensity*PI*
                pow(St.UnheatedRadius,3.0)*4.0)/3.0 ){
                St.Liquid = false; St.Gas = true; 
                            St.Temperature = St.SuperBoilingTemp + 
                                (St.VapourEnergy-Ec.LatentVapour*St.Mass)
                                /St.HeatCapacity;
                std::cout << "\n\n***** Sample has Boiled THIS STEP! *****\n";
            }
        }else{  //!< Else it has vapourised!
            St.Liquid = false; St.Gas = true; 
            std::cout << "\n\n***** WARNING! SAMPLE IS GASEOUS *****\n";
        }
    }

    if( St.Mass <= MinMass ){ //!< Lower limit for mass of dust
        St.Liquid = false;
        St.Gas = true;
    }

    M2_Debug( "\n**** Temp change = " << EnergyIn/(St.Mass*St.HeatCapacity) 
        << "\n**** Ein = " << EnergyIn << "\n**** Mass = " << St.Mass 
        << "\n**** Cv = " << St.HeatCapacity << "\nSt.Temperature = " 
        << St.Temperature);
    M2_Debug( "\n**** St.Liquid = " << St.Liquid << "\n**** St.Gas = " 
        << St.Gas);

    assert(St.Temperature > 0);
}

void Matter::update(){
    M_Debug("\tIn Matter::update()\n\n");
    //!< ORDER DEPENDENT UPDATES:
    M_Debug("\t"); update_dim();
    
    //!< ORDER INDEPENDENT UPDATES:
    if( ConstModels[2] == 'v'|| ConstModels[2] == 'V')  update_heatcapacity();
    else if( ConstModels[2] == 'c'||ConstModels[2] == 'C')  
        St.HeatCapacity = St.HeatCapacity;
    else if( ConstModels[2] == 's'|| ConstModels[2] == 'S') 
        St.HeatCapacity = 0.56; //!< This is an arbitrary fixed number
    else{
        std::cout << "\nError! In Matter::update(...)\n"
            << "Invalid input for HeatCapacity Model.";
        assert( strchr("vVcCsS",ConstModels[2]) );  
    }

    //!< Calculate the critical dimensionless angular velocity for rotational
    //!< breakup condition
    double CriticalRotation = 0.56*sqrt(8*Ec.SurfaceTension
        /(St.Density*pow(St.Radius,3)));
    double FudgeFactor = 0.6;
    M2_Debug("\nCriticalRotation = " << CriticalRotation << "\n");
    
    //!< Determine if rotational breakup has occured, if it's a liquid and the
    //!< model is turned on
    if( St.RotationalFrequency > FudgeFactor*CriticalRotation && St.Liquid 
        && ( ConstModels[4] == 'r' || ConstModels[4] == 'b' ) ){
        std::cout << "\nRotational Breakup has occured!";

        if( St.Mass/2 < MinMass ){
            St.Gas = true;
            std::cout << "\nMass too small to support breakup, vaporisation"
                << " assumed...";
        }else{
            St.Mass = St.Mass/2;    // Assume mass is half of previous...
            St.Breakup = true;
            std::cout << "\nMass has split in two...";
        }
    } 
    //!< Call update functions for models with variable constants
    M_Debug("\t"); update_emissivity();
    M_Debug("\t"); update_vapourpressure();
    M_Debug("\t"); update_boilingtemp();
};

void Matter::update_models(char emissivmodel, char linexpanmodel, 
    char heatcapacitymodel, char boilingmodel, char breakupmodel){
    M_Debug("\tIn Matter::update_models(...)\n\n");
    ConstModels[0] = emissivmodel;
    ConstModels[1] = linexpanmodel;
    ConstModels[2] = heatcapacitymodel;
    ConstModels[3] = boilingmodel;
    ConstModels[4] = breakupmodel;
};

void Matter::update_models(std::array<char,CM> &constmodels){
    M_Debug("\tIn Matter::update_models(std::array<char,CM> &constmodels)\n\n");
    ConstModels = constmodels;
}

void Matter::update_mass(double LostMass){ 
    M_Debug("\tIn Matter::update_mass(double LostMass)\n\n");
    M2_Debug("\nMassLoss = " << LostMass  << "\nSt.Mass = " << St.Mass 
        << "\nRadius = " << St.Radius << "\nDensity = " << St.Density);

    //!< Subtract lost mass
    St.Mass -= LostMass;

    //!< If it's a liquid, account for the lost latent energy
    if( St.Liquid ){
        St.FusionEnergy = Ec.LatentFusion*St.Mass;
    }

    //Pause();
    //!< Check if it has fallen below lower mass limit
    if( St.Mass <= MinMass ){ //!< Lower limit for mass of dust
        St.Liquid = false;
        St.Gas = true;
    }
}

void Matter::update_temperature(double EnergyIn){

    M_Debug("\tIn Matter::update_temperature(double EnergyIn = " << EnergyIn 
        << " kJ)\n\n");

    //!< Assert temperature changes smaller than 10% of current temperature
    if(abs(EnergyIn/(St.Mass*St.HeatCapacity)) > St.Temperature*0.1){
        M2_Debug("\n\nSt.Temperature = " << St.Temperature << "\nSt.Mass = " 
            << St.Mass << "\nSt.HeatCapacity = " << St.HeatCapacity 
            << "\nEnergyIn = " << EnergyIn);
        std::cerr << "\nEnergyIn = " << EnergyIn << "\nSt.Mass = " << St.Mass 
            << "\nSt.HeatCapacity = " << St.HeatCapacity 
            << "\nError, EnergyIn/(St.Mass*St.HeatCapacity) = " 
            << EnergyIn/(St.Mass*St.HeatCapacity) << "\nTemperature = " 
            << St.Temperature;
    }

    //!< Assert temperature change is less than 10%
    assert(abs(EnergyIn/(St.Mass*St.HeatCapacity)) < St.Temperature*0.1); 

    //!< Update state
    M_Debug("\t"); update_state(EnergyIn);
    M2_Debug("\nSt.Temperature = " << St.Temperature);

    //!< If it's not a gas, check temperature is within acceptable range
    if(!St.Gas) 
        assert(St.Temperature <= St.SuperBoilingTemp 
            && St.Temperature >= 0.0 );
};

void Matter::update_motion(const threevector &ChangeInPosition,
    const threevector &ChangeInVelocity, double Rotation){
    M_Debug("\tIn Matter::update_motion(...)\n\n");
    
    //!< Calculate new position
    St.DustPosition = St.DustPosition + ChangeInPosition;

    //!< Check that position is sensible
    if( St.DustPosition.getx() == 0.0 ){
        static bool runOnce = true;
        std::string Warning = "Dust Radial Position <= 0.0! Angular";
        Warning += " position poory defined!";
        WarnOnce(runOnce,Warning);
    }

    //!< Update linear and rotational velocity
    St.DustVelocity = St.DustVelocity + ChangeInVelocity;
    St.RotationalFrequency = St.RotationalFrequency + Rotation;
};

void Matter::update_charge(double charge, double potential, double deltat, 
    double deltas){
    M_Debug("\tIn Matter::update_charge(double potential, ...\n\n");
    St.Potential = potential;
    St.DeltaTherm = deltat;
    St.DeltaSec = deltas;

    //!< If the grain is in fact positive ...
    if ( (St.DeltaSec + St.DeltaTherm) >= 1.0 || St.Potential < 0.0 ){ 
        St.Positive = true;
    }else{
        St.Positive = false;
    }
    
    //!< Determine whether electrostatic breakup has occured
    if( fabs(charge) > 8*PI*sqrt(epsilon0*Ec.SurfaceTension*pow(St.Radius,3)) 
        && St.Liquid && (ConstModels[4] == 'e' || ConstModels[4] == 'b') ){
        std::cout << "\nElectrostatic Breakup! Instantaneous vaporisation"
            << " assumed.";
        St.Gas = true;
    }
}
