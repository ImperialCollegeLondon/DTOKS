/** @file Lithium.cpp
 *  @brief Implementation of Lithium class with element constants
 *
 *  Constructors for beryllium class and element specific dependencies of heat
 *  capacity, vapour pressure and size on temperature.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#include "Lithium.h"
#include "Constants.h"

/** @brief Constant structure defining true physical constants of Lithium
 */
const struct ElementConsts LithiumConsts = {
    'L',      //!< Specifies the element
    453.65,   //!< K, Melting temperature at atmospheric pressure
    1603.0,   //!< K, Boiling temperature at atmospheric pressure
    520.0,    //!< kJ/mol, Bond Energy, 
              //!< Estimated as being equal to Latent Vapour Energy
    2.9,      //!< ev, Work Function, 
              //!< Taken from old DTOKS and matched with wikipedia
    1.0,      //!< kJ/(m^2 K), HeatTransfer coefficient, 
              //!< http://www.engineeringtoolbox.com/
              //!< overall-heat-transfer-coefficients-d_284.html
    0.006941, //!< AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
    432.28,   //!< kJ/kg, Latent Fusion Energy From Wikipedia
    19593.71, //!< kJ/kg, Latent Vapour Energy From Wikipedia
    0.35,     //!< N/m, Surface Tension B. Keene, 1993, 
              //!< Review of data for the surface tension of pure metals
    534,      //!< (kg/m^3) from Wikipedia, Room temperature density
    0.0848    //!< kW/m K at 20 degrees celsius
};

Lithium::Lithium():Matter(LithiumConsts){
    E_Debug("\n\nIn Lithium::Lithium():Matter(&LithiumConsts)\n\n");
    set_defaults();
    update();
}

Lithium::Lithium(double radius):Matter(radius,LithiumConsts){
    E_Debug("\n\nIn Lithium::Lithium(double radius):..\n\n");
    set_defaults();
    update();
}

Lithium::Lithium(double radius, double tempin)
    :Matter(radius,tempin,LithiumConsts){
    E_Debug("\n\nIn Lithium::Lithium(double radius, double tempin):...\n\n");
    set_defaults();
    update_state(0.0);
    update_models('c','c','c','y','n');
    update();
    E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius 
        << "\nSt.Density = " << St.Density << "\nSt.Volume = " << St.Volume);
}

Lithium::Lithium(double radius, double tempin, std::array<char,CM> &constmodels)
    :Matter(radius,tempin,LithiumConsts){
    E_Debug("\n\nIn Lithium::Lithium(double radius, double tempin, ...)\n\n\t");
    set_defaults();

    update_state(0.0);
    update_models(constmodels);
    update();
    E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius 
        << "\nSt.Density = " << St.Density << "\nSt.Volume = " << St.Volume);
}

Lithium::Lithium(double radius, double tempin, std::array<char,CM> &constmodels, 
    const threevector& position, const threevector& velocity)
    :Matter(radius,tempin,LithiumConsts){
    E_Debug("\n\nIn Lithium::Lithium(double radius, double tempin, ...)\n\n\t");
    set_defaults();

    update_state(0.0);
    update_models(constmodels);
    update_motion(position,velocity,0.0);
    update();
    E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius 
        << "\nSt.Density = " << St.Density << "\nSt.Volume = " << St.Volume);
}

void Lithium::set_defaults(){
    E_Debug("\tIn Lithium::set_defaults()\n\n");

    //!< https://www.engineersedge.com/
    //!< materials/specific_heat_capacity_of_metals_13259.htm
    St.HeatCapacity = 4.169;              //!< kJ/(kg-K)         
    St.Emissivity = 0.1;                  //!< http://www.fusion.ucla.edu/APEX/
                                          //!< meeting15/Apex4_01-tanaka.pdf
    St.SuperBoilingTemp = Ec.BoilingTemp; //!< K, At any pressure
    St.Density = 534;                     //!< kg/m^3, from Wikipedia
    update_models('c','c','c','y','n');
    
//  St.ThermConduct = 0.163;        // kW/m K at 20 degrees celsius
}

void Lithium::update_heatcapacity(){
    E_Debug("\tIn Lithium::update_heatcapacity()\n\n");
    //!< Temperature dependant heat capacity for Lithium taken from:
    //!< D. Harry W., Lewis Reserch Cent. 24 (1968), pg 8, figure 4
    if( St.Temperature < Ec.MeltingTemp ){ 
        static bool runOnce = true;
        std::string WarningMessage = "In Lithium::update_heatcapacity():\n";
        WarningMessage += "Extending heat capacity model outside temperature";
        WarningMessage += " range! T < 300K";
        WarnOnce(runOnce,WarningMessage);

        St.HeatCapacity = 4169-0.2427*St.Temperature+
            1.045e-3*pow(St.Temperature,2.0);
    }else if( St.Temperature >= Ec.MeltingTemp 
        && St.Temperature <= Ec.BoilingTemp ){ 
        St.HeatCapacity = 4169-0.2427*St.Temperature+
            1.045e-3*pow(St.Temperature,2.0);
    }else{
        static bool runOnce = true;
        std::string WarningMessage = "In Lithium::update_heatcapacity():\n";
        WarningMessage += "Extending heat capacity model outside temperature";
        WarningMessage += " T > Ec.BoilingTemp";
        WarnOnce(runOnce,WarningMessage);

        St.HeatCapacity = 4169-0.2427*St.Temperature+
        1.045e-3*pow(St.Temperature,2.0);
    }

    E1_Debug("\n\nTemperature is : " << St.Temperature << "\nSt.Gas = " 
        << St.Gas << "\nSt.Liquid = " << St.Liquid << "\nCv of Solid: " 
        << St.HeatCapacity/Ec.AtomicMass << "[kJ/(kg K)]"; );
    //!< Conversion J/(Kg K) to kJ/( kg K ),
    St.HeatCapacity = (St.HeatCapacity /1000); 
}

void Lithium::update_radius(){
    E_Debug("\tIn Lithium::update_radius()\n\n");
    
    double DensityTemp = 562.0-0.1*St.Temperature;
    St.LinearExpansion = 1.0+pow(534/DensityTemp,1.0/3.0);
    St.Radius=St.UnheatedRadius*St.LinearExpansion;
    E1_Debug("\nTemperature = " << St.Temperature 
        << "\n\nSt.LinearExpansion = " << St.LinearExpansion << "\nSt.Radius = "
        << St.Radius);
    assert(St.Radius>0); // Assert radius is positive   
}

void Lithium::update_vapourpressure(){
    E_Debug("\tIn Lithium::update_vapourpressure()\n\n");
    St.VapourPressure = probe_vapourpressure(St.Temperature);
}

double Lithium::probe_vapourpressure(double Temperature)const{
    double VapourPressure(0.0);
    //!< Model being used:
    //!< http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
    //!< High temperature model:
    //!< https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680018893.pdf
    if( Temperature < Ec.MeltingTemp ){
        static bool runOnce = true;
        std::string WarningMessage = "In Lithium::update_vapourpressure():\n";
        WarningMessage += "Extending model outside temperature range!";
        WarningMessage += " (from 298K< to 0K<)";
        WarnOnce(runOnce,WarningMessage);

        VapourPressure = 101325*pow(10,5.667 - 8310/Temperature); 
    }else if( Temperature >= Ec.MeltingTemp && Temperature < 800 ){
        VapourPressure = 101325*pow(10,5.055 - 8023/Temperature); 
    }else if( Temperature >= 800 ){
        
        VapourPressure = 101325*pow(10,10.015-8064.5/Temperature);
    }else{
        std::cerr << "\nError! Negative Temperature in Lithium::" 
            << "probe_vapourpressure(double "<< Temperature <<")";
    }
    
    return VapourPressure;
}