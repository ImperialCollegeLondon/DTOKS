/** @file Tungsten.cpp
 *  @brief Implementation of Tungsten class with element constants
 *
 *  Constructors for beryllium class and element specific dependencies of heat
 *  capacity, vapour pressure and size on temperature.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#include "Tungsten.h"
#include "Constants.h"

/** @brief Constant structure defining true physical constants of Tungsten
 */
const struct ElementConsts TungstenConsts = {
    'W',          //!< Specifies the element
    3422,         //!< K, Melting temperature at atmospheric pressure
    5555,         //!< K, Boiling temperature at atmospheric pressure
    774,          //!< kJ/mol, Bond Energy, 
                  //!< Estimated as being equal to Latent Vapour Energy
    3.4,          //!< ev, Work Function, 
                  //!< Taken from DTOKS and matched with wikipedia
    5.7,          //!< kJ/(m^2 K), HeatTransfer coefficient, 
                  //!< http://www.engineeringtoolbox.com/
                  //!< overall-heat-transfer-coefficients-d_284.html
    0.18384,      //!< AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
    35.3/0.18384, //!< kJ/kg, Latent Fusion Energy From Wikipedia
    774/0.18384,  //!< kJ/kg, Latent Vapour Energy From Wikipedia
    2.333,        //!< N/m, Surface Tension B. Keene, 1993, 
                  //!< Review of data for the surface tension of pure metals
    19211,        //!< kg/m^3, at 298.15K from: P. Hidnert and W.T. Sweeney, 
                  //!< Sci. Pap. Bur. Stand. 20, 483 (1925). 
    0.163         //!< kW/m K, at 20 degrees celsius
};

Tungsten::Tungsten():
Matter(TungstenConsts){
    E_Debug("\n\nIn Tungsten::Tungsten():Matter(&TungstenConsts)\n\n");
    set_defaults();
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Tungsten::Tungsten(double radius):
Matter(radius,TungstenConsts){
    E_Debug("\n\nIn Tungsten::Tungsten(double radius):"
        << "Matter(radius,TungstenConsts)\n\n");
    set_defaults();
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Tungsten::Tungsten(double radius, double tempin):
Matter(radius,tempin,TungstenConsts){
    E_Debug("\n\nIn Tungsten::Tungsten(double radius, double tempin):"
        << "Matter(radius,tempin,TungstenConsts)\n\n");
    set_defaults();
    update_state(0.0);
    update_models('c','c','c','y','n');
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Tungsten::Tungsten(double radius, double tempin, 
std::array<char,CM> &constmodels):
Matter(radius,tempin,TungstenConsts){
    E_Debug("\n\nIn Tungsten::Tungsten(double radius, double tempin, "
        << "std::array<char,CM> &constmodels):"
        << "Matter(radius,tempin,TungstenConsts))\n\n");
    set_defaults();

    update_state(0.0);
    update_models(constmodels);
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Tungsten::Tungsten(double radius, double tempin, 
std::array<char,CM> &constmodels, const threevector& position, 
const threevector& velocity):
Matter(radius,tempin,TungstenConsts){
    E_Debug("\n\nIn Tungsten::Tungsten(double radius, double tempin, "
        << "std::array<char,CM> &constmodels, const threevector& position, "
        << "const threevector& velocity):"
        << "Matter(radius,tempin,TungstenConsts))\n\n");
    set_defaults();

    update_state(0.0);
    update_models(constmodels);
    update_motion(position,velocity,0.0);
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

void Tungsten::set_defaults(){
    E_Debug("\n\n\tIn Tungsten::set_defaults()");

    //!< http://www.engineersedge.com/
    //!< materials/specific_heat_capacity_of_metals_13259.html
    //!< http://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
    St.HeatCapacity = 0.13398;            //!< kJ/(kg-K)         
    St.Emissivity = 0.04;                 //!< Arb, 
    St.SuperBoilingTemp = Ec.BoilingTemp; //!< K, At any pressure
    St.Density = 19600;                   //!< (kg/m^3) from Wikipedia
    E_Debug("\n\n\t"); update_models('c','c','c','y','n');
    
//  St.ThermConduct = 0.163;        // kW/m K at 20 degrees celsius
}

void Tungsten::update_heatcapacity(){
    E_Debug("\n\n\tIn Tungsten::update_heatcapacity()");

    if( St.Temperature < 300 ){ 

        static bool runOnce = true;
        std::string WarningMessage = "In Tungsten::update_heatcapacity():\n";
        WarningMessage += "Extending heat capacity model outside temperature";
        WarningMessage += " range! T < 300K";
        WarnOnce(runOnce,WarningMessage);

        St.HeatCapacity = 24.943 - 7.72e4*pow(St.Temperature,-2) + 
            2.33e-3*St.Temperature + 1.18e-13*pow(St.Temperature,4);
    }else if( St.Temperature > 300 && St.Temperature < Ec.MeltingTemp ){ 
        //!< http://nvlpubs.nist.gov/nistpubs/jres/75a/jresv75an4p283_a1b.pdf
        St.HeatCapacity = 24.943 - 7.72e4*pow(St.Temperature,-2) +
            2.33e-3*St.Temperature + 1.18e-13*pow(St.Temperature,4);
    }else{
        //!< http://webbook.nist.gov/cgi/inchi?ID=C7440337&Mask=2

        double t = St.Temperature/1000;
                St.HeatCapacity = 35.56404 -1.551741e-7*t + 2.915253e-8*pow(t,2)
                -1.891725e-9*pow(t,3)-4.107702e-7*pow(t,-2);
    }

    E1_Debug("\n\tTemperature is : " << St.Temperature << "\n\tSt.Gas = " 
        << St.Gas << "\n\tSt.Liquid = " << St.Liquid << "\n\tCv of Solid: " 
        << St.HeatCapacity/Ec.AtomicMass << "[kJ/(kg K)]"; );
    //!< Conversion kJ/(mol K) to kJ/( kg K ), AtomicMass [kg mol^-1]
    St.HeatCapacity = (St.HeatCapacity /(1000 * Ec.AtomicMass)); 
}

void Tungsten::update_radius(){
    E_Debug("\n\n\tIn Tungsten::update_radius()");
    if( St.Temperature > 173 && St.Temperature <= 1500 ){
        static bool runOnce = true;
        std::string WarningMessage = "In Tungsten::update_radius():\n";
        WarningMessage += "Extending model outside temperature range!";
        WarningMessage += " (from 738K to 1500K)";
        WarnOnce(runOnce,WarningMessage);
        St.LinearExpansion = 1+(4.28*St.Temperature)*1e-6;
    }else if( St.Temperature > 1500 && St.Temperature < Ec.MeltingTemp ){
        St.LinearExpansion = 1+3.9003e-4+
            1.3896*1e-3-8.2797*1e-7*St.Temperature +
            4.0557*1e-9*pow(St.Temperature,2)-
            1.2164*1e-12*pow(St.Temperature,3)+
            1.7034*1e-16*pow(St.Temperature,4);

    }else if( St.Temperature == Ec.MeltingTemp ){ //!< Model while it's melting
        //!< http://nvlpubs.nist.gov/nistpubs/jres/75a/jresv75an4p283_a1b.pdf
        St.LinearExpansion = 1.02105 +
            0.03152*St.FusionEnergy/(Ec.LatentFusion*St.Mass);

    }else if( St.Temperature > Ec.MeltingTemp 
        && St.Temperature <= St.SuperBoilingTemp){
        St.LinearExpansion = pow(1.18+6.20*1e-5*(St.Temperature-3680)+
            3.23*1e-8*pow((St.Temperature-3680),2),(1./3.));
    }
    St.Radius=St.UnheatedRadius*St.LinearExpansion;
    E1_Debug("\n\tTemperature = " << St.Temperature 
        << "\n\tSt.LinearExpansion = "<< St.LinearExpansion 
        << "\n\tSt.Radius = " << St.Radius);
    assert(St.Radius>0);  
}

void Tungsten::update_vapourpressure(){
    E_Debug("\n\n\tIn Tungsten::update_vapourpressure()");
    St.VapourPressure = probe_vapourpressure(St.Temperature);
}

double Tungsten::probe_vapourpressure(double Temperature)const{
    E_Debug("\n\n\tIn Tungsten::probe_vapourpressure(double Temperature)"
        << "const");
    double VapourPressure(0.0);
    //!< Vapour pressure model:
    //!< // http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
    if( Temperature < 298 ){
        static bool runOnce = true;
        std::string WarningMessage = "In Tungsten::update_vapourpressure():\n";
        WarningMessage += "Extending model outside temperature range!";
        WarningMessage += " (from 298K< to 0K<)";
        WarnOnce(runOnce,WarningMessage);

        VapourPressure = 101325*pow(10,2.945 - 44094/Temperature +
            1.3677*log10(Temperature)); 
    }else if( Temperature >= 298 && Temperature < 2500 ){
        VapourPressure = 101325*pow(10,2.945 - 44094/Temperature +
        1.3677*log10(Temperature));
    }else if( Temperature >= 2500 ){
        //!< E. R. Plante and A. B. Sessoms
        VapourPressure = 101325*pow(10,7.871-45385/Temperature);
    }else{
        std::cerr << "\nError! Negative Temperature in" 
            << "Tungsten::probe_vapourpressure(double Temperature)";
    }
    
    return VapourPressure;
}
