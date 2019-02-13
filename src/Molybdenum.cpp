/** @file Molybdenum.cpp
 *  @brief Implementation of Molybdenum class with element constants
 *
 *  Constructors for beryllium class and element specific dependencies of heat
 *  capacity, vapour pressure and size on temperature.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#include "Molybdenum.h"
#include "Constants.h"

/** @brief Constant structure defining true physical constants of Molybdenum
 */
const struct ElementConsts MolybdenumConsts = {
    'M',     //!< Specifies the element
    2896.0,  //!< K, Melting temperature at atmospheric pressure
    4912.0,  //!< K, Boiling temperature at atmospheric pressure
    681.0,   //!< kJ/mol, Bond Energy, 
             //!< Estimated as being equal to Latent Vapour Energy
    4.2,     //!< ev, Work Function, Taken from DTOKS and matched with wikipedia
    1.0,     //!< kJ/(m^2 K), HeatTransfer coefficient, 
             //!< http://www.engineeringtoolbox.com
             //!< /overall-heat-transfer-coefficients-d_284.html
    0.09595, //!< AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
    390.62,  //!< kJ/kg, Latent Fusion Energy From Wikipedia
    6232.41, //!< kJ/kg, Latent Vapour Energy From Wikipedia
    2.24,    //!< N/m, Surface Tension B. Keene, 1993, 
             //!< Review of data for the surface tension of pure metals
    10280,   //!< (kg/m^3) from Wikipedia, Room temperature density
    0.138    //!< kW/m K at 20 degrees celsius
};

Molybdenum::Molybdenum():
Matter(MolybdenumConsts){
    E_Debug("\n\nIn Molybdenum::Molybdenum():Matter(MolybdenumConsts)\n\n");
    set_defaults();
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Molybdenum::Molybdenum(double radius):
Matter(radius,MolybdenumConsts){
    E_Debug("\n\nIn Molybdenum::Molybdenum(double radius):"
        << "Matter(radius,MolybdenumConsts)\n\n");
    set_defaults();
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Molybdenum::Molybdenum(double radius, double tempin):
Matter(radius,tempin,MolybdenumConsts){
    E_Debug("\n\nIn Molybdenum::Molybdenum(double radius, double tempin):"
        << "Matter(radius,tempin,MolybdenumConsts)\n\n");
    set_defaults();
    update_state(0.0);
    update_models('c','c','c','y','n');
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Molybdenum::Molybdenum(double radius, double tempin, 
std::array<char,CM> &constmodels):
Matter(radius,tempin,MolybdenumConsts){
    E_Debug("\n\nIn Molybdenum::Molybdenum(double radius, double tempin, "
        << "std::array<char,CM> &constmodels):"
        << "Matter(radius,tempin,MolybdenumConsts)\n\n");
    set_defaults();

    update_state(0.0);
    update_models(constmodels);
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

Molybdenum::Molybdenum(double radius, double tempin, 
std::array<char,CM> &constmodels, const threevector& position, 
const threevector& velocity):
Matter(radius,tempin,MolybdenumConsts){
    E_Debug("\n\nIn Molybdenum::Molybdenum(double radius, double tempin, "
        << "std::array<char,CM> &constmodels, const threevector& position, "
        << "const threevector& velocity):"
        << "Matter(radius,tempin,MolybdenumConsts))\n\n");
    set_defaults();

    update_state(0.0);
    update_models(constmodels);
    update_motion(position,velocity,0.0);
    update();
    E1_Debug("\n\tMass after = " << St.Mass << "\n\tRadius After = " 
        << St.Radius << "\n\tSt.Density = " << St.Density << "\n\tSt.Volume = "
        << St.Volume);
}

void Molybdenum::set_defaults(){
    E_Debug("\n\n\tIn Molybdenum::set_defaults()");

    //!< http://www.engineersedge.com/
    //!< materials/specific_heat_capacity_of_metals_13259.html
    St.HeatCapacity = 0.27716616;         //!< kJ/(kg-K)         
    St.Emissivity = 0.1;                  //!< Arb, 
                                          //!< http://www.engineeringtoolbox.com
                                          //!< /emissivity-coefficients-d_447
                                          //!< .html
    St.SuperBoilingTemp = Ec.BoilingTemp; //!< K, At any pressure
    St.Density = 10280;                   //!< kg/m^3, from Wikipedia
    update_models('c','c','c','y','n');
    

//  St.ThermConduct = 0.163;        // kW/m K at 20 degrees celsius
}

void Molybdenum::update_heatcapacity(){
    E_Debug("\n\n\tIn Molybdenum::update_heatcapacity()");

    double Temperatures[104] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 
        18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 
        110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 
        250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 650, 700, 
        750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 
        2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500, 2550, 
        2600, 2650, 2700, 2750, 2800, 2850, 2896};

    double HeatCapacities[104] = {0.0021, 0.0047, 0.0071, 0.0092, 0.0132, 
        0.0172, 0.0219, 0.0275, 0.0333, 0.0444, 0.067, 0.092, 0.108, 0.129, 
        0.188, 0.255, 0.494, 0.891, 1.43, 2.14, 2.89, 3.76, 4.69, 5.72, 6.81, 
        7.89, 9.03, 9.96, 10.93, 11.83, 12.62, 13.45, 14.88, 16.08, 17.11, 
        17.98, 18.78, 19.45, 20.03, 20.56, 21.01, 21.44, 21.81, 22.15, 22.48, 
        22.76, 22.99, 23.21, 23.40, 23.59, 23.76, 23.92, 24.59, 25.10, 25.49, 
        25.83, 26.14, 26.41, 26.67, 26.89, 27.14, 27.36, 27.58, 27.81, 28.03, 
        28.27, 28.51, 28.78, 29.05, 29.36, 29.68, 30.04, 30.38, 30.75, 31.13, 
        31.53, 31.96, 32.39, 32.83, 33.28, 33.76, 34.29, 34.77, 35.31, 35.88, 
        36.43, 37.09, 37.71, 38.36, 39.03, 39.74, 40.45, 41.24, 42.03, 42.86, 
        43.79, 44.72, 45.72, 46.82, 47.97, 49.15, 50.49, 51.86, 53.29};
    double MinDiff(1000.0);
    int MinIndex(0);
    for(int i(0); i < 104; i ++){
        double Diff = Temperatures[i] - St.Temperature;
        if( Diff < MinDiff ){
            MinDiff = Diff;
            MinIndex = i;
        }
    }

    //<! Conversion J/(mol K) to kJ/( kg K ), AtomicMass [kg mol^-1]
    St.HeatCapacity = HeatCapacities[MinIndex]/(1000 * Ec.AtomicMass); 
}

void Molybdenum::update_radius(){
    E_Debug("\n\n\tIn Molybdenum::update_radius()");
    
    double Temperatures[104] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 16, 
        18, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
        110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 
        250, 260, 270, 280, 290, 300, 350, 400, 450, 500, 550, 600, 650, 700, 
        750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950,
        2000, 2050, 2100, 2150, 2200, 2250, 2300, 2350, 2400, 2450, 2500, 2550,
        2600, 2650, 2700, 2750, 2800, 2850, 2896};

    double Extensions[104] = {0.0026, 0.0053, 0.0083, 0.0115, 0.0152, 0.0195, 
        0.0244, 0.0301, 0.0366, 0.0442, 0.063, 0.087, 0.101, 0.116, 0.153, 
        0.198, 0.351, 0.572, 0.876, 1.28, 1.79, 2.42, 3.11, 3.75, 4.42, 5.06, 
        5.63, 6.24, 6.79, 7.34, 7.81, 8.35, 9.26, 10.02, 10.70, 11.33, 11.77, 
        12.27, 12.71, 13.11, 13.41, 13.73, 13.97, 14.19, 14.40, 14.55, 14.71, 
        14.81, 14.94, 15.04, 15.14, 15.21, 15.61, 15.93, 16.21, 16.41, 16.62, 
        16.81, 17.02, 17.24, 17.48, 17.65, 17.94, 18.16, 18.40, 18.64, 18.97, 
        19.25, 19.55, 19.91, 20.35, 20.72, 21.15, 21.60, 21.99, 22.44, 22.86, 
        23.34, 23.82, 24.35, 24.84, 25.38, 25.86, 26.35, 26.96, 27.61, 28.24, 
        28.94, 29.71, 30.51, 31.38, 32.3, 33.18, 34.02, 34.95, 35.98, 36.96, 
        37.99, 39.12, 40.32, 41.64, 42.94, 44.38, 45.80};

    double MinDiff(1000.0);
    int MinIndex(0);
    for(int i(0); i < 104; i ++){
        double Diff = Temperatures[i] - St.Temperature;
        if( Diff < MinDiff ){
            MinDiff = Diff;
            MinIndex = i;
        }
    }

    St.LinearExpansion = pow(Extensions[MinIndex]*1e-6,1.0/3.0);
    St.Radius=St.UnheatedRadius*St.LinearExpansion;

    E1_Debug("\n\tTemperature = " << St.Temperature 
        << "\n\tSt.LinearExpansion = " << St.LinearExpansion 
        << "\n\tSt.Radius = " << St.Radius);
    assert(St.Radius>0);
}

void Molybdenum::update_vapourpressure(){
    E_Debug("\n\n\tIn Molybdenum::update_vapourpressure()");
    St.VapourPressure = probe_vapourpressure(St.Temperature);
}

double Molybdenum::probe_vapourpressure(double Temperature)const{
    E_Debug("\n\n\tIn Molybdenum::probe_vapourpressure(double Temperature)"
        << "const");
    double VapourPressure(0.0);
    
    if( Temperature < Ec.MeltingTemp ){
        VapourPressure = 101325*pow(10,11.529 -34626/Temperature -
            1.1331*log10(Temperature));
    }else{
        VapourPressure = 101325*pow(10,11.529 -34626/Temperature -
            1.1331*log10(Temperature));
        static bool runOnce = true;
        std::string WarningMessage = "In Molybdenum::probe_vapourpressure():\n";
        WarningMessage += "Extending model outside temperature range!";
        WarningMessage += " ( from St.MeltingTemp > to St.BoilingTemp > )";
        WarnOnce(runOnce,WarningMessage);
        std::cerr << "\nError! Negative Temperature in "
            << "Molybdenum::probe_vapourpressure(double Temperature)";
    }
    return VapourPressure;
}