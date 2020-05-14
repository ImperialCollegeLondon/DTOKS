#include "Constants.h"
#include "Functions.h"

// Solve the Modified orbital motion limited potential for large dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/
// Thomas-D-2016-PhD-thesis.pdf
double solveMOML(double TemperatureRatio, double MassRatio, double Ionization,
    double HeatCapacityRatio ){

    return TemperatureRatio/Ionization-
        LambertW(sqrt(2*PI*TemperatureRatio*
        (1+HeatCapacityRatio*TemperatureRatio))*
        exp(TemperatureRatio/Ionization))/Ionization+
        log(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio)/2;
}

double MOMLattractedCurrent(double Ti, double Te, double HeatCapacityRatio, 
    double MassRatio, double Potential){

    double SheathPotential = -Potential*echarge*Te/echarge - (1.0/2.0)*
        log((2.0*PI/MassRatio)*(1.0+HeatCapacityRatio*Ti/Te));

    return 1-echarge*SheathPotential/(echarge*Ti);
}

