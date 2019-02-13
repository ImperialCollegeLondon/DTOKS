#include "Functions.h"
#include "Constants.h"

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/
// Thomas-D-2016-PhD-thesis.pdf
double solveSMOML(double TemperatureRatio, double MassRatio, double Ionization,
double HeatCapacityRatio, double PlasmaFlowSpeed){
	PlasmaFlowSpeed = PlasmaFlowSpeed/sqrt(2);
	double s1 = Ionization*sqrt(PI)*(1+2*pow(PlasmaFlowSpeed,2))*
        erf(PlasmaFlowSpeed)/(4*PlasmaFlowSpeed)+exp(-pow(PlasmaFlowSpeed,2))/2;
	double s2 = Ionization*Ionization*sqrt(PI)*erf(PlasmaFlowSpeed)/
        (2*PlasmaFlowSpeed);
	
	return TemperatureRatio*s1/s2-LambertW(sqrt(2*PI*TemperatureRatio*
        (1+HeatCapacityRatio*TemperatureRatio))*exp(TemperatureRatio*s1/s2)/s2)+
        log(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio)/2;
}

