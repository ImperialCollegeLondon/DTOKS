#include "Functions.h"
#include "Constants.h"

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double solveSOML(double TemperatureRatio, double PlasmaFlowSpeed, double MassRatio){
	PlasmaFlowSpeed = PlasmaFlowSpeed/sqrt(2);
	double s1 = sqrt(PI)*(1+2*pow(PlasmaFlowSpeed,2))*erf(PlasmaFlowSpeed)/(4*PlasmaFlowSpeed)+exp(-pow(PlasmaFlowSpeed,2))/2;
	double s2 = sqrt(PI)*erf(PlasmaFlowSpeed)/(2*PlasmaFlowSpeed);

	return TemperatureRatio*s1/s2-LambertW(sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio*s1/s2)/s2);
}

