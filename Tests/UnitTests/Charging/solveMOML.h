#include "LambertW.h"
#include "Constants.h"

// Solve the Modified orbital motion limited potential for large dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double solveMOML(double TemperatureRatio, double HeatCapacityRatio, double MassRatio){

	return TemperatureRatio-LambertW(sqrt(2*PI*TemperatureRatio*(1+HeatCapacityRatio*TemperatureRatio))*exp(TemperatureRatio))+log(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio)/2;
}

