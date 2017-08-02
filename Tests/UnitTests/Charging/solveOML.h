#include "LambertW.h"
#include "Constants.h"

// Solve the orbital motion limited potential for large dust grains.
// See drews Thesis, pg 48
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double solveOML(double TemperatureRatio, double Ionization, double MassRatio){

	return TemperatureRatio/Ionization-LambertW(sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/Ionization)/Ionization);
}

