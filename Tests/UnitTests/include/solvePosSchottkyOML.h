#include "Functions.h"
#include "Constants.h"

// Solve the orbital motion limited potential for small dust grains accounting for electron emission.
// WARNING: This is only valid for positively charged dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double solvePosSchottkyOML(double DustTemperature, double Ti, double Te, double n0, double Dsec){

	double Wf = 3.4*echarge;
	double TemperatureRatio = Ti/Te;
	double Z = 1.0; 		// Ionization
	double vi = sqrt(Kb*Ti/(2*PI*Mp));
	double ve = sqrt(Kb*Te/(2*PI*Me));

	double Arg = TemperatureRatio*exp(TemperatureRatio)*(vi+Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature))/(Z*echarge*n0))/(ve*(1-Dsec));
	double Coeff = TemperatureRatio;
	return -Coeff + LambertW(Arg);
}
