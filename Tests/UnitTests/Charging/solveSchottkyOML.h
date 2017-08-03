#include "LambertW.h"
#include "Constants.h"

// Solve the Modified orbital motion limited potential for large dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double solveSchottkyOML(double DustTemperature, double Ti, double Te){


	double Wf = 3.4*echarge;
	double n0 = 10e19;
	double TemperatureRatio = Ti/Te;
	double Z = 1.0; 		// Ionization
	double vi = sqrt(Kb*Ti/(2*PI*Mp));
	double ve = sqrt(Kb*Te/(2*PI*Me));

	double Arg = TemperatureRatio*exp(TemperatureRatio/Z)*(ve-Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature))/(Z*echarge*n0))/(Z*vi);
	double Coeff = TemperatureRatio/Z;
//	std::cout << "\nArg = " << Arg;
	if( Arg < 0 ){
		// Swap vi and ve over and change sign of electron emission...
		Arg = TemperatureRatio*exp(TemperatureRatio)*(vi+Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature))/(Z*echarge*n0))/(ve);
		Coeff = TemperatureRatio;
		return Coeff + LambertW(Arg);
	}
	return Coeff - LambertW(Arg);
}

