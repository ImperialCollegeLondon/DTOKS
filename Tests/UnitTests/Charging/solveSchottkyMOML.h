#include "LambertW.h"
#include "Constants.h"

// Solve the Modified orbital motion limited potential for large dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double solveSchottkyMOML(double DustTemperature, double Te, double Ti, double HeatCapacityRatio, double MassRatio){


	double Wf = 3.4*echarge;
	double n0 = 10e19;
	double TemperatureRatio = Ti/Te;
	double cs = sqrt(echarge*(Te+HeatCapacityRatio*Ti)/Mp);
//f(DustTemperature > 1000 ){
//		std::cout << "\nRichardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature)))/(n0*echarge) = " << Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature))/(n0*echarge);
//		std::cout << "\nexp(-Wf/(Kb*DustTemperature)) = " << exp(-Wf/(Kb*DustTemperature)) << "\n-Wf/(Kb*DustTemperature) = " << -Wf/(Kb*DustTemperature);
//		std::cout << "\nArg = " << sqrt(MassRatio*TemperatureRatio)*exp(-TemperatureRatio)*(1/sqrt(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio))-(Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature)))/(n0*echarge);

//		std::cin.get();
//	
	double Arg =  sqrt(MassRatio*TemperatureRatio)*exp(-TemperatureRatio)*(1/sqrt(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio))-(Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature)))/(n0*echarge*cs);
	if( Arg < 0.0 )	return 0.0;
	else	return TemperatureRatio-LambertW(sqrt(MassRatio*TemperatureRatio)*exp(-TemperatureRatio)*(1/sqrt(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio))-(Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature)))/(n0*echarge*cs))-log((1/sqrt(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio))-(Richardson*pow(DustTemperature,2)*exp(-Wf/(Kb*DustTemperature)))/(n0*echarge*cs));
}

