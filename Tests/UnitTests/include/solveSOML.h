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

/*
double solvePosSOML(double guess, double TemperatureRatio, double PlasmaFlowSpeed, double MassRatio){
	PlasmaFlowSpeed = PlasmaFlowSpeed/sqrt(2);
	double Z = 1.0;
	double uip = PlasmaFlowSpeed+sqrt(fabs(Potential));
	double uim = PlasmaFlowSpeed-sqrt(fabs(Potential));

	// Ion Current Term 1 and Term 2
	double Ti1 = (PlasmaFlowSpeed+(1.0/(2.0*PlasmaFlowSpeed))*(1-2.0*Z*Potential/TemperatureRatio))
			*(erf(uip)+erf(uim));
	double Ti2 = (1.0/(PlasmaFlowSpeed*sqrt(PI)))*(uip*exp(-uim*uim)+uim*exp(-uip*uip));
	double Ecurr = (1.0/sqrt(Z*Z*TemperatureRatio*MassRatio))*(1.0+Potential);

	double x1 = guess;
	do{
		guess = x1;

		double fx = Ti1+Ti2-Ecurr;

		// Attempted to differentiate Ion Current Term 1 and Term 2... This becomes very complex
		double fxprime = (1.0/sqrt(PI*x1))*exp(-pow(PlasmaFlowSpeed+sqrt(x1),2.0))*exp(-uim*uim)+ 
				(1.0/sqrt(PI*x1))*exp(-pow(PlasmaFlowSpeed-sqrt(x1),2.0))*exp(-uip*uip)+
				
				-(1.0/sqrt(Z*Z*TemperatureRatio*MassRatio)); // d Ecurr / dPotential
		x1 = guess - fx/fxprime;

	}while(fabs(guess-x1)>1e-4);// >1e-2){

	return TemperatureRatio*s1/s2-LambertW(sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio*s1/s2)/s2);
}
*/
