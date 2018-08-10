#include "Constants.h"
#include "Functions.h"
#include <math.h>

// Solve the Modified orbital motion limited potential for large emitting dust grains.
// See the paper by Minas and Nikoleta, equation (1) and (2)
// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
double solveMOMLEM(double TemperatureRatio, double MassRatio, double EmittedDensityRatio,
			double EmittedTemperatureRatio, double guess){

	double a = 1.0/EmittedDensityRatio;
	double b = 1.0/TemperatureRatio;
	double c = log(sqrt(TemperatureRatio/MassRatio)+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
	double d = EmittedTemperatureRatio;
	// Uncomment following line to compare MOMLEM results with figure (1) of paper:
	// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
	

//	std::cout << "\na = " << a;
//	std::cout << "\nb = " << b;
//	std::cout << "\nc = " << c;
//	std::cout << "\nd = " << d;

	double x1 = guess;
	do{

		guess = x1;

		double fx = a*exp(-guess/b)*(1+erf(-sqrt(-guess/b)))-exp(guess)*(1+erf(sqrt(guess-c)))*a
				-exp((guess-c)/d)*(1-erf(sqrt((guess-c)/d)));
		double fxprime = -b*a*exp(-guess/b)*(1+erf(-sqrt(-guess/b)))+a/(sqrt(PI)*b*sqrt(-guess/b))
				-exp(guess)*(1+erf(sqrt(guess-c)))+exp(c)/(sqrt(PI)*sqrt(guess-c))
				-(1.0/d)*exp((guess-c)/d)*(1-erf(sqrt((guess-c)/d)))
				+1.0/(sqrt(PI)*d*sqrt((guess-c)/d));
		x1 = guess - fx/fxprime;
//		std::cout << "\nguess = " << guess;
//		std::cout << "\nx1 = " << x1;
//		std::cout << "\nfx = " << fx;
//		std::cout << "\nfxprime = " << fxprime; std::cin.get();

	}while(fabs(guess-x1)>1e-2);// >1e-2){


	return guess;
}

