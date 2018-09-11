#include "Constants.h"
#include "Functions.h"
#include "math.h"

// BIBHAS Charging Test:
// This test is designed to find the floating potential for arbitary sized dust grain
// The calculation follows the work by R. DE Bibhas,
// see R. DE Bibhas, Astrophys. Space Sci. 30, (1974).
double solveBIBHAS(double Ti, double Te, double MassRatio, double Ionization, double DebyeLength, double guess){ 
	double s = 1.0+DebyeLength; 
	double Radius = 1.0; // Normalised
	double x1 = guess;
	double HeatCapacityRatio = 1.0;

//	double r = Radius*(1.0-0.17/sqrt(Ionization))*sqrt(1.0/(2.0*Ionization))*(1.0+sqrt(1.0+2.0*Ionization));
//	double Pot_m = (Ionization*echarge*echarge/r - echarge*echarge*Radius/(2.0*(r*r-Radius*Radius))
//			+ehcarge*echarge*Radius/(2*r*r))/(echarge*Te);

	double Coeff = Ionization*sqrt(Ti/(Te*MassRatio))*(s*s/(Radius*Radius));
	double a = Coeff;
	double b = Coeff*((s*s-1.0)/(s*s));
	do{
		guess = x1;
		double SheathPotential = guess-(1.0/2.0)*log((1.0+HeatCapacityRatio*(Ti/Te))*(2.0*PI)/MassRatio);
		double fx = a-b*exp(Ti/(Te*(s*s-1.0*1.0)))-exp(SheathPotential);
		double fxprime = -(b/(s*s-1.0*1.0))*exp(Ti/(Te*(s*s-1.0*1.0)))-exp(SheathPotential);
		x1 = guess - fx/fxprime;

//		double fx = a-b*exp(guess/(s*s-1.0*1.0))-exp(guess);
//		double fxprime = -(b/(s*s-1.0*1.0))*exp(guess*Ti/(Te*(s*s-1.0*1.0)))-exp(guess);
//		x1 = guess - fx/fxprime;

	}while(fabs(guess-x1)>1e-4);// >1e-2){

/*
	if( Ti != Te ){
		std::cerr << "\nWarning! Mopel is invalid for Ti != Te";
	}
	do{
		guess = x1;
		double u_d = sqrt(echarge*Ti/(2.0*PI*MassRatio*Me))*exp(-guess)*(s*s/(Radius*Radius))
			*(1.0-((s*s-Radius*Radius)/(s*s))*exp(guess*Radius*Radius*Ti/(Te*(s*s-Radius*Radius))));
		double u_dprime = -1.0*u_d-sqrt(echarge*Ti/(2.0*PI*MassRatio*Me))
			*exp(-guess)*exp(guess*Radius*Radius*Ti/(Te*(s*s-Radius*Radius)));

//		std::cout << "\nu_dprime = " << u_dprime;
//		std::cout << "\nexp(-x) = " << exp(-guess);
//		std::cout << "\nexp(x*a^2/(s^2-a^2)) = " << exp(guess*Radius*Radius*Ti/(Te*(s*s-Radius*Radius)));
//		std::cout << "\nU = " << sqrt(Kb*Ti/(2.0*PI*MassRatio*Me));
		
		double p = sqrt(MassRatio*Me/(2.0*echarge*Ti))*u_d;
		double p_prime = sqrt(MassRatio*Me/(2.0*echarge*Ti))*u_dprime;
		double q = sqrt((s*s/(s*s-Radius*Radius))*MassRatio*Me/(2.0*echarge*Ti))*u_d;
		double q_prime = sqrt((s*s/(s*s-Radius*Radius))*MassRatio*Me/(2.0*echarge*Ti))*u_dprime;

//		std::cout << "\nsqrt(MassRatio*Me/(2.0*echarge*Ti)) = " << sqrt(MassRatio*Me/(2.0*echarge*Ti));
//		std::cout << "\nu_d = " << u_d;
//		std::cout << "\nu_d_lim = " << sqrt(echarge*Ti/(2.0*PI*MassRatio*Me))*exp(-guess)*(1.0-guess);
//		std::cout << "\n\nPot = " << guess;
//		std::cout << "\nPot_lim = " << 1.11;
//		std::cout << "\n\nP = " << p;
//		std::cout << "\nQ = " << q;
//		std::cout << "\n\nP' = " << p_prime;
//		std::cout << "\nQ' = " << q_prime;

		double Coeff = Ionization*sqrt(Ti/(Te*MassRatio))*(s*s/(Radius*Radius));
		double a = Coeff*(exp(-p*p)+sqrt(PI)*p*(1.0+erf(p)));
		double b = Coeff*(((s*s-Radius*Radius)/(s*s))*exp(Radius*Radius*guess*Ti/(Te*(s*s-Radius*Radius)))
				*(exp(-q*q)+sqrt(PI)*q*(1.0+erf(q))));

		double a_prime = Coeff*(sqrt(PI)*p_prime*(1.0+erf(p)));
		double b_prime = Coeff*(((s*s-Radius*Radius)/(s*s))*exp(Radius*Radius*guess*Ti/(Te*(s*s-Radius*Radius)))
				*sqrt(PI)*q_prime*(1.0+erf(q))+b*Radius*Radius*Ti/(Coeff*Te*(s*s-Radius*Radius)));
//		double b_prime = Coeff*(((s*s-Radius*Radius)/(s*s))*exp(Radius*Radius*guess*Ti/(Te*(s*s-Radius*Radius)))
//				*sqrt(PI)*q_prime*(1.0+erf(q))
//				+(Radius*Radius*Ti/(Te*s*s))*exp(Radius*Radius*guess*Ti/(Te*(s*s-Radius*Radius)))
//				*(exp(-q*q)+sqrt(PI)*q*(1.0+erf(q))));

		double fx = a-b-exp(guess);
		double fxprime = a_prime-b_prime-exp(guess);
		x1 = guess - fx/fxprime;

	}while(fabs(guess-x1)>1e-4);// >1e-2){
*/

	return guess;
}

