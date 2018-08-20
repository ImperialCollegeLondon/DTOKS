#include "Constants.h"
#include "Functions.h"

// Solve the Potential for a sphere in a collisionless magnetoplasma following
// L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).
// Assuming that Lambda_s >> r_p, i.e dust grain is small.
double solvePHL(double Phi, double Beta, double TemperatureRatio, double MassRatio, double AtomicNumber, double DebyeLength){ 
	double z = Beta/(1+Beta);
	double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+1.150*z*z*z*z*z;

	double eta = -(Phi/Beta)*(1+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
	double w(1.0);
	// Ensure that w is well defined
	if( Beta == 0.0 || eta == -1.0 ){
		w = 1.0;
	}else if( std::isnan(eta) ){
		std::cout << "\nWarning! w being set to 1.0 (Assuming high B field limit) but Phi/Beta is nan while Beta != 0.";
		w = 1.0;
	}else{
		w = eta/(1+eta);
	} 
	
	double A = 0.678*w+1.543*w*w-1.212*w*w*w;

	double Potential = TemperatureRatio/AtomicNumber
		-LambertW((A+(1.0-A)*i_star)*sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/AtomicNumber)/AtomicNumber);

	// Iterate until we converge on the correct current
	while( fabs(Phi - Potential) >= 0.005 ){
		Phi = solvePHL(Potential, Beta, TemperatureRatio, MassRatio, AtomicNumber, DebyeLength);
	
		eta = -(Phi/Beta)*(1+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
		double w(1.0);
		if( Beta == 0.0 || eta == -1.0 ){
			w = 1.0;
		}else if( std::isnan(eta) ){
			std::cout << "\nWarning! w being set to 1.0 (Assuming high B field limit) but Phi/Beta is nan while Beta != 0.";
			w = 1.0;
		}else{
			w = eta/(1+eta);
		} 

	        A = 0.678*w+1.543*w*w-1.212*w*w*w;
	
	        Potential = TemperatureRatio/AtomicNumber
                -LambertW((A+(1.0-A)*i_star)*sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/AtomicNumber)/AtomicNumber);
	}
	return Potential;
}

// Calculate the electron current to a sphere
// L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).
// Assuming that Lambda_s >> r_p, i.e dust grain is small.
double solveI_curr(double Phi, double Beta, double TemperatureRatio, double MassRatio, double AtomicNumber, double DebyeLength){ 
	double z = Beta/(1+Beta);
	double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+1.150*z*z*z*z*z;

	double eta = -(Phi/Beta)*(1+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
	double w(1.0);
	if( Beta == 0.0 || eta == -1.0 ){
		w = 1.0;
	}else if( std::isnan(eta) ){
		std::cout << "\nWarning! w being set to 1.0 (Assuming high B field limit) but Phi/Beta is nan while Beta != 0.";
		w = 1.0;
	}else{
		w = eta/(1+eta);
	} 

	double A = 0.678*w+1.543*w*w-1.212*w*w*w;

	return (A+(1.0-A)*i_star)*exp(Phi);
}

