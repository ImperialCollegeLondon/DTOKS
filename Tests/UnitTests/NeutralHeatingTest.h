#include <iostream>
#include "Constants.h"
#include "MathHeader.h" 

double NeutralHeating(double NeutralDensity, double NeutralTemp, double radius){
	return 4*PI*pow(radius,2)*2*NeutralDensity*sqrt(Kb*NeutralTemp/(2*PI*Mp))*NeutralTemp*Kb;
}

// Test orbital motion limited theory
void NeutralHeatingTest(){
	clock_t begin = clock();

	double ConvertKelvsToeV(8.621738e-5);

	threevector PlasmaVel(0.0,0.0,0.0);
	threevector DustVel(0.0,0.0,10.0);
	threevector RelativeVelocity = PlasmaVel-DustVel;
	double ElectronTemp = 10/ConvertKelvsToeV;	// 10ev converted to Kelvin
	double Radius = 1e-6;				// metres, radius of dust
	double Potential = 1;

	for( double Density(0.01e18); Density < 100e18; Density *= 2 ){
		for( double Temp(0.01/ConvertKelvsToeV); Temp < 100/ConvertKelvsToeV; Temp *= 2 ){
			double NeutralHeat = NeutralHeating(Density,Temp,Radius);
			std::cout << Density << "\t" << Temp*ConvertKelvsToeV << "\t"  << Radius << "\t" << NeutralHeat << "\n";
		}
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIonNeutralDrag UnitTest completed in " << elapsd_secs << "s\n";
}
