#include <iostream>
#include "Constants.h"
#include "MathHeader.h" 

// Test orbital motion limited theory
void NeutralHeatingTest(){
	clock_t begin = clock();

	double ConvertKelvsToeV(8.621738e-5);

	for( double Density(0.01e18); Density < 100e18; Density *= 2 ){
		for( double Temp(0.01/ConvertKelvsToeV); Temp < 100/ConvertKelvsToeV; Temp *= 2 ){
			double NeutralPower = 2*Density*sqrt(Kb*Temp/Mp)*Temp*Kb;	// J / m^2
			std::cout << Density << "\t" << Temp << "\t" << NeutralPower << "\n";
		}
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
	
	std::cout << "\n\n*****\n\nNeutralHeating UnitTest completed in " << elapsd_secs << "s\n";
}
