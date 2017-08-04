#include "solveSMOML.h"
#include <iostream>

// Test Shifted modified orbital motion limited theory
void SMOMLTest(){
	clock_t begin = clock();

//	double U = 1.0; // Plasma Flow speed
	double Gamma = 3.0;	// Heat Capacity Ratio
	for( double TiTe(0.01); TiTe < 50; TiTe *= 1.5 ){ // Loop over temperature ratios
		for( double U(0.1); U < 3; U *= 1.5 ){ // Loop over plasma flow speed
			double Potential = solveSMOML(TiTe,Gamma,U,Mp/Me);
			std::cout << "\n" << TiTe << "\t" << U << "\t" << Mp/Me << "\t" << Potential;
		}
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nSMOML UnitTest completed in " << elapsd_secs << "s\n";
}
