#include "solveSOML.h"
#include <iostream>

// Test Shifted orbital motion limited theory
void SOMLTest(){
	clock_t begin = clock();

//	double U = 1.0; // Plasma Flow speed
//	for( double TiTe(0.01); TiTe < 10; TiTe *= 1.5 ){ // Loop over temperature ratios
	double TiTe=1.0;
		for( double U(0.1); U < 5; U *= 1.2 ){ // Loop over plasma flow speed
			double Potential = solveSOML(TiTe,U*sqrt(2.0),Mp/Me);
			std::cout << "\n" << TiTe << "\t" << U << "\t" << Mp/Me << "\t" << Potential;
		}
//	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nSOML UnitTest completed in " << elapsd_secs << "s\n";
}
