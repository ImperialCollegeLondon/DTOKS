#include "solveOML.h"
#include <iostream>

// Test orbital motion limited theory
void OMLTest(){
	clock_t begin = clock();

	double Z = 1.0; // Ionization state of plasma
	for( double TiTe(0.01); TiTe < 100; TiTe *= 2 ){
			double Potential = solveOML(TiTe,Z,Mp/Me);
			std::cout << "\n" << TiTe << "\t" << Z << "\t" << Mp/Me << "\t" << Potential;
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nOML UnitTest completed in " << elapsd_secs << "s\n";
}
