#include "solveMOML.h"
#include <iostream>

// Test orbital motion limited theory
void MOMLTest(){
	clock_t begin = clock();

	double Gamma = 3.0; // Ionization state of plasma
	for( double TiTe(0.01); TiTe < 100; TiTe *= 2 ){
			double Potential = solveMOML(TiTe,Gamma,Mp/Me);
			std::cout << "\n" << TiTe << "\t" << Gamma << "\t" << Mp/Me << "\t" << Potential;
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nMOML UnitTest completed in " << elapsd_secs << "s\n";
}
