#include "solveSchottkyMOML.h"
#include <iostream>

// Test orbital motion limited theory
void SchottkyMOMLTest(){
	clock_t begin = clock();

	double Gamma = 3.0;
	for( double Td(300); Td < 5555; Td += 1 ){			// Loop over Temperature in K
		for( double Te(1); Te < 2; Te *= 2 ){			// Loop over Temperature in ev
			for( double Ti(1); Ti < 2; Ti *= 2 ){		// Loop over Temperature in ev
				double Potential = solveSchottkyMOML(Td,Te,Ti,Gamma,Mp/Me);
				std::cout << "\n" << Td << "\t" << Te << "\t" << Ti << "\t" << Mp/Me << "\t" << Potential;
			}
		}
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nMOML UnitTest completed in " << elapsd_secs << "s\n";
}
