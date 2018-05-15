#include "Constants.h"
#include "solveNegSchottkyOML.h"
#include "solvePosSchottkyOML.h"
#include <iostream>

// Plot results using something
void SchottkyOMLTest(){
	clock_t begin = clock();

	double converteVtoK(11600);
	double Potential = -0.1;

	double ne(1e18), ni(1e18);
	double Ti(10e4), Te(10e4);

	for(double Td(300); Td < 5500; Td ++){ // Loop over temperatures
		double Sec = sec(Td/converteVtoK,'W'); 
		for( double Ti(10e4); Ti < 10e6; Ti *= 2 ){
			for( double Te(10e4); Te < 10e6; Te *= 2 ){
				Potential = solveNegSchottkyOML(Ti,Te,Td,ni,ne,Sec,Potential);
				if ( Potential < 0 ){

//					std::cout << "\n" << Td << "\t" << Ti << "\t" << Te << "\t" << Potential;
					Potential = solvePosSchottkyOML(Td,Ti,Te,ne,Sec); // Quasi-neutrality assumed
					Potential = -Potential;
				}

				std::cout << "\n" << Td << "\t" << Ti << "\t" << Te << "\t" << Potential;
			}
		}
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nUnitTest 6 completed in " << elapsd_secs << "s\n";
}
