#include "Constants.h"
#include "solveFlippedOML.h"
#include <iostream>

// Plot results using something
void FlippedOMLTest(){
	clock_t begin = clock();

	double converteVtoK(11600);

	double Potential = 0;

	double ne(1e18), ni(1e18);

	for(double Td(300); Td < 5e3; Td ++){ // Loop over temperatures
		double Sec = sec(Td/converteVtoK,'w'); 
		for( double Ti(10e4); Ti < 10e6; Ti *= 1.1 ){
			for( double Te(10e4); Te < 10e6; Te *= 1.1 ){
				double gammae = ne*exp(Potential)*sqrt(Kb*Te/(2*PI*Me));
				Potential = solveFlippedOML(Ti,Te,Td,ni,ne,Sec,Potential);
				std::cout << "\n" << Td << "\t" << Te/Ti << "\t" << Sec << "\t" << Potential;
			}
		}
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
//	std::cout << "\n\n*****\n\nUnitTest 4 completed in " << elapsd_secs << "s\n";
}
