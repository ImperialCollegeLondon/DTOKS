#include "Functions.h"
#include "Constants.h"
#include <iostream>
#include <math.h>

// This test prints the value of the 'effective yield' from the Richardson-Dushmann formula (Without Schottky correction).
// This can also be readily compared to the publication
// Replicating work of "Dust in tokamaks: An overview of the physical model of the dust in tokamaks code"
// Bacharis, Minas Coppins, Michael Allen, John E.
// Page 3

// Plot results using matlab
void DeltaThermTest(){
	clock_t begin = clock();
	double electronDensity = 1e18;
	double NormalisedPotential = 2;
	for(double T(280); T < 5e3; T ++){ // Loop over temperatures
		double gammae = electronDensity*exp(-NormalisedPotential)*sqrt(echarge*10/(2*PI*Me));
		std::cout << "\n" << T << ", " 
			<< Richardson*pow(T,2)*exp(-(4.55*echarge)/(Kb*T))/(echarge*gammae);
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
//	std::cout << "\n\n*****\n\nUnitTest 3 completed in " << elapsd_secs << "s\n";
}
