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
	double electronDensity = 1e15;	// Density in m^(-3)
	double Pot = 1.0;

	// WorkFunctions, from https://en.wikipedia.org/wiki/Work_function#cite_note-12, 
        // CRC Handbook of Chemistry and Physics version 2008, p. 12–114:
        // W = 4.32 - Be = 4.98 - C = 5.0 - Fe = 4.67

	double Wf = 1.0;	// Work function in eV
	double Te = 1;			// Electron Temperature in eV
//      for(double Wf(1.0); Wf < WorkFunction; Wf ++ ){
//      for(double Pot(0.0); Pot < NormalisedPotential; Pot ++ ){
	for(double T(280); T < 5.5e3; T ++){ // Loop over temperatures
		std::cout << "\n" << T << ", "
			<< Richardson*pow(T,2)*exp(-(Wf*echarge)/(Kb*T))/(echarge* electronDensity*exp(-Pot)*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-(Wf*echarge)/(Kb*T))/(echarge* electronDensity*exp(-(Pot+1.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-(Wf*echarge)/(Kb*T))/(echarge* electronDensity*exp(-(Pot+4.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-((Wf+1.0)*echarge)/(Kb*T))/(echarge* electronDensity*exp(-Pot)*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-((Wf+1.0)*echarge)/(Kb*T))/(echarge* electronDensity*exp(-(Pot+1.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-((Wf+1.0)*echarge)/(Kb*T))/(echarge* electronDensity*exp(-(Pot+4.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-((Wf+4.0)*echarge)/(Kb*T))/(echarge* electronDensity*exp(-Pot)*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-((Wf+4.0)*echarge)/(Kb*T))/(echarge* electronDensity*exp(-(Pot+1.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-((Wf+4.0)*echarge)/(Kb*T))/(echarge* electronDensity*exp(-(Pot+4.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-((Wf*echarge)/(Kb*T)-Pot))/(echarge*electronDensity*exp(-Pot)*sqrt(echarge*Te/(2*PI*Me))) << ", " 
			<< Richardson*pow(T,2)*exp(-((Wf*echarge)/(Kb*T)-(Pot+1.0)))/(echarge*electronDensity*exp(-(Pot+1.0))*sqrt(echarge*Te/(2*PI*Me))) << ", " 
			<< Richardson*pow(T,2)*exp(-((Wf*echarge)/(Kb*T)-(Pot+4.0)))/(echarge*electronDensity*exp(-(Pot+4.0))*sqrt(echarge*Te/(2*PI*Me))) << ", " 
			<< Richardson*pow(T,2)*exp(-(((Wf+1.0)*echarge)/(Kb*T)-Pot))/(echarge*electronDensity*exp(-Pot)*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-(((Wf+1.0)*echarge)/(Kb*T)-(Pot+1.0)))/(echarge*electronDensity*exp(-(Pot+1.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-(((Wf+1.0)*echarge)/(Kb*T)-(Pot+4.0)))/(echarge*electronDensity*exp(-(Pot+4.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-(((Wf+4.0)*echarge)/(Kb*T)-Pot))/(echarge*electronDensity*exp(-Pot)*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-(((Wf+4.0)*echarge)/(Kb*T)-(Pot+1.0)))/(echarge*electronDensity*exp(-(Pot+1.0))*sqrt(echarge*Te/(2*PI*Me))) << ", "
			<< Richardson*pow(T,2)*exp(-(((Wf+4.0)*echarge)/(Kb*T)-(Pot+4.0)))/(echarge*electronDensity*exp(-(Pot+4.0))*sqrt(echarge*Te/(2*PI*Me)));

	}
//      }
//      }

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
//	std::cout << "\n\n*****\n\nUnitTest 3 completed in " << elapsd_secs << "s\n";
}
