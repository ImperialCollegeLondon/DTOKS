#include "solveMOMLWEM.h"
#include <iostream>

// Test orbital motion limited theory
void MOMLWEMTest(){
	clock_t begin = clock();

	double HeatCapacity = 5.0/3.0; // Ionization state of plasma
	double Ionization = 1.0; // Ionization state of plasma
//	for( double TiTe(0.01); TiTe < 100; TiTe *= 1.5 ){
	for( double TiTe(1.0); TiTe < 1.4; TiTe *= 1.5 ){
		for( double DeltaTot(0.01); DeltaTot < 20; DeltaTot *= 1.1 ){
			double Potential = solveMOMLWEM(TiTe,Ionization,Mp/Me,HeatCapacity,DeltaTot);
//			std::cout << "\n" << TiTe << "\t" << Ionization << "\t" << Mp/Me << "\t" << HeatCapacity 
//				<< "\t" << DeltaTot << "\t" << Potential;
		}
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nMOMLWEM UnitTest completed in " << elapsd_secs << "s\n";
}
