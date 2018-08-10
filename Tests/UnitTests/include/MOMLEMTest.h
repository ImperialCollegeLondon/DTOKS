#include "solveMOMLEM.h"
#include <iostream>

// Test orbital motion limited theory
void MOMLEMTest(){
	clock_t begin = clock();

	double guess = -1.0;	// Initial guess at potential
	double EmittedDensityRatio = 0.3;
	double EmittedTemperatureRatio = 0.5;
	double TiTe = 1.0;
//	for( double TiTe(0.01); TiTe < 100; TiTe *= 1.5 ){
	for( EmittedTemperatureRatio = 0.02; EmittedTemperatureRatio < 0.7; EmittedTemperatureRatio += 0.01 ){
		double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio,guess);
		std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
			<< "\t" << EmittedTemperatureRatio << "\t" << guess << "\t" << Potential 
			<< "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
	}
	EmittedDensityRatio = 0.3;
	EmittedTemperatureRatio = 0.5;
	TiTe = 1.0;

	std::cout << "\n";
	for( TiTe = 0.5; TiTe < 7.0; TiTe += 0.1 ){
		double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio,guess);
		std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
			<< "\t" << EmittedTemperatureRatio << "\t" << guess << "\t" << Potential
			<< "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
	}
	EmittedDensityRatio = 0.3;
	EmittedTemperatureRatio = 0.5;
	TiTe = 1.0;

	std::cout << "\n";
	for( EmittedDensityRatio = 0.01; EmittedDensityRatio < 0.5; EmittedDensityRatio += 0.01 ){
		double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio,guess);
		std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
			<< "\t" << EmittedTemperatureRatio << "\t" << guess << "\t" << Potential
			<< "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nMOMLEM UnitTest completed in " << elapsd_secs << "s\n";
}
