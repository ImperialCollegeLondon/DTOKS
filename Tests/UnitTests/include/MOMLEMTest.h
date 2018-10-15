#include "solveMOMLEM.h"
#include <iostream>
#include <limits>

// Test orbital motion limited theory
void MOMLEMTest(){
	clock_t begin = clock();

	double guess = -1.0;	// Initial guess at potential
	double EmittedDensityRatio = 0.3;
	double EmittedTemperatureRatio = 0.5;
	double TiTe = 1.0;
//	for( double TiTe(0.01); TiTe < 100; TiTe *= 1.5 ){

  	std::cout.precision(std::numeric_limits<double>::max_digits10);
	double CritVal = FindCriticalVal(TiTe,EmittedDensityRatio,Mp/Me,'d');
//	std::cout << "\n" << CritVal;

	exploreWellCase(TiTe,0.2,EmittedTemperatureRatio,1.0/(Mp/Me));
	for( EmittedTemperatureRatio = 0.2; EmittedTemperatureRatio < 0.201; EmittedTemperatureRatio += 0.01 ){
//		double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio);

//		std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
//			<< "\t" << EmittedTemperatureRatio << "\t" << Potential 
//			<< "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
//		std::cout << "\t" << Potential;
	}
	EmittedTemperatureRatio = 0.5;
/*
	std::cout << "\n";
	CritVal = FindCriticalVal(EmittedDensityRatio,EmittedTemperatureRatio,Mp/Me,'t');
	std::cout << "\n" << CritVal;
	for( TiTe = 0.5; TiTe < 7.0; TiTe += 0.1 ){
		double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio);
//		std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
//			<< "\t" << EmittedTemperatureRatio << "\t" << guess << "\t" << Potential
//			<< "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
//		std::cout << "\t" << Potential;
	}
	TiTe = 1.0;

	std::cout << "\n";
	CritVal = FindCriticalVal(EmittedTemperatureRatio,TiTe,Mp/Me,'c');
	std::cout << "\n" << CritVal;
	for( EmittedDensityRatio = 0.01; EmittedDensityRatio < 0.5; EmittedDensityRatio += 0.01 ){
		double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio);
//		std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
//			<< "\t" << EmittedTemperatureRatio << "\t" << guess << "\t" << Potential
//			<< "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
//		std::cout << "\t" << Potential;
	}
*/
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nMOMLEM UnitTest completed in " << elapsd_secs << "s\n";
}
