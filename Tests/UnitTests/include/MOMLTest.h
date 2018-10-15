#include "solveMOML.h"
#include <iostream>

// Test Modified orbital motion limited theory
void MOMLTest(){
	clock_t begin = clock();

	double Gamma = 1.0; // Ionization state of plasma

//	double Current0 = MOMLattractedCurrent(1.0,1.0,Gamma,Mp/Me,0.0);

//	std::cout << 0.0 << "\t" << Gamma << "\t" << Current0 << "\t" << Current0 << "\n";
	for( double TiTe(0.01); TiTe < 100; TiTe *= 1.5 ){
//		double Current = MOMLattractedCurrent(1.0,1.0,Gamma,Mp/Me,TiTe);
		for( double Gamma(0.01); Gamma < 100; Gamma *= 2 ){
			double Potential = solveMOML(TiTe,Gamma,Mp/Me);

			std::cout << "\n" << TiTe << "\t" << Gamma << "\t" << Mp/Me << "\t" << Potential;
		}
//		std::cout << TiTe << "\t" << Gamma << "\t" << Current << "\t" << Current << "\n";
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nMOML UnitTest completed in " << elapsd_secs << "s\n";
}
