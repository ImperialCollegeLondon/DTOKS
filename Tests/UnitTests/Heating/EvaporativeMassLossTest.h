#include <iostream>
#include "Constants.h"
#include "MathHeader.h" 

// Test orbital motion limited theory
void EvaporativeMassLossTest(){
	clock_t begin = clock();

	double ConvertKelvsToeV(8.621738e-5);

	double AmbientPressure = 0.0;
	double AtomicMass = 0.18384; 	// kg/mol

	for( double Temp(0.1/ConvertKelvsToeV); Temp < 100/ConvertKelvsToeV; Temp *= 1.1 ){
		double VapourPressure = 101325*pow(10,2.945 - 44094/Temp + 1.3677*log10(Temp)); // http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf

		double EvapMassLoss = AvNo*(VapourPressure-AmbientPressure)/
			sqrt(2*PI*AtomicMass*R*Temp)*(AtomicMass/AvNo); // Kg/m^2

		std::cout << Temp << "\t" << VapourPressure << "\t" <<  EvapMassLoss << "\n";
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
	
	std::cout << "\n\n*****\n\nEvaporativeCooling UnitTest completed in " << elapsd_secs << "s\n";
}
