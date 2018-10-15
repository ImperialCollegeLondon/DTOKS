#include <iostream>
#include "Constants.h"
#include "MathHeader.h" 

// Test orbital motion limited theory
void EvaporativeCoolingTest(){
	clock_t begin = clock();

	double ConvertKelvsToeV(8.621738e-5);

//	double BondEnergy = 774;	// kJ/mol, Tungsten
	double BondEnergy = 320.3;	// kJ/mol, Beryllium
//	double BondEnergy = 13.810;	// kJ/mol, Iron
//	double BondEnergy = 345;	// kJ/mol, Graphite
	double AmbientPressure = 0.0;
//	double AtomicMass = 0.18384; 	// kg/mol, Tungsten
	double AtomicMass = 0.009012182;	// kg/mol, Beryllium
//	double AtomicMass = 0.055845;		// kg/mol, Iron
//	double AtomicMass = 0.0120107;		// kg/mol, Graphite//	double MeltingTemp = 3422;	// K, Tungsten 
//	double BoilingTemp = 5555;	// K, Tungsten 
	double MeltingTemp = 1560;	// K, Beryllium 
	double BoilingTemp = 2742;	// K, Beryllium 
//	double MeltingTemp = 1811;	// K, Iron
//	double BoilingTemp = 3134;	// K, Iron
//	double MeltingTemp = 4500;	// K, Graphite
//	double BoilingTemp = 400;	// K, Graphite

	for( double Temp(MeltingTemp); Temp < BoilingTemp; Temp += 1 ){
		// Tungsten, http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
//		double VapourPressure = 101325*pow(10,2.945 - 44094/Temp + 1.3677*log10(Temp));
		// Beryllium, http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
		double VapourPressure = pow(10,5.786 -15731/Temp);
		// Iron 
//		double VapourPressure = pow(10,11.353 - 19574/Temp);
		// Graphite
//		double VapourPressure = 0.0;

		double MaxwellEnergy = (3*Kb*Temp/(2*1000)); // Converted to kJ.

		double EvapPower = 1000*AvNo*(VapourPressure-AmbientPressure)/
                        sqrt(2*PI*AtomicMass*R*Temp)*(MaxwellEnergy+BondEnergy/AvNo); // J/m^2

		std::cout << Temp << "\t" << VapourPressure << "\t" <<  EvapPower << "\n";
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
	
	std::cout << "\n\n*****\n\nEvaporativeCooling UnitTest completed in " << elapsd_secs << "s\n";
}
