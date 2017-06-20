#include <ctime>			// Time program
#include "ConstantHeatingTest.h"
#include "ConstantHeatEmissivTest.h"
#include "ConstantElectronPlasmaHeatingTest.h"
#include "ConstantPlasmaHeatingTest.h"
#include "ConstantPlasmaHeatingNeutralRecombTest.h"

//#include "CompareConstEmissivTest.h"

int main(){

	char Element[4]={'W','B','F','G'}; // Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
	int out(0);
	for(int i(0); i < 4; i ++){
		std::cout << "\nElement["<<i<<"] is = ";
		std::cout << Element[i];

		// Test bout 1,
		// Integration test one for all materials
		// Test if constant heating is comparable to analytic result
		std::cout << "\nRunning Test 1...\n";
		out = ConstantHeatingTest(Element[i]);	
		std::cout << "\nFinished Running Test 1!\n";
		// Test bout 2,
		// Integration test two for all materials
		// Test if constant heating with thermal radiation is comparable to analytic result

//		std::cout << "\nRunning Test 2...\n";
//		out = ConstantHeatEmissivTest(Element[i]);
//		std::cout << "\nFinished Running Test 2!\n";

		// Test bout 3,
		// Integration test three for all materials
		// Test if constant heating with thermal radiation with non-const emissivity is comparable to analytic result
//		double Emissivs[4]={0.98397,9.73852e-06,0.240571,0.33886}; // Final temperature (TEST 3) Emissivities
//		double Emissivs[4]={0.04,0.18,0.2,0.70}; // DTOKS emissivities
//		out = CompareConstEmissivTest(Element[i],Emissivs[i]);

		// Test bout 4,
		// Integration test four for all materials
		// Test if constant heating with plasma heating and thermal radiation is comparable to analytic result
		// Note, emissivity is once again a const
//		out = ConstantElectronPlasmaHeatingTest(Element[i]);

		// Test bout 5,
		// Integration test four for all materials
		// Test if constant heating with plasma heating and thermal radiation
		// Note, emissivity is once again a const
//		out = ConstantPlasmaHeatingTest(Element[i]);

		// Test bout 6,
		// Integration test four for all materials
		// Test if constant heating with plasma heating, thermal radiation and neutral recombination
		// Note, emissivity is once again a const
//		out = ConstantPlasmaHeatingNeutralRecombTest(Element[i]);

		// Other terms:
		// TEE has no analytical solution as the solution involves exponential integrals
		// SEE has no result in standard mathematical functions
		// Qvap has extremely complicated analytical solutions involving the imaginary error function due to Antoinne eq.
		// Since ion flux depends on electron yield (i.e TEE and SEE), no general analytical solution for negative grains
		// For the case of constant potential, Q_(i) is constant

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";
		std::cout << "\n\n*****\n\n\n";

		std::cout << "\n\n*****\n\n\n";	
	}
	return 0;
}
