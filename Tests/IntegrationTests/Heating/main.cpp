// ## Heating:
// This contains tests which aim to establish the Temperature of the dust grain after a finite
// number of time steps.  Comparisons have been made to the expected analytic result.
// Those without analytical result have not been considered.
// The dust is started with the following parameters:
// Size = 1e-6;               // m
// Potential = 1;             // Normalised Potential
// Temp = 280;                // K
// In each of these tests, the 'Vapourise()' Heating model command is called. This means the
// heat equation is solved using a 4th order Runge-Kutta method until either:
// 1) The sample reaches thermal equilibrium
// 2) The sample vapourises (Mass < 10e-25)

#include <ctime>			// Time program
#include "ConstantHeatingTest.h"
#include "ConstantHeatEmissivTest.h"
#include "ConstantElectronPlasmaHeatingTest.h"
#include "ConstantPlasmaHeatingTest.h"
#include "ConstantPlasmaHeatingNeutralRecombTest.h"
#include "CompareConstEmissivTest.h"

int main(){

	char Element[4]={'W','B','F','G'}; // Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
	int out(0);
	for(int i(0); i < 4; i ++){
		std::cout << "\nElement["<<i<<"] is = ";
		std::cout << Element[i];

// 		Heating test 1, ConstantHeatingTest :
//		Test if constant heating power is comparable to analytic result
//		The equation being solved is:
//		dT/dt = Power / (Mass * HeatCapacity)
//		where Power is a constant, C
		std::cout << "\nRunning Test 1...\n";
		out = ConstantHeatingTest(Element[i]);	
		std::cout << "\nFinished Running Test 1!\n";

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";

// 		Heating test 2, ConstantHeatEmissivTest :
//		Test if constant heating with thermal radiation is comparable to analytic result
//		The equation being solved is:
//		dT/dt = Power / (Mass * HeatCapacity)
//		where Power is:
//		Power = C - SurfaceArea*sigma*epsilon*(Temperature^4-AmbientTemperature^4)
		std::cout << "\nRunning Test 2...\n";
		out = ConstantHeatEmissivTest(Element[i]);
		std::cout << "\nFinished Running Test 2!\n";

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";

// 		Heating test 3, ConstantElectronPlasmaHeatingTest : NOTE THIS DEVIATES AND I DON'T KNOW WHY EXACTLY! 
Test if constant heating with plasma heating and thermal radiation is comparable to analytic result
The equation being solved is:
dT/dt = Power / (Mass * HeatCapacity)
where Power is:
Power = C + SurfaceArea*(2*ElectronFlux*ElectronTemperature*Kb+2*NeutralFlux*NeutralTemperature*Kb)
		std::cout << "\nRunning Test 3...\n";
		out = ConstantElectronPlasmaHeatingTest(Element[i]);
		std::cout << "\nFinished Running Test 3!\n";

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";

// 		Heating test 4, ConstantPlasmaHeatingTest : NOTE THIS DEVIATES AND I DON'T KNOW WHY EXACTLY! 
//		Test if constant heating with plasma heating and thermal radiation is comparable to analytic result
//		The equation being solved is:
//		dT/dt = Power / (Mass * HeatCapacity)
//		where Power is:
//		Power = C + SurfaceArea*(2*ElectronFlux*ElectronTemperature*Kb+2*NeutralFlux*NeutralTemperature*Kb + IonFluxPower)
		std::cout << "\nRunning Test 4...\n";
		out = ConstantPlasmaHeatingTest(Element[i]);
		std::cout << "\nFinished Running Test 4!\n";

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";

// 		Heating test 5, ConstantPlasmaHeatingNeutralRecombTest : NOTE THIS DEVIATES AND I DON'T KNOW WHY EXACTLY! SAME PROBLEM
//		Test if constant heating with plasma heating, thermal radiation and neutral recomb is comparable to analytic result
//		The equation being solved is:
//		dT/dt = Power / (Mass * HeatCapacity)
//		where Power is:
//		Power = C + SurfaceArea*(2*ElectronFlux*ElectronTemperature*Kb+2*NeutralFlux*NeutralTemperature*Kb + IonFluxPower
//		                        + echarge*IonFlux*14.7 )
		std::cout << "\nRunning Test 5...\n";
		out = ConstantPlasmaHeatingNeutralRecombTest(Element[i]);
		std::cout << "\nFinished Running Test 5!\n";

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";
		std::cout << "\n\n*****\n\n\n";


// 		Heating test 6, CompareConstEmissivTest : NOTE THIS TEST IS SUPPOSED TO FAIL! 
//		Test if constant heating with thermal radiation with non-const emissivity is comparable to analytic result
//		The equation being solved is:
//		dT/dt = Power / (Mass * HeatCapacity)
//		where Power is:
//		Power = C - SurfaceArea*sigma*epsilon*(Temperature^4-AmbientTemperature^4)
		double Emissivs[4]={0.98397,9.73852e-06,0.240571,0.33886}; // Final temperature (TEST 3) Emissivities
		//		double Emissivs[4]={0.04,0.18,0.2,0.70}; // DTOKS emissivities
		std::cout << "\nRunning Test 6...\n";
//		out = CompareConstEmissivTest(Element[i],Emissivs[i]);
		std::cout << "\nFinished Running Test 6!\n";
		
		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";


		// Explanation for Other missing tests of terms:
		// TEE has no analytical solution as the solution involves exponential integrals
		// SEE has no result in standard mathematical functions
		// Qvap has extremely complicated analytical solutions involving the imaginary error function due to Antoinne eq.
		// Since ion flux depends on electron yield (i.e TEE and SEE), no general analytical solution for negative grains
		// For the case of constant potential, Q_(i) is constant

		std::cout << "\n\n*****\n\n\n";	
	}
	return 0;
}
