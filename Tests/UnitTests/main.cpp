#include "BackscatterTest.h"

// CHARGING TESTS
#include "DeltaSecTest.h"
#include "DeltaThermTest.h"
#include "ChargingTimescales.h"
#include "DTOKSchargingTest.h"
#include "DTOKSwellchargingTest.h"
#include "SchottkyOMLTest.h"
#include "OMLTest.h"
#include "MOMLTest.h"
#include "SOMLTest.h"
#include "SMOMLTest.h"

// FORCE TESTS
#include "HybridIonDrag.h"
#include "FortovIonDrag.h"

int main(){

	// Backscatter Unit Test:
	// This test prints the values of the fraction of backscattered energy and the fraction of back scattered particles
	// as calculated by the backscatter function from DTOKS. This can be readily compared to the results published in
	// the DTOKS papers
//	BackscatterTest();

	// ***** 	CHARGING TESTS		***** //

	// Delta Sec Unit Test:
	// This test prints the value of the empirical function calculating the yield due to secondary electron emission.
	// The result can be readily compared to the publication
	// Replicating work of "Dust in tokamaks: An overview of the physical model of the dust in tokamaks code"
	// Bacharis, Minas Coppins, Michael Allen, John E.
	// Page 2 & 3
//	DeltaSecTest();

	// Delta Therm Unit Test:
	// This test prints the value of the 'effective yield' from the Richardson-Dushmann formula (Without Schottky correction).
	// This can also be readily compared to the publication
        // Replicating work of "Dust in tokamaks: An overview of the physical model of the dust in tokamaks code"
        // Bacharis, Minas Coppins, Michael Allen, John E.
        // Page 2 & 3
//	DeltaThermTest();

	// Charging Timescale Test:
	// This test is used to verify that the timestep as calculated by Krasheninnikovs is always smaller
	// than the electron plasma frequency. In practice, this is found to not be perfectly true, there exist
	// extreme conditions where this is not the case.
	// Time step based on the formulation by Krasheninnikov
	// Smirnov, R. D., Pigarov, A. Y., Rosenberg, M., Krasheninnikov, S. I., & Mendis, D. a. (2007). 
	// Modelling of dynamics and transport of carbon dust particles in tokamaks. 
	// Plasma Physics and Controlled Fusion, 49(4), 347–371.	
//	ChargingTimescales();
	
	// DTOKS Charging Test:
	// This test output the potential as calculated by the DTOKS solution to the OML equation.
	// The form of this is such that it depends on the sign of the potential and the magnitude of the total electron emission.
	// This test was used to show the discontinuity in the potential when the total electron emission yield approaches 1
	// Two different conditional formulations of the problem are made and their differences highlighted by this test
//	DTOKSchargingTest();

	// This test outputs the potential as calculated by part of the DTOKS solution to the OML equation, considering only a 
	// scenario with a potential well. This removes the discontinuities observed in the previous model.
//	DTOKSwellchargingTest();
	
	// OML Charging Test:
	// This test is designed to find the floating potential for small dust grains in a stationary plasma following OML theory.
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
	OMLTest();

	// MOML Charging Test:
	// This test is designed to find the floating potential for large dust grains in a stationary plasma following MOML theory.
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
//	MOMLTest();

	// SOML Charging Test:
	// This test is designed to find the floating potential for small dust grains in a flowing plasma following SOML theory.
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
//	SOMLTest();

	// SMOML Charging Test:
	// This test is designed to find the floating potential for large dust grains in a flowing plasma following SMOML theory.
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
//	SMOMLTest();

	// Schottky OML Charging Test:
	// This test was made to see what the solution is for electron emission with Schottky correction where the potential of the
	// dust grain is accounted for. The minimisation of the positive solution in C++ is not stable and gives an incorrect 
	// answer. Switching to matlab minimisation function, some weird things happen but, in principle, I showed that the 
	// function could be minimised.
//	SchottkyOMLTest();

	// SchottkyMOML Charging Test: DOESN'T WORK!
	// This test is designed to find the floating potential for large negative dust grains with electron emission 
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
//	SchottkyMOMLTest();


	// ***** 	FORCE TESTS		***** //
	// Hybrid Ion Drag Test
	// This test is designed to test the magnitude of the HybridIonDrag force as formulated in the paper given below.
	// The Hybrid Ion Drag model is a function of the plasma species temperature ratio, the ion mach number, ion density,
	// the dust grain potential and the dust grain radius
	// Khrapak, S. A., Ivlev, A. V., Zhdanov, S. K., & Morfill, G. E. (2005). Hybrid approach to the ion drag force. 
	// Physics of Plasmas, 12(4), 1–8. https://doi.org/10.1063/1.1867995 
//	HybridIonDragTest();

	// Fortov et al./DTOKS Ion Drag Test
	// This test is designed to test the magnitude of the drag force as formulated in the paper given below by fortov.
	// The Hybrid Ion Drag model is a function of the temperature of electrons and ions, the ion mach number, ion density,
	// the dust grain potential and the dust grain radius
	// Fortov, V. E., Ivlev, A. V., Khrapak, S. A., Khrapak, A. G., & Morfill, G. E. (2005).
	//  Complex (dusty) plasmas: Current status, open issues, perspectives. Physics Reports, 421(1–2), 1–103. 
	// https://doi.org/10.1016/j.physrep.2005.08.007
//	FortovIonDragTest();

	return 0;
}
