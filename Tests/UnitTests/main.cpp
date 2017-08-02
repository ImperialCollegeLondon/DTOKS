#include "BackscatterTest.h"
#include "DeltaSecTest.h"
#include "DeltaThermTest.h"
#include "ChargingTimescales.h"
#include "DTOKSchargingTest.h"
#include "SchottkyOMLTest.h"
#include "OMLTest.h"
#include "MOMLTest.h"
#include "SOMLTest.h"

int main(){

	// Backscatter Unit Test:
	// This test prints the values of the fraction of backscattered energy and the fraction of back scattered particles
	// as calculated by the backscatter function from DTOKS. This can be readily compared to the results published in
	// the DTOKS papers
//	BackscatterTest();

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
	
	// Schottky OML Charging Test:
	// This test was made to see what the solution is for electron emission with Schottky correction where the potential of the
	// dust grain is accounted for. The minimisation of the positive solution in C++ is not stable and gives an incorrect 
	// answer. Switching to matlab minimisation function, some weird things happen but, in principle, I showed that the 
	// function could be minimised.
//	SchottkyOMLTest();

	// OML Charging Test:
	// This test is designed to find the floating potential for small dust grains in a stationary plasma following OML theory.
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
//	OMLTest();

	// MOML Charging Test:
	// This test is designed to find the floating potential for large dust grains in a stationary plasma following MOML theory.
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
//	MOMLTest();

	// SOML Charging Test:
	// This test is designed to find the floating potential for small dust grains in a flowing plasma following SOML theory.
	// This employs an approximate series expansion to the Lambert W function to find the floating potential
	SOMLTest();
	return 0;
}
