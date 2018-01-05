// # UnitTests
// This directory contains small independant scripts which are designed to calculate one
// specific value or the value of an equation or function over a range of parameters. The
// tests are split into three main categories.
// Categories :            'Charging Tests', 'Force Tests', 'Heat Tests'
// 
// There is also one uncategorised test called 'BackscatterTest'.


#include "BackscatterTest.h"
#include "DeltaSecTest.h"
#include "DeltaThermTest.h"
#include "ChargingTimescales.h"
#include "DTOKSchargingTest.h"
#include "SchottkyOMLTest.h"
#include "OMLTest.h"
#include "MOMLTest.h"
#include "SOMLTest.h"
#include "SMOMLTest.h"
#include "IonNeutralDragTest.h"
#include "NeutralHeatingTest.h"
#include "EvaporativeCoolingTest.h"
#include "EvaporativeMassLossTest.h"

int main(){

	// Backscatter Unit Test:
	// This test prints the values of the fraction of backscattered energy and the fraction of back scattered particles
	// as calculated by the backscatter function from DTOKS. This can be readily compared to the results published in
	// "Dust in tokamaks: An overview of the physical model of the dust in tokamaks code." Physics of Plasmas, 17(4).
	// Bacharis, M., Coppins, M., & Allen, J. E. (2010). 
//	BackscatterTest();

	// ***** CHARGING TESTS ***** //

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

	// ***** ***** ***** ***** //

	// ***** FORCE TESTS ***** //

	// IonNeutralDrag Force Test: 
	// This test is designed to compare the relative magnitudes of the ion and neutral drag force.
	// The Ion Drag force is the one that was originally used in DTOKS while the neutral drag is from DUSTT,
	// See Pigarov, A. Y., Krasheninnikov, S. I., Soboleva, T. K. and Rognlien, T. D. (2005) 
	// ‘Dust-particle transport in tokamak edge plasmas’, Physics of Plasmas, 12(12), pp. 1–15. 
//	IonNeutralDragTest();

	// ***** ***** ***** ***** //

	// ***** HEAT TESTS ***** //

	// Neutral Heating Test: 
	// This test is designed to test the neutral heating and show the variation in magnitude for 
	// the neutral heating term
//	NeutralHeatingTest();

	// Evaporative Cooling Test: 
	// This test is designed to test the evaporative cooling as a function of temperature,
	// Here it is shown only for Tungsten 
//	EvaporativeCoolingTest();

	// Evaporative Mass Loss Test: 
	// This test is designed to test the evaporative mass loss as a function of temperature,
	// Here it is shown only for Tungsten 
	EvaporativeMassLossTest();

	// NeutralHeating Test: 
	// This test is designed to asses the heating power of neutral heating over a range of parameters
	// Heating Power is in W/m^2
 	NeutralHeatingTest();

	return 0;
}
