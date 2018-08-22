#include "Constants.h"
#include <iostream>

// This test outputs the potential as calculated by part of the DTOKS solution to the OML equation, considering a potential well.
void HybridIonDragTest(){
	clock_t begin = clock();

	// TEST TO CALCULATE THE DTOKS FLOATING POTENTIAL FOR CONSTANT ELECTRON DENSITY AND ELECTRON TEMPERATURE

	double Radius = 1e-6; 		// Dust radius in m
	double IonDensity = 1e18; 	// Ion Density in m^-3
	double PlasmaVel = 1e4;		// PlasmaVelocity in m s^-^1
	for( double u(0.01); u < 100; u *= 1.1){	// Normalised ion flow velocity
		for( double TeTi(0.01); TeTi < 100; TeTi *= 1.1 ){	// temperature ratio Te/Ti
			for( double potential(0.0); potential < 5; potential += 0.1 ){	// Normalised dust potential

				double z = potential*4.0*PI*epsilon0;
				double CoulombLogarithm = 17.0;		// Approximation of coulomb logarithm	
			
				double Coefficient = sqrt(2*PI)*pow(Radius,2.0)*IonDensity*Mp*pow(PlasmaVel,2.0);
				double term1 = sqrt(PI/2.0)*erf(u/sqrt(2))*(1.0+u*u+(1.0-(1.0/(u*u)))
						*(1.0+2*TeTi*z)+4*z*z*TeTi*TeTi*CoulombLogarithm/(u*u));
				double term2 = (1.0/u)*exp(-u*u/2.0)*(1.0+2.0*TeTi*z+u*u-4*z*z*TeTi*TeTi*CoulombLogarithm);	
				
				double HybridDrag = Coefficient*(term1+term2);

				std::cout << "\n" << u << "\t" << TeTi << "\t" << potential << "\t" << HybridDrag;
			}
		}
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nHybridIonDrag UnitTest completed in " << elapsd_secs << "s\n";
}
