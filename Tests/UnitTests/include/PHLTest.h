#include "solvePHL.h"
#include <iostream>

// Test Patacchini, Hutchinson and Lapenta charging of a sphere in a collisionless magnetoplasma
// L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).
// Assuming that Lambda_s >> r_p, i.e dust grain is small.
void PHLTest(){
	clock_t begin = clock();

	double Beta = 3.0; // Ionization state of plasma
	double Ti = 1.0; // Ion Temperature in eV
	double Te = 1.0; // Electron Temperature in eV
	double Phi = 1.0; // Initial guess at potential
	double Lambda = 0.0; // Debye length relative to probe radius
	double AtomicNumber = 1.0; // Charge state of element. 1.0 = First ionisation
	double AtomicMass = 1.0*Mp; // Atomic Mass Of Element
	for( double Beta(0); Beta < 1000; Beta += 0.1 ){	// Loop over magnetic field strengths
		double Potential = solvePHL(Phi,Beta,Ti/Te,AtomicMass/Me,AtomicNumber,Lambda);
		std::cout << "\n" << Phi << "\t" << Beta*sqrt(Me*AtomicNumber*AtomicNumber*Te/(Ti*AtomicMass)) << "\t" << Ti/Te << "\t" << AtomicMass/Me << "\t" << AtomicNumber << "\t" << Lambda << "\t" << Potential;
		std::cout << "\t" << (1.0-(AtomicNumber/(Ti/Te))*Potential) 
			<< "\t" << (1.0/sqrt(Ti/Te))*sqrt(AtomicMass/Me)*solveI_curr(Potential,Beta,Ti/Te,AtomicMass/Me,1.0,Lambda);
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nPHL UnitTest completed in " << elapsd_secs << "s\n";
}
