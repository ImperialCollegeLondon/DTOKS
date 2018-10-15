#include "solveTH.h"
#include <iostream>

// TH Charging Test:
// This test is designed to find the floating potential for small dust grains in magnetised plasmas 
// The calculation follows the work by Drew Thomas and Josh Holgate and implements a semi-empirical model
// see D. M. Thomas and J. T. Holgate, ArXiv Prepr. (2016). for details
void THTest(){
	clock_t begin = clock();

	double Beta = 3.0; // Ionization state of plasma
	double Ti = 1.0; // Ion Temperature in eV
	double Te = 1.0; // Electron Temperature in eV
	double Phi = 1.0; // Initial guess at potential
	double AtomicNumber = 1.0; // Charge state of element. 1.0 = First ionisation
	double AtomicMass = 1.0*Mp; // Atomic Mass Of Element
//	double Lambda = 1.1;

	for( double Lambda(15.0); Lambda > 1.5; Lambda -= 12.0 ){	// Loop over magnetic field strengths
		for( double Beta(0); Beta < 20; Beta += 0.0005 ){	// Loop over magnetic field strengths

//		for( double Lambda(1.0/15.0); Lambda < 6.0/15.0; Lambda += 4.0/15.0 ){	// Loop over magnetic field strengths
			double Potential = solveTH(Beta,Ti/Te,AtomicMass/Me,AtomicNumber,Lambda);
			std::cout << "\n" << Beta << "\t" << Potential << "\t" << Ti/Te << "\t" << AtomicMass/Me
					<< "\t" << AtomicNumber << "\t" << Lambda;
		}
		std::cout << "\n\n";
	}
	double Lambda = 1.0;
	for( double Beta(0); Beta < 20; Beta += 0.001 ){	// Loop over magnetic field strengths
		double Potential = solveTH(Beta,Ti/Te,AtomicMass/Me,AtomicNumber,Lambda);
		std::cout << "\n" << Beta << "\t" << Potential << "\t" << Ti/Te << "\t" << AtomicMass/Me
				<< "\t" << AtomicNumber << "\t" << Lambda;
	}
	std::cout << "\n\n";

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nTH UnitTest completed in " << elapsd_secs << "s\n";
}
