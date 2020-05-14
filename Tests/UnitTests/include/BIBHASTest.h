#include "solveBIBHAS.h"
#include <iostream>

// BIBHAS Charging Test:
// This test is designed to find the floating potential for arbitary sized dust grain
// The calculation follows the work by R. DE Bibhas,
// see R. DE Bibhas, Astrophys. Space Sci. 30, (1974).
void BIBHASTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\n#BIBHASTest\t\tVariables : " << Variables;

    double Beta = 3.0; // Ionization state of plasma
    double Ti = 1.0; // Ion Temperature in eV
    double Te = 1.0; // Electron Temperature in eV
    double Phi = 1.0; // Initial guess at potential
    double AtomicNumber = 1.0; // Charge state of element. 1.0 = First ionisation
    double AtomicMass = 1.0*Mp; // Atomic Mass Of Element
//    double Lambda = 1.1;
    std::cout << "\n\n#chi\ttau\tMi/Me\tAtomic Number\tLambda\n";
    for( double Lambda(0.0001); Lambda < 100.0; Lambda *= 1.01 ){    // Loop over magnetic field strengths
//        for( double Ti(0.01); Ti <= 1.0; Ti *= 10 ){
            double Potential = solveBIBHAS(Ti,Te,AtomicMass/Me,AtomicNumber,Lambda,-2.5);
            std::cout << "\n" << Potential << "\t" << Ti/Te << "\t" << AtomicMass/Me
                    << "\t" << AtomicNumber << "\t" << Lambda;
            std::cout << "\n\n";
//        }
    }
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nBIBHAS UnitTest completed in " << elapsd_secs << "s\n";
}
