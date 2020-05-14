#include "Functions.h"
#include <iostream>

// This test prints the value of the Maxwellian funciton for different energies and temperatures

// Plot results using GNUplot
void MaxwellianTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\n#MaxwellianTest\t\tVariables : " << Variables;
    // Replicate results of Minas' paper:
    for( double Energy(0.01); Energy < 100; Energy *= 1.1 ){
        for( double Temperature(0.01); Temperature < 100; Temperature *= 1.1 ){
            std::cout << "\n" << Energy << " " << Temperature << " " <<  maxwellian(Energy,Temperature);
        }
    }
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
//    std::cout << "\n\n*****\n\nUnitTest 2 completed in " << elapsd_secs << "s\n";
}
