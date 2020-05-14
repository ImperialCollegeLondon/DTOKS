#include "Functions.h"
#include "Constants.h"
#include <iostream>

// This test prints the value of the empirical function calculating the yield due to secondary electron emission.
// The result can be readily compared to the publication

// Replicating work of "Dust in tokamaks: An overview of the physical model of the dust in tokamaks code"
// Bacharis, Minas Coppins, Michael Allen, John E.
// Page 2 & 3

// Plot results using matlab
void DeltaSecTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\n#DeltaSecTest\t\tVariables : " << Variables;
    // Replicate results of Minas' paper:
    std::cout << "\n\n#Te (eV)\tdelta_sec\n";
    for(double i(1); i < 10e3; i ++) 
        std::cout << "\n" << i << ", " <<  sec(i,'w'); // Convert from K to ev
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
//    std::cout << "\n\n*****\n\nUnitTest 2 completed in " << elapsd_secs << "s\n";
}
