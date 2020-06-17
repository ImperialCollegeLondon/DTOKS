#include <iostream>
#include "Constants.h"
#include "MathHeader.h" 

// Test orbital motion limited theory
void NeutralHeatingTest(unsigned int Variables){
    clock_t begin = clock();

    std::cout << "\n#NeutralHeatingTest\t\tVariables : " << Variables;
    double ConvertKelvsToeV(8.621738e-5);
    
    int i(0), j(0);
    for( double Density(0.01e18); Density < 100e18; Density *= 1.5 ){
        i = i + 1;
        for( double Temp(0.01/ConvertKelvsToeV); Temp < 100/ConvertKelvsToeV; Temp *= 1.5 ){
            if( i== 1 ) j = j + 1;
            double NeutralPower = 2*Density*sqrt(Kb*Temp/Mp)*Temp*Kb;    // J / m^2
            std::cout << Density << "\t" << Temp*ConvertKelvsToeV << "\t" << NeutralPower << "\n";
        }
    }

    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
    
    std::cout << "\n\n*****\n\nNeutralHeating UnitTest completed in " << elapsd_secs << "s\n" << i << " : " << j << "\n";
}
