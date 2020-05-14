#include <iostream>
#include "Constants.h"
#include "MathHeader.h" 

// Test orbital motion limited theory
void NeutralRecombTest(unsigned int Variables){
    clock_t begin = clock();

    std::cout << "\n#NeutralRecombinationTest\t\tVariables : " << Variables;
    double ConvertKelvsToeV(8.621738e-5);
    
    int i(0), j(0);
    for( double Density(0.01e18); Density < 100e18; Density *= 1.5 ){
        i = i + 1;
        for( double DustTemp(0.01/ConvertKelvsToeV); DustTemp < 100/ConvertKelvsToeV; DustTemp *= 1.5 ){
            for( double Ti(1/ConvertKelvsToeV); Ti < 10000/ConvertKelvsToeV; Ti *= 1.5 ){
                if( i== 1 ) j = j + 1;
                double NeutralPower = (14.7*echarge - 2.0*Kb*DustTemp)*Density*sqrt(Kb*Ti/(2.0*PI*Mp));
                std::cout << Density << "\t" << DustTemp*ConvertKelvsToeV << "\t" << Ti*ConvertKelvsToeV << "\t" << NeutralPower << "\n";
            }
        }
    }

    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
    
    std::cout << "\n\n*****\n\nNeutralRecomb UnitTest completed in " << elapsd_secs << "s\n" << i << " : " << j << "\n";
}
