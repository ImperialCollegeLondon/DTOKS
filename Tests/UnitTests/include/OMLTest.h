#include "solveOML.h"
#include <iostream>

// Test orbital motion limited theory
void OMLTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\nOMLTest\t\tVariables : " << Variables;
    std::cout << "\n\nTiTe\tMi/Me\tZ\tPotential";

    double Mi=Mp;
    double Z=1.0;
    for( double TiTe(0.001); TiTe < 200; TiTe *= 2 ){

        for( Mi=Mp; Mi < (Variables-1)*10*Mp; Mi *= 2 ){
            for( Z=1.0; Z < (Variables-2)*50; Z ++ ){
                double Potential = solveOML(TiTe,Mi/Me,Z);
                std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z  << "\t"
                    << Potential;
            }
            double Potential = solveOML(TiTe,Mi/Me,Z);
            std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z  << "\t"
                    << Potential;
        }
        double Potential = solveOML(TiTe,Mi/Me,Z);
        std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z  << "\t"
                    << Potential;
    }
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nOML UnitTest completed in " << elapsd_secs 
        << "s\n";
}
