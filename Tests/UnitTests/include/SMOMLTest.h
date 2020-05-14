#include "solveSMOML.h"
#include <iostream>

// Test Shifted modified orbital motion limited theory
void SMOMLTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\n#SMOMLTest\t\tVariables : " << Variables;
    std::cout << "\n\nTiTe\tMi/Me\tZ\tU\tPotential";

    double Mi=Mp;
    double Z= 1.0;
    double U = 1.0; // Plasma Flow speed
    double Gamma = 3.0; // Heat Capacity Ratio
    // Loop over temperature ratios
    for( double TiTe(0.01); TiTe < 50; TiTe *= 1.5 ){ 
        for( Mi=Mp; Mi < (Variables-1)*10*Mp; Mi *= 2 ){
            for( Z=1.0; Z < (Variables-2)*50; Z ++ ){
                // Loop over plasma flow speed
                for( double U(0.1); U < (Variables-3)*3; U *= 1.5 ){ 
                    double Potential = solveSMOML(TiTe,Mi/Me,Z,Gamma,U);
                    std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
                        << "\t" << U << "\t" << Potential;
                }
                double Potential = solveSMOML(TiTe,Mi/Me,Z,Gamma,U);
                std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
                    << "\t" << U << "\t" << Potential;
            }
            double Potential = solveSMOML(TiTe,Mi/Me,Z,Gamma,U);
            std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
                << "\t" << U << "\t" << Potential;
        }
        double Potential = solveSMOML(TiTe,Mi/Me,Z,Gamma,U);
        std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
            << "\t" << U << "\t" << Potential;
    }
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nSMOML UnitTest completed in " << elapsd_secs 
        << "s\n";
}
