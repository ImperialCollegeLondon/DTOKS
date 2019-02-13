#include "solveSOML.h"
#include <iostream>

// Test Shifted orbital motion limited theory
void SOMLTest(double Variables){
    clock_t begin = clock();
    std::cout << "\nSOMLTest\t\tVariables : " << Variables;
    std::cout << "\n\nTiTe\tMi/Me\tZ\tU\tPotential";

    double U = 1.0; // Plasma Flow speed
    double TiTe=1.0;
    double Mi=Mp;
    double Z=1.0;
    // Loop over temperature ratios
    for( double TiTe(0.01); TiTe < 10; TiTe *= 1.5 ){ 
        for( Mi=Mp; Mi < (Variables-1)*10*Mp; Mi *= 2 ){
            for( Z=1.0; Z < (Variables-2)*50; Z ++ ){
                // Loop over plasma flow speed
                for( double U(0.1); U < (Variables-3)*5; U *= 1.2 ){ 
                    double Potential = solveSOML(TiTe,Mi/Me,Z,U*sqrt(2.0));
                    std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z
                        << "\t" << U << "\t" << Potential;
                }
                double Potential = solveSOML(TiTe,Mi/Me,Z,U*sqrt(2.0));
                std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z
                    << "\t" << U << "\t" << Potential;
            }
            double Potential = solveSOML(TiTe,Mi/Me,Z,U*sqrt(2.0));
            std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z
                << "\t" << U << "\t" << Potential;
        }
        double Potential = solveSOML(TiTe,Mi/Me,Z,U*sqrt(2.0));
        std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z
            << "\t" << U << "\t" << Potential;
    }
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nSOML UnitTest completed in " << elapsd_secs 
        << "s\n";
}
