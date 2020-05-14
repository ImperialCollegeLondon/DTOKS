#include "solveMOMLWEM.h"
#include <iostream>

// Test orbital motion limited theory
void MOMLWEMTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\n#MOMLWEMTest\t\tVariables : " << Variables;
    std::cout << "\n\nTiTe\tMi/Me\tZ\tDeltaTot\tPotential";
    
    double HeatCapacity = 5.0/3.0; // Ionization state of plasma
    double Z = 1.0; // Ionization state of plasma
    double Mi=Mp;
    double DeltaTot=0.01;
//  for( double TiTe(0.01); TiTe < 100; TiTe *= 1.5 ){
    for( double TiTe(0.01); TiTe < 1.4; TiTe *= 1.5 ){
        for( Mi=Mp; Mi < (Variables-1)*10*Mp; Mi *= 2 ){
            for( Z=1.0; Z < (Variables-2)*50; Z ++ ){
                for( DeltaTot=0.01; DeltaTot < (Variables-3)*20; 
                    DeltaTot *= 1.1 ){

                    double Potential = solveMOMLWEM(TiTe,Mi/Me,Z,
                        HeatCapacity,DeltaTot);
                    std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" 
                        << Z << "\t" << HeatCapacity << "\t" 
                        << DeltaTot << "\t" << Potential;
                }
                double Potential = solveMOMLWEM(TiTe,Mi/Me,Z,
                    HeatCapacity,DeltaTot);
                std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" 
                    << Z << "\t" << HeatCapacity << "\t" 
                    << DeltaTot << "\t" << Potential;
            }
            double Potential = solveMOMLWEM(TiTe,Mi/Me,Z,
                HeatCapacity,DeltaTot);
            std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" 
                << Z << "\t" << HeatCapacity << "\t" 
                << DeltaTot << "\t" << Potential;
        }
        double Potential = solveMOMLWEM(TiTe,Mi/Me,Z,
            HeatCapacity,DeltaTot);
        std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" 
            << Z << "\t" << HeatCapacity << "\t" 
            << DeltaTot << "\t" << Potential;
    }
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nMOMLWEM UnitTest completed in " << elapsd_secs << "s\n";
}
