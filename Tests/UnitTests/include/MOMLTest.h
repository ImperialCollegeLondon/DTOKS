#include "solveMOML.h"
#include <iostream>

// Test Modified orbital motion limited theory
void MOMLTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\n#MOMLTest\t\tVariables : " << Variables;
    std::cout << "\n\nTiTe\tMi/Me\tZ\tGamma\tPotential";

//  double Current0 = MOMLattractedCurrent(1.0,1.0,Gamma,Mp/Me,0.0);
    double Gamma = 1.0; // Ionization state of plasma
    double Mi=Mp;
    double Z=1.0;
    for( double TiTe(0.001); TiTe < 200; TiTe *= 2 ){
        for( Mi=Mp; Mi < (Variables-1)*10*Mp; Mi *= 2 ){
            for( Z=1.0; Z < (Variables-2)*50; Z ++ ){
//              double Current = MOMLattractedCurrent(1.0,1.0,Gamma,Mp/Me,TiTe);
                for( double Gamma(0.01); (Variables-3)*Gamma < 100; 
                    Gamma *= 2 ){
                    double Potential = solveMOML(TiTe,Mi/Me,Z,Gamma);
                    std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
                        << "\t" << Gamma << "\t" << Potential;
                }
                double Potential = solveMOML(TiTe,Mi/Me,Z,Gamma);
                std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
                        << "\t" << Gamma << "\t" << Potential;
            }
            double Potential = solveMOML(TiTe,Mi/Me,Z,Gamma);
            std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
                        << "\t" << Gamma << "\t" << Potential;
        }
        double Potential = solveMOML(TiTe,Mi/Me,Z,Gamma);
        std::cout << "\n" << TiTe << "\t" << Mi/Me << "\t" << Z 
                        << "\t" << Gamma << "\t" << Potential;
    }
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nMOML UnitTest completed in " << elapsd_secs 
        << "s\n";
}
