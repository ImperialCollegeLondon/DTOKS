#include "Constants.h"
#include "DTOKSsolveOML.h"
#include <iostream>

// This test outputs the potential as calculated by part of the DTOKS solution to the OML equation, considering a potential well.
void DTOKSwellchargingTest(unsigned int Variables){
    clock_t begin = clock();
    std::cout << "\n#DTOKSwellcharging\t\tVariables : " << Variables;
    double converteVtoK(11600);

    // TEST TO CALCULATE THE DTOKS FLOATING POTENTIAL FOR CONSTANT ELECTRON DENSITY AND ELECTRON TEMPERATURE

    double Te = 1;         // Electron Temp in ev DO NOT CHANGE THIS VALUE. CHECK DTOKSsolveOML.h FIRST
    double Td = 300;     // Dust Temp in K
    double ne = 1e18;     // Electron Density in m^-3
    
    double Potential = DTOKSsolveOML(0.0, 0.01, Te, Potential)-Td*Kb/(echarge*Te);

    for( double Td(300); Td < 5500; Td += 10){    
        for( double TiTe(0.01); TiTe < 100; TiTe *= 1.1 ){
            double Sec = sec(Td/converteVtoK,'w');
            double gammae = ne*exp(Potential)*sqrt(echarge*Te/(2*PI*Me));
            double Therm = Richardson*pow(Td,2)*exp(-(4.55*echarge)/(Kb*Td))/(echarge*gammae);

            Potential = DTOKSsolveOML(0.0, TiTe, Te, Potential)-Td*Kb/(echarge*Te);

            std::cout << "\n" << Td << "\t" << TiTe << "\t" << (Sec+Therm) << "\t" << Potential;
        }
    }

    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nDTOKSwellcharging UnitTest completed in " << elapsd_secs << "s\n";
}
