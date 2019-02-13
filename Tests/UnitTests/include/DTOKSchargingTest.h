#include "Constants.h"
#include "DTOKSsolveOML.h"
#include <iostream>

// This test output the potential as calculated by the DTOKS solution to the OML
// equation.
// The form of this is such that it depends on the sign of the potential and the
// magnitude of the total electron emission.
// This test was used to show the discontinuity in the potential when the total
// electron emission yield approaches 1
// Two different conditional formulations of the problem are made and their
// differences highlighted by this test
void DTOKSchargingTest(double Variables){
    clock_t begin = clock();
    std::cout << "\nOMLTest\t\tVariables : " << Variables;
    std::cout << "\n\nTiTe\tMi/Me\tZ\tTd\tPotential";

    double converteVtoK(11600);

    // TEST TO CALCULATE THE DTOKS FLOATING POTENTIAL FOR 
    // CONSTANT ELECTRON DENSITY AND ELECTRON TEMPERATURE

    // Electron Temp in ev DO NOT CHANGE THIS VALUE. CHECK DTOKSsolveOML.h FIRST
    double Te = 1;      
    double Td = 300;    // Dust Temp in K
    double ne = 1e18;   // Electron Density in m^-3
    double Potential = 2;   // Normalised potential
    double Z = 1.0;     // Mean Ionization
    
    double Sec = sec(Td/converteVtoK,'w');
    double gammae = ne*exp(Potential)*sqrt(echarge*Te/(2*PI*Me));
    double Therm = Richardson*pow(Td,2)*exp(-(4.55*echarge)/(Kb*Td))/
        (echarge*gammae);
    if( (Sec + Therm) >= 1.0 ){ 
        Potential = DTOKSsolveOML( 0.0, 0.01, Te, Potential) -
            Td*Kb/(echarge*Te);
    }else{ // If the grain is negative...
        Potential = DTOKSsolveOML( Sec + Therm, 0.01, Te, Potential);
        if( Potential < 0.0 ){ // But if it's positive
            // But! If it's now positive, our assumptions must be wrong!
            // So now we assume it's positive and calculate the potential with 
            // a well.
            Potential = DTOKSsolveOML(0.0, 0.01, Te, Potential)-
                Td*Kb/(echarge*Te);
        }
    }

    for( double TiTe(0.01); TiTe < 100; TiTe *= 1.1 ){
        for( double Td(300); Td < (Variables-1)*5500; Td += 10){  
        
            double Sec = sec(Td/converteVtoK,'w');
            double gammae = ne*exp(Potential)*sqrt(echarge*Te/(2*PI*Me));
            double Therm = Richardson*pow(Td,2)*exp(-(4.55*echarge)/(Kb*Td))/
                (echarge*gammae);
            // BEGINING STRUCTURE ONE : NEW CODE
//          if( (Sec + Therm) < 1.0 ){ // DTOKSsolveOML only defined for deltatot < 1.0
//              Potential = DTOKSsolveOML( Sec + Therm, TiTe, Te, Potential); 
//          }else{ // If the grain is in fact positive ...
//              Potential = DTOKSsolveOML( Sec + Therm, TiTe, Te, Potential);
//              if( Potential < 0.0 ){
//                  Potential = DTOKSsolveOML(0.0, TiTe, Te, Potential)-Kb*Td/(echarge*Te);
//              }
//          }

            // BEGINNING STRUCTURE TWO : OLD CODE
            if( (Sec + Therm) >= 1.0 ){ 
                Potential = DTOKSsolveOML( 0.0, TiTe, Te, Potential) - 
                    Td*Kb/(echarge*Te);
            }else{ // If the grain is negative...
                Potential = DTOKSsolveOML( Sec + Therm, TiTe, Te, Potential);
                if( Potential < 0.0 ){ // But if it's positive
                    // But! If it's now positive, our assumptions must be wrong!
                    // So now we assume it's positive and calculate the 
                    // potential with a well.
                    Potential = DTOKSsolveOML(0.0, TiTe, Te, Potential)-
                    Td*Kb/(echarge*Te);
                }
            }

            std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << Z << "\t" 
                << Td << "\t" << (Sec+Therm) << "\t" << Potential;
        }
        Potential = DTOKSsolveOML( Sec + Therm, TiTe, Te, Potential);
        if( Potential < 0.0 ){ // But if it's positive
            // But! If it's now positive, our assumptions must be wrong!
            // So now we assume it's positive and calculate the potential with 
            // a well.
            Potential = DTOKSsolveOML(0.0, TiTe, Te, Potential)-
                Td*Kb/(echarge*Te);
        }
        std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << Z << "\t" 
                << Td << "\t" << (Sec+Therm) << "\t" << Potential;
    }

    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nDTOKSCharging UnitTest completed in " << elapsd_secs << "s\n";
}
