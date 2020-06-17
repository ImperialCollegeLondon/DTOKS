#include "solveMOMLEM.h"
#include <iostream>
#include <limits>

// Test orbital motion limited theory
void MOMLEMTest(unsigned int Variables){
    clock_t begin = clock();

    std::cout << "\n#MOMLEMTest\t\tVariables : " << Variables;
    /* CASE 1
    REPLICATING FIGURE 3.
    Finding Critical Emmitted temperature ratio with fixed:
        Emitted density ratio         = 0.3
        Temperature Ratio         = 1.0
    */
    double guess = -1.0;    // Initial guess at potential
    double EmittedDensityRatio = 0.3;
    double TiTe = 1.0;
//    for( double TiTe(0.01); TiTe < 100; TiTe *= 1.5 ){

      std::cout.precision(std::numeric_limits<double>::max_digits10);
    std::cout << "\n* CRITICAL EMITTED TEMPERATURE RATIO *";
    double CritVal = FindCriticalVal(TiTe,EmittedDensityRatio,Mp/Me,'d');
    std::cout << "\n" << CritVal;

    double EmittedTemperatureRatio = CritVal;
    exploreWellCase(TiTe,0.2,EmittedTemperatureRatio,1.0/(Mp/Me));
    for( EmittedTemperatureRatio = CritVal; EmittedTemperatureRatio < 0.201; EmittedTemperatureRatio += 0.01 ){
//        double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio);

//        std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
//            << "\t" << EmittedTemperatureRatio << "\t" << Potential 
//            << "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
    }

    /* CASE 2
    REPLICATING FIGURE 2.
    Finding Critical Emitted Density Ratio with fixed:
        Emitted temperature ratio     = 0.5
        Temperature Ratio         = 1.0
    */
    TiTe = 1.0;

//    std::cout << "\n* CRITICAL EMITTED DENSITY RATIO *";
//    CritVal = FindCriticalVal(TiTe,EmittedTemperatureRatio,Mp/Me,'c');
//    std::cout << "\n" << CritVal;
/*
    for( TiTe = 0.5; TiTe < 7.0; TiTe += 0.1 ){
        double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio);
//        std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
//            << "\t" << EmittedTemperatureRatio << "\t" << guess << "\t" << Potential
//            << "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
//        std::cout << "\t" << Potential;
    }
*/
    /* CASE 3
    REPLICATING FIGURE 4.
    Finding Critical Temperature Ratio with fixed:
        Emitted density ratio         = 0.3
        EittedTemperature Ratio     = 0.5
    */
    EmittedTemperatureRatio = 0.5;

//    std::cout << "\n* CRITICAL TEMPERATURE RATIO *";
//    CritVal = FindCriticalVal(EmittedDensityRatio,EmittedTemperatureRatio,Mp/Me,'t');
//    std::cout << "\n" << CritVal;
/*
    for( EmittedDensityRatio = 0.01; EmittedDensityRatio < 0.5; EmittedDensityRatio += 0.01 ){
        double Potential = solveMOMLEM(TiTe,Mp/Me,EmittedDensityRatio,EmittedTemperatureRatio);
//        std::cout << "\n" << TiTe << "\t" << Mp/Me << "\t" << EmittedDensityRatio 
//            << "\t" << EmittedTemperatureRatio << "\t" << guess << "\t" << Potential
//            << "\t" << log(sqrt(TiTe/(Mp/Me))+EmittedDensityRatio*sqrt(EmittedTemperatureRatio));
//        std::cout << "\t" << Potential;
    }
*/
    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
//    std::cout << "\n\n*****\n\nMOMLEM UnitTest completed in " << elapsd_secs << "s\n";
}
