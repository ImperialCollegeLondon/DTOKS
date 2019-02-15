/** @file Beryllium.cpp
 *  @brief Tests of individual component models of DTOKSU
 *
 * This directory contains Tests of individual component models of DTOKSU, meaning the
 * Heating Model, Charging Model and Force Model. The purpose of these tests is to compare
 * the difference in result between the original models present in DTOKS with the four
 * additional models since I began work on the code. These namely are:
 * 
 * Neutral drag, Neutral Heating, Evaporative Cooling (and Mass loss), Variable Heat capacity,
 * Variable Emissivity
 * 
 * Additionally, there is a test to compare the difference between the original DTOKS heating
 * model with and without these newer models. The following standard plasma conditions are
 * utilised throughout these tests unless otherwise stated:
 * 
 * DustSize         = 1e-6;         // m
 * Temp             = 280;          // K
 * Potential        = 1;            // Normalised Potential
 * NeutralDensity   = 1e18;         // m^-3, Neutral Density
 * IonDensity       = 1e18;         // m^-3, Ion Density
 * ElectronDensity  = 1e18;         // m^-3, Electron Density
 * NeutralTemp      = 10*1.16e4;    // K, Neutral Temperature, convert from eV
 * IonTemp          = 10*1.16e4;    // K, Ion Temperature, convert from eV
 * ElectronTemp     = 10*1.16e4;    // K, Electron Temperature, convert from eV
 * 
 * For each pair of tests, two different text files are generated with the force data within
 * them. Both are named after the name of the function. In the case of the 'neutraldragtest'
 * function, they are "Data/NeutralDragOn.txt" and "Data/NeutralDragOff.txt".

 * NOTE! These tests must be run from the /DTOKS-U/ directory since the file paths to the plasma
 * data and emissivity data files are relative not global.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */


#include <ctime>            // Time program

#include "NeutralDragTest.h"
#include "NeutralHeatingTest.h"
#include "EvaporativeCoolingTest.h"
#include "VariableHeatCapacityTest.h"
#include "VariableEmissivityTest.h"
#include "BeforeAfterHeatingTest.h"

static void show_usage(std::string name){
    std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
    << "\n\nOptions:\n"
    << "\t-h,--help\t\tShow this help message\n\n"
    << "\t-m,--mode MODE\tstring variable defining the Test. Available Options:\n"
    << "\t\tNeutralDragTest          : Test neutral drag force\n"
    << "\t\tNeutralHeatingTest       : Test neutral heating\n"
    << "\t\tEvaporativeCoolingTest   : Test effect of evaporative cooling\n"
    << "\t\tVariableHeatCapacityTest : Test impact of variable heat capacity\n"
    << "\t\tVariableEmissivityTest   : Test impact of variable emissivity \n"
    << "\t\tBeforeAfterHeatingTest   : Test impact of heating\n\n";
}

template<typename T> int InputFunction(int &argc, char* argv[], int &i, 
    std::stringstream &ss0, T &Temp){
    if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i+=1;
        ss0 << argv[i]; // Increment 'i' get the argument as the next argv[i].
        ss0 >> Temp;
        ss0.clear(); ss0.str("");
        return 0;
    }else{ // Uh-oh, there was no argument to the destination option.
        std::cerr << "\noption requires argument." << std::endl;
        return 1;
    }
}


int main(int argc, char* argv[]){
    if( argc == 1 ){
        std::cout << "\nNo Mode Selected! Program Exiting\n\n";
        return 0;
    }
    char Element[5]={'W','B','F','G','D'};// Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium, (F) : Iron or (D) : Deuterium.
    int out(0);

    // Determine user input for testing mode
    std::string Test_Mode("");
    std::vector <std::string> sources;
    std::stringstream ss0;
    for (int i = 1; i < argc; ++i){ // Read command line input
        std::string arg = argv[i];
        if     ( arg == "--help"    || arg == "-h" ){   
            show_usage( argv[0]); return 0;         
        }
        else if( arg == "--mode"    || arg == "-m" )    
            InputFunction(argc,argv,i,ss0,Test_Mode);
        else{
            sources.push_back(argv[i]);
        }
    }

    for(int i(0); i < 1; i ++){

        std::cout << "\nElement["<<i<<"] is = ";
        std::cout << Element[i];

//      Model Test 1, Neutral Drag Test:
//      This test evaluates the effect of introducing the neutral drag force into DTOKS. The force
//      model is used with an arbitrary initial position (1.0,0.0,0.0). The initial velocity of
//      the dust is set to (2.0,3.0,5.0) with the plasma velocity being set to
//      PlasmaVelocity   = (0.0,0.0,0.5*sqrt(Kb*NeutralTemp/Mp));
//      The equation for neutral drag is given by:
//      Fnd = (PlasmaVelocity-SampleVelocity)*Mp*sqrt(4*PI)*NeutralFlux*PI*pow(SampleRadius,2)/SampleMass
//      The numerical solution over 100,000 time steps is calculated. The results for with and
//      without the neutral drag force can readily be compared.
        if( Test_Mode == "NeutralDragTest" ){
            out = NeutralDragTest(Element[i],false);
            out = NeutralDragTest(Element[i],true);
        }

//      Model Test 2, Neutral Heating Test:
//      This test evaluates the effect of including Neutral heating in the heating model. The
//      neutral heating term is given by:
//      NeutralHeatFlux = SurfaceArea*2*NeutralFlux*NeutralTemperature*Kb; // Watts     
//      The heating model is until the sample vapourises or reaches thermal equilibrium.
        else if( Test_Mode == "NeutralHeatingTest" ){
            out = NeutralHeatingTest(Element[i],false);
            out = NeutralHeatingTest(Element[i],true);
        }

//      Model Test 3, Evaporative Cooling Test:
//      This test evaluates the effect of energy and mass loss within the heating model. The mass
//      loss provides an additional end condition for the simulation below the boiling temperature
//      due to vaporisation. The equation for evaporative cooling is given by the Hertz-Knudsen
//      equation which is
//      EvaporativeFlux = (StickCoeff*SampleSurfaceArea*AvNo*(SampleVapourPressure-AmbientPressure))/
//                              sqrt(2*PI*SampleAtomicMass*R*DustTemperature);
//      which gives the rate of mass loss as
//      MassLoss = EvaporationFlux*SampleAtomicMass/AvNo
//      and the rate of energy loss
//      HeatLoss = EvaporationFlux*(3*Kb*DustTemperature/2)+SampleBondEnergy/AvNo)
//      The heating model is ran until the sample vapourises or reaches thermal equilibrium.
        else if( Test_Mode == "EvaporativeCoolingTest" ){
            out = EvaporativeCoolingTest(Element[i],false);
            out = EvaporativeCoolingTest(Element[i],true);
        }

//      Model Test 4, Variable Heat Capacity:
//      This test evaluates the effect of allowing Heat capacity to vary as a function of
//      temperature. For each of the different materials, a different empirical relation has been
//      used to define the variation with temperature. Most of these relations are found at:
//      Chase, M. (1998). NIST-JANAF Thermochemical Tables, 4th Edition.
//      Journal of Physical and Chemical Reference Data, Monograph 9.
//      These are compared to the 'default' constant values of heat capacity which are set by the
//      constructors of each element where, by design, a function called '<elementname>_defaults()'
        else if( Test_Mode == "VariableHeatCapacityTest" ){
//      is called by the constructor to establish their value.
            out = VariableHeatCapacityTest(Element[i],false);
            out = VariableHeatCapacityTest(Element[i],true);
        }

//      Model Test 5, Variable Emissivity Test:
//      This test evaluates the effect of allowing Emissivity to vary as a function of
//      temperature and size. For each of the different materials, a set of tabulated data has
//      been generated from the Wiscombe Mie scattering Fortran code, see reference
//      Wiscombe, W. J. (1980). Improved Mie scattering algorithms. Applied Optics, 19(9), 1505–9.
//      This data is not available in the github repository but can be made available locally
//      if placed under the heading:
//      /Models/EmissivityData/EmissivityData<ElementName>/Temp_<Temperature(K)>.txt
//      for example ...
//      /Models/EmissivityData/EmissivityDataIron/Temp_1383.txt
//      where Temp_1383.txt contains two collum data for a temperature of 1383K of comma delimited
//      data. The first collum is the size of the dust grain in metres and the second collum is
//      the emissivity.
//      These are compared to the 'default' constant values of emissivity which are set by the
//      constructors of each element where, by design, a function called '<elementname>_defaults()'
//      is called by the constructor to establish their value.

        else if( Test_Mode == "VariableEmissivityTest" ){
            out = VariableEmissivityTest(Element[i],false);
            out = VariableEmissivityTest(Element[i],true);
        }

//      Model Test 6, Before After Heating Test:
//      This test evaluates the total effect of the introduction of all four of the new models
//      to the heating model. This includes the variational emissivity and heat capacity, plus the
//      addition of the neutral heating and evaporative cooling terms. 
//      These are compared to the 'default value' versions with the additional two terms in the heat
//      equation turned off.
        else if( Test_Mode == "BeforeAfterHeatingTest" ){
            out = BeforeAfterHeatingTest(Element[i],false);
            out = BeforeAfterHeatingTest(Element[i],true);
        }else
            std::cout << "\n\nInput not recognised! Exiting program.\n";
        std::cout << "\n\n*****\n"; 
    }
    return 0;
}
