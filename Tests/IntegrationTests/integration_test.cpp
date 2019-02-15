#include <ctime>            // Time program
#include <sstream>          // std::stringstream
#include <string>           // std::string

#include "ChargeTest.h" // CHARGING TESTS
#include "ForceTest.h"  // FORCE TESTS
#include "HeatTest.h"   // HEATING TESTS

static void show_usage(std::string name){
    std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
    << "\n\nOptions:\n"
    << "\t-h,--help\t\tShow this help message\n\n"
    << "\t-m,--mode MODE\tstring variable defining the Test to be run. Available Options:\n"
    << "\t\tDTOKSOML                   : Test the charging model is producing the expected result\n"
    << "\t\tGravity                    : Test constant gravitational force is reproduced by analytical model\n"
    << "\t\tEFieldGravity              : Test constant Lorentz force (Electric field with no magnetic field) and gravity\n"
    << "\t\tBField                     : Test constant Magnetic force (no Electric field with magnetic field) \n"
    << "\t\tLorentz                    : Test full Lorentz force; Electric field and magnetic field\n"
    << "\t\tLorentzGravity             : Test full Lorentz force (Electric field and magnetic field) and gravity\n"
    << "\t\tConstantHeating            : Test constant heating is comparable to analytic result\n"
    << "\t\tThermalRadiation           : Test constant heating with thermal radiation is comparable to analytic result\n"
    << "\t\tElectronHeatingRadiation   : Test constant heating with electron heating and thermal radiation\n"
    << "\t\tPlasmaHeatingRadiation     : Test constant heating with plasma heating and thermal radiation\n"
    << "\t\tPlasmaHeatingNeutralRecomb : Test constant heating with plasma heating, thermal radiation and neutral recombination\n\n";
}

template<typename T> int InputFunction(int &argc, char* argv[], int &i, std::stringstream &ss0, T &Temp){
    if (i + 1 < argc) { // Make sure we aren't at the end of argv!
        i+=1;
        ss0 << argv[i]; // Increment 'i' so we don't get the argument as the next argv[i].
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
    // Determine user input for testing mode
    std::string Test_Mode("");
    double accuracy(0.5);
    std::vector <std::string> sources;
    std::stringstream ss0;
    for (int i = 1; i < argc; ++i){ // Read command line input
        std::string arg = argv[i];
        if     ( arg == "--help"    || arg == "-h" ){   show_usage( argv[0]); return 0;         }
        else if( arg == "--mode"    || arg == "-m" )    InputFunction(argc,argv,i,ss0,Test_Mode);
        else if( arg == "--accuracy"|| arg == "-a" )    InputFunction(argc,argv,i,ss0,accuracy);
        else{
            sources.push_back(argv[i]);
        }
    }

    int NumOfElements = 4;
    char Element[NumOfElements]={'W','B','F','G'}; // Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron
    std::cout.precision(10);
    int out(0);
    for(int i(0); i < NumOfElements; i ++){


        std::cout << "\nRunning " << Test_Mode << " Test.\nTesting Element[" 
            << i << "]=" << Element[i] << "\nAccuracy = " << accuracy << "\n";
        // Charging Test
        if( Test_Mode == "DTOKSOML" ){
            // Test bout 1,
            // Integration test one for all materials
            // Test if the charging model is producing the expected result
            // This test is very rudimentary and just checks that the charging model matches DTOKS's original format
            out = ChargeTest(Element[i]);   
        }else if( Test_Mode == "Gravity" ){
            // Test bout 2,
            // Integration test two for all materials
            // Test if constant gravitational force is reproduced by analytical model
            // Produce expected results over 100 or so steps
//          out = ConstantGravityTest(Element[i]);
            out = ForceTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "EFieldGravity" ){
            // Test bout 3,
            // Integration test three for all materials
            // Test if constant Lorentz force (Electric field with no magnetic field) and gravity
            // Produce expected results over 100 or so steps
            out = ForceTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "BField" ){
            // Test bout 4, 
            // Integration test four for all materials
            // Test if a constant magnetic field force
            // Produces expected results over 100 or so steps
            out = ForceTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "Lorentz" ){
            // Test bout 5, 
            // Integration test four for all materials
            // Test if full Lorentz force (Electric field and magnetic field) and gravity
            // Produce expected results over 100 or so steps
            out = ForceTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "LorentzGravity" ){
            // Test bout 6,
            // Integration test four for all materials
            // Test if full Lorentz force (Electric field and magnetic field) and gravity
            // Produce expected results over 100 or so steps
            out = ForceTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "ConstantHeating" ){
            // Test bout 7,
            // Integration test one for all materials
            // Test if constant heating is comparable to analytic result
            out = HeatTest(Element[i],Test_Mode,accuracy);  
        }else if( Test_Mode == "ThermalRadiation" ){
            // Test bout 8, NOTE THIS DEVIATES AND I DON'T KNOW WHY EXACTLY! 
            // Integration test two for all materials
            // Test if constant heating with thermal radiation is comparable to analytic result
            out = HeatTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "ElectronHeatingRadiation" ){
            // Test bout 9, NOTE THIS DEVIATES AND I DON'T KNOW WHY EXACTLY! 
            // Integration test four for all materials
            // Test if constant heating with plasma heating and thermal radiation is comparable to analytic result
            // Note, emissivity is once again a const
            out = HeatTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "PlasmaHeatingRadiation" ){
            // Test bout 10, NOTE THIS DEVIATES AND I DON'T KNOW WHY EXACTLY! 
            // Integration test four for all materials
            // Test if constant heating with plasma heating and thermal radiation
            // Note, emissivity is once again a const
            out = HeatTest(Element[i],Test_Mode,accuracy);
        }else if( Test_Mode == "PlasmaHeatingNeutralRecomb" ){
            // Test bout 10, NOTE THIS DEVIATES AND I DON'T KNOW WHY EXACTLY! SAME PROBLEM
            // Integration test four for all materials
            // Test if constant heating with plasma heating, thermal radiation and neutral recombination
            // Note, emissivity is once again a const
            out = HeatTest(Element[i],Test_Mode,accuracy);
        }else{
            std::cout << "\n\nInput not recognised! Exiting program.\n";
            break;
        }
        std::cout << "\nFinished Running Test " << i+1 << "/" << NumOfElements << "!\n";    
        
        // Other Heating terms not tested:
        // TEE has no analytical solution as the solution involves exponential integrals
        // SEE has no result in standard mathematical functions
        // Qvap has extremely complicated analytical solutions involving the imaginary error function due to Antoinne eq.
        // Since ion flux depends on electron yield (i.e TEE and SEE), no general analytical solution for negative grains
        // For the case of constant potential, Q_(i) is constant


        if( out == 1 ) std::cout << "\n# PASSED!";
        if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
        if( out == -1 ) std::cout << "\n# FAILED!";

        std::cout << "\n\n*****\n\n\n"; 
    }
    return 0;
}
