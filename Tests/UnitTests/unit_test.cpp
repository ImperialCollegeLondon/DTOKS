#include <sstream>
#include <vector>
#include <string>

// FUNCTIONS TESTS
#include "BackscatterTest.h"
#include "DeltaSecTest.h"
#include "DeltaThermTest.h"
#include "MaxwellianTest.h"

// HEATING TESTS
#include "EvaporativeCoolingTest.h"
#include "EvaporativeMassLossTest.h"
#include "NeutralHeatingTest.h"

// CHARGING TESTS
#include "ChargingTimescales.h"
#include "DTOKSchargingTest.h"
#include "DTOKSwellchargingTest.h"
#include "SchottkyOMLTest.h"
#include "OMLTest.h"
#include "MOMLTest.h"
#include "MOMLWEMTest.h"
#include "MOMLEMTest.h"
#include "SOMLTest.h"
#include "SMOMLTest.h"
#include "PHLTest.h"
#include "THTest.h"
#include "BIBHASTest.h"

// FORCE TESTS
#include "HybridIonDrag.h"
#include "FortovIonDrag.h"
#include "NeutralDrag.h"

static void show_usage(std::string name){
    std::cerr << "Usage: int main(int argc, char* argv[]) <option(s)> SOURCES"
    << "\n\nOptions:\n"
    << "\t-h,--help\t\tShow this help message\n\n"
    << "\t-v,--variable\t\tThe number of variables to loop over; Ti/Te, Mi/Me, "
    << "Z, (U,lambda,Beta_i)\n\n"
    << "\t-m,--mode MODE\tstring variable defining the Test to be run. Availabl"
    << "e Options:\n"
    << "\t\tBackscatter    : fraction of backscattered energy and back scattere"
    << "d particles\n"
    << "\t\tDeltaSec       : empirical function calculating the yield due to se"
    << "condary electron emission\n"
    << "\t\tDeltaTherm     : value of the 'effective yield' from the Richardson"
    << "-Dushmann formula\n"
    << "\t\tMaxwellian     : value of the Maxwellian function for different val"
    << "ues of temperature and energy\n"
    << "\t\tEvapCooling    : heat loss due to evaporation\n"
    << "\t\tEvapMassLoss   : mass loss due to evaporation\n"
    << "\t\tNeutralHeating : heat gained from neutral collisions\n"
    << "\t\tDTOKScharging  : output the potential as calculated by the DTOK"
    << "S solution to the OML equation\n"
    << "\t\tDTOKSWell      : potential as calculated by the DTOKS solution "
    << "to the OML equation with a well\n"
    << "\t\tOML            : floating potential for small dust grains in a "
    << "stationary plasma following OML theory\n"
    << "\t\tMOML           : floating potential for large dust grains in a "
    << "stationary plasma following MOML theory\n"
    << "\t\tMOMLWEM        : floating potential for large emitting dust gra"
    << "ins in a stationary plasma \n"
    << "\t\tMOMLEM         : same as previous but as calculated by MOML-EM,"
    << "\n\t\t\tsee N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018)\n"
    << "\t\tSOML           : floating potential for small dust grains in a "
    << "flowing plasma following SOML theory.\n"
    << "\t\tSMOML          : floating potential for large dust grains in a "
    << "flowing plasma following SMOML theory.\n"
    << "\t\tSchottkyOML    : solution for the electron emission with Schott"
    << "ky correction\n"
    << "\t\tPHL            : floating potential for small dust grains in co"
    << "llisionless weakly magnetised plasmas.\n"
    << "\t\tTH             : floating potential for small dust grains magne"
    << "tised plasmas semi-empirical from pot.\n"
    << "\t\tBIBHAS         : floating potential for arbitary sized dust gra"
    << "in.\n"
    << "\t\tHybridIonDrag  : magnitude of the HybridIonDrag force, see http"
    << "s://doi.org/10.1063/1.1867995\n"
    << "\t\tFortovIonDrag  : magnitude of the ion drag force, see https://d"
    << "oi.org/10.1016/j.physrep.2005.08.007\n"
    << "\t\tNeutralDrag    : magnitude of the neutral drag force\n\n";

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

    // Determine user input for testing mode
    std::string Test_Mode("");
    unsigned int VariableNum(1);
    std::vector <std::string> sources;
    std::stringstream ss0;
    for (int i = 1; i < argc; ++i){ // Read command line input
        std::string arg = argv[i];
        if     ( arg == "--help"    || arg == "-h" ){   
            show_usage( argv[0]); return 0;         
        }else if( arg == "--variable"|| arg == "-v" )    
            InputFunction(argc,argv,i,ss0,VariableNum);
        else if( arg == "--mode"    || arg == "-m" )    
            InputFunction(argc,argv,i,ss0,Test_Mode);
        else{
            sources.push_back(argv[i]);
        }
    }
    assert(VariableNum > 0 && VariableNum <=6);

    // *****    FUNCTIONS TESTS     ***** //
    // Backscatter Unit Test:
    // This test prints the values of the fraction of backscattered energy and 
    // the fraction of back scattered particles
    // as calculated by the backscatter function from DTOKS. This can be 
    // readily compared to the results published in
    // the DTOKS papers
    if( Test_Mode == "Backscatter" )
        BackscatterTest();

    // Delta Sec Unit Test:
    // This test prints the value of the empirical function calculating the 
    // yield due to secondary electron emission.
    // The result can be readily compared to the publication
    // Replicating work of "Dust in tokamaks: An overview of the physical model
    // of the dust in tokamaks code"
    // Bacharis, Minas Coppins, Michael Allen, John E.
    // Page 2 & 3
    else if( Test_Mode == "DeltaSec" )
        DeltaSecTest();

    // Delta Therm Unit Test:
    // This test prints the value of the 'effective yield' from the 
    // Richardson-Dushmann formula (Without Schottky correction).
    // This can also be readily compared to the publication
    // Replicating work of "Dust in tokamaks: An overview of the physical model
    // of the dust in tokamaks code"
    // Bacharis, Minas Coppins, Michael Allen, John E.
    // Page 2 & 3
    else if( Test_Mode == "DeltaTherm" )
    DeltaThermTest();

    // Maxwellian Unit Test:
    // This test prints the value of the Maxwellian function for different 
    // values of temperature and energy.
    // The results are plotted in 3D, with the expected maxwellian distribution
    // recovered for a fixed temperature or Energy
    else if( Test_Mode == "Maxwellian" )
        MaxwellianTest();

    // *****    HEATING TESTS       ***** //
    else if( Test_Mode == "EvapCooling" )
        EvaporativeCoolingTest();
    else if( Test_Mode == "EvapMassLoss" )
        EvaporativeMassLossTest();
    else if( Test_Mode == "NeutralHeating" )
        NeutralHeatingTest();


    // *****    CHARGING TESTS      ***** //
    // Charging Timescale Test:
    // This test is used to verify that the timestep as calculated by 
    // Krasheninnikovs is always smaller
    // than the electron plasma frequency. In practice, this is found to not be
    // perfectly true, there exist
    // extreme conditions where this is not the case.
    // Time step based on the formulation by Krasheninnikov
    // Smirnov, R. D., Pigarov, A. Y., Rosenberg, M., Krasheninnikov, S. I., 
    // & Mendis, D. a. (2007). 
    // Modelling of dynamics and transport of carbon dust particles in tokamaks. 
    // Plasma Physics and Controlled Fusion, 49(4), 347–371.    
//      else if( Test_Mode == "ChargingTime" )
//      ChargingTimescales();
    
    // DTOKS Charging Test:
    // This test output the potential as calculated by the DTOKS solution to
    // the OML equation.
    // The form of this is such that it depends on the sign of the potential
    // and the magnitude of the total electron emission.
    // This test was used to show the discontinuity in the potential when the
    // total electron emission yield approaches 1
    // Two different conditional formulations of the problem are made and their
    // differences highlighted by this test
    else if( Test_Mode == "DTOKScharging" )
        DTOKSchargingTest(VariableNum);

    // This test outputs the potential as calculated by part of the DTOKS
    // solution to the OML equation, considering only a 
    // scenario with a potential well. This removes the discontinuities
    // observed in the previous model.
    else if( Test_Mode == "DTOKSWell" )
        DTOKSwellchargingTest();
    
    // OML Charging Test:
    // This test is designed to find the floating potential for small dust
    // grains in a stationary plasma following OML theory.
    // This employs an approximate series expansion to the Lambert W function
    // to find the floating potential
    else if( Test_Mode == "OML" )
        OMLTest(VariableNum);

    // MOML Charging Test:
    // This test is designed to find the floating potential for large dust
    // grains in a stationary plasma following MOML theory.
    // This employs an approximate series expansion to the Lambert W function
    // to find the floating potential
    else if( Test_Mode == "MOML" )
        MOMLTest(VariableNum);

    // MOMLwEM Charging Test:
    // This test is designed to find the floating potential for large emitting
    // dust grains in a stationary plasma 
    // following MOML With electron emission theory, see: N. Rizopoulou and M.
    // Bacharis, Phys. Plasmas 25, (2018). 
    // Equation (1) and (2)
    // This employs an approximate series expansion to the Lambert W function
    // to find the floating potential
    else if( Test_Mode == "MOMLWEM" )
        MOMLWEMTest(VariableNum);

    // MOML-EM Charging Test:
    // This test is designed to find the floating potential for large emitting
    // dust grains in a stationary plasma 
    // following MOML-EM theory, see: N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018). &
    // N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
    // This employs a Newton Rhapson method to solve the equations of the two papers finding a kinetic potential well
    else if( Test_Mode == "MOMLEM" )
        MOMLEMTest();


    // SOML Charging Test:
    // This test is designed to find the floating potential for small dust grains in a flowing plasma following SOML theory.
    // This employs an approximate series expansion to the Lambert W function to find the floating potential
    else if( Test_Mode == "SOML" )
        SOMLTest(VariableNum);

    // SMOML Charging Test:
    // This test is designed to find the floating potential for large dust grains in a flowing plasma following SMOML theory.
    // This employs an approximate series expansion to the Lambert W function to find the floating potential
    else if( Test_Mode == "SMOML" )
        SMOMLTest(VariableNum);

    // Schottky OML Charging Test:
    // This test was made to see what the solution is for electron emission with Schottky correction where the potential of the
    // dust grain is accounted for. The minimisation of the positive solution in C++ is not stable and gives an incorrect 
    // answer. Switching to matlab minimisation function, some weird things happen but, in principle, I showed that the 
    // function could be minimised.
    else if( Test_Mode == "SchottkyOML" )
        SchottkyOMLTest();

    // SchottkyMOML Charging Test: DOESN'T WORK!
    // This test is designed to find the floating potential for large negative dust grains with electron emission 
    // This employs an approximate series expansion to the Lambert W function to find the floating potential
//      else if( Test_Mode == "SchottkyMOML" )
//      SchottkyMOMLTest();

    // PHL Charging Test:
    // This test is designed to find the floating potential for small dust grains in collisionless weakly magnetised plasmas.
    // The calculation follows the work by Patacchini et al and implements a semi-empirical model
    // see L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007). for details
    else if( Test_Mode == "PHL" )
        PHLTest();

    // TH Charging Test:
    // This test is designed to find the floating potential for small dust grains in magnetised plasmas 
    // The calculation follows the work by Drew Thomas and Josh Holgate and implements a semi-empirical model
    // see D. M. Thomas and J. T. Holgate, ArXiv Prepr. (2016). for details
    else if( Test_Mode == "TH" )
        THTest();

    // BIBHAS Charging Test:
    // This test is designed to find the floating potential for arbitary sized dust grain
    // The calculation follows the work by R. DE Bibhas,
    // see R. DE Bibhas, Astrophys. Space Sci. 30, (1974).
    else if( Test_Mode == "BIBHAS" )
        BIBHASTest();




    // *****    FORCE TESTS     ***** //
    // Hybrid Ion Drag Test
    // This test is designed to test the magnitude of the HybridIonDrag force as formulated in the paper given below.
    // The Hybrid Ion Drag model is a function of the plasma species temperature ratio, the ion mach number, ion density,
    // the dust grain potential and the dust grain radius
    // Khrapak, S. A., Ivlev, A. V., Zhdanov, S. K., & Morfill, G. E. (2005). Hybrid approach to the ion drag force. 
    // Physics of Plasmas, 12(4), 1–8. https://doi.org/10.1063/1.1867995 
    else if( Test_Mode == "HybridIonDrag" )
        HybridIonDragTest();

    // Fortov et al./DTOKS Ion Drag Test
    // This test is designed to test the magnitude of the drag force as formulated in the paper given below by fortov.
    // The Hybrid Ion Drag model is a function of the temperature of electrons and ions, the ion mach number, ion density,
    // the dust grain potential and the dust grain radius
    // Fortov, V. E., Ivlev, A. V., Khrapak, S. A., Khrapak, A. G., & Morfill, G. E. (2005).
    //  Complex (dusty) plasmas: Current status, open issues, perspectives. Physics Reports, 421(1–2), 1–103. 
    // https://doi.org/10.1016/j.physrep.2005.08.007
    else if( Test_Mode == "FortovIonDrag" )
        FortovIonDragTest();
    // Neutral Drag test
    // This test is designed to test the magnitude of the neutral drag force relative to the Ion drag force.
    // This neutral drag force is formulated by the OML flux for uncharged species to a sphere
    else if( Test_Mode == "NeutralDrag" )
        IonNeutralDragTest();
    else
        std::cout << "\n\nInput not recognised! Exiting program\n.";

    return 0;
}
