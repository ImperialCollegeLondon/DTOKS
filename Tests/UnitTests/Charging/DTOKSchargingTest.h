#include "Constants.h"
#include "DTOKSsolveOML2.h"
#include "DTOKSsolveOML3.h"
#include <iostream>

// This test output the potential as calculated by the DTOKS solution to the OML equation.
// The form of this is such that it depends on the sign of the potential and the magnitude of the total electron emission.
// This test was used to show the discontinuity in the potential when the total electron emission yield approaches 1
// Two different conditional formulations of the problem are made and their differences highlighted by this test
void DTOKSchargingTest(){
	clock_t begin = clock();

	double converteVtoK(11600);
	// TEST TO COMPARE TWO DIFFERENT WAYS OF APPROACHING THE DTOKS OML PROBLEM

	double Potential = 2;
	double Potential2 = 2;

	for(double Td(300); Td < 5e3; Td ++){ // Loop over temperatures
		double Sec = sec(Td/converteVtoK,'w'); 
		for( double ne(1e18); ne < 1e20; ne *= 1.1 ){
			for( double Te(1); Te < 100; Te *= 1.1 ){

				double gammae = ne*exp(Potential)*sqrt(echarge*Te/(2*PI*Me));
				for( double Ti(1); Ti < 100; Ti *= 1.1 ){
					double Therm = Richardson*pow(Td,2)*exp(-(4.55*echarge)/(Kb*Td))/(echarge*gammae);

					// BEGINING STRUCTURE ONE : NEW CODE
					if( (Sec + Therm) < 1.0 ){ // DTOKSsolveOML2 only defined for deltatot < 1.0
						Potential = DTOKSsolveOML2( Sec + Therm, Ti, Te, Potential); 
					}else{ // If the grain is in fact positive ...
						Potential = DTOKSsolveOML2( Sec + Therm, Ti, Te, Potential);
						if( Potential < 0.0 ){
							Potential = DTOKSsolveOML2(0.0, Ti, Te, Potential)-Kb*Td/(echarge*Te);
						}
					}
					// BEGINNING STRUCTURE TWO : OLD CODE
					if( (Sec + Therm) >= 1.0 ){ 
						Potential2 = DTOKSsolveOML2( 0.0, Ti, Te, Potential2) - Td*Kb/(echarge*Te);
					}else{ // If the grain is negative...
						Potential2 = DTOKSsolveOML2( Sec + Therm, Ti, Te, Potential2);
						if( Potential2 < 0.0 ){ // But if it's positive
							// But! If it's now positive, our assumptions must be wrong!
							// So now we assume it's positive and calculate the potential with a well.
							Potential2 = DTOKSsolveOML2(0.0, Ti, Te, Potential2)-Td*Kb/(echarge*Te);
						}
					}
					// Proof that they're identical... They are mostly but not always
					if( Potential != Potential2 ){
						std::cout << "\n" << Td << "\t" << ne << "\t" << Ti/Te << "\t" << Potential << "\t" << Potential2;
					}
				}
			}
		}
	}


	// TEST TO CALCULATE THE DTOKS FLOATING POTENTIAL FOR CONSTANT ELECTRON DENSITY AND ELECTRON TEMPERATURE
/*
	double Te = 1; 		// Electron Temp in ev
	double Td = 300; 	// Dust Temp in K
	double ne = 1e18; 	// Electron Density in m^-3
	double Potential = 2;	// Normalised potential
	
	double Sec = sec(Td/converteVtoK,'w');
 	double gammae = ne*exp(Potential)*sqrt(echarge*Te/(2*PI*Me));
	double Therm = Richardson*pow(Td,2)*exp(-(4.55*echarge)/(Kb*Td))/(echarge*gammae);
	if( (Sec + Therm) >= 1.0 ){ 
		Potential = DTOKSsolveOML3( 0.0, 0.01, Potential) - Td*Kb/(echarge*Te);
	}else{ // If the grain is negative...
		Potential = DTOKSsolveOML3( Sec + Therm, 0.01, Potential);
		if( Potential < 0.0 ){ // But if it's positive
			// But! If it's now positive, our assumptions must be wrong!
			// So now we assume it's positive and calculate the potential with a well.
			Potential = DTOKSsolveOML3(0.0, 0.01, Potential)-Td*Kb/(echarge*Te);
		}
	}

	for( double Td(300); Td < 5500; Td += 10){	
		for( double TiTe(0.01); TiTe < 100; TiTe *= 1.1 ){
			double Sec = sec(Td/converteVtoK,'w');
			double gammae = ne*exp(Potential)*sqrt(echarge*Te/(2*PI*Me));
			double Therm = Richardson*pow(Td,2)*exp(-(4.55*echarge)/(Kb*Td))/(echarge*gammae);
			// BEGINING STRUCTURE ONE : NEW CODE
//			if( (Sec + Therm) < 1.0 ){ // DTOKSsolveOML3 only defined for deltatot < 1.0
//				Potential = DTOKSsolveOML3( Sec + Therm, TiTe, Potential); 
//			}else{ // If the grain is in fact positive ...
//				Potential = DTOKSsolveOML3( Sec + Therm, TiTe, Potential);
//				if( Potential < 0.0 ){
//					Potential = DTOKSsolveOML3(0.0, TiTe, Potential)-Kb*Td/(echarge*Te);
//				}
//			}

			// BEGINNING STRUCTURE TWO : OLD CODE
			if( (Sec + Therm) >= 1.0 ){ 
				Potential = DTOKSsolveOML3( 0.0, TiTe, Potential) - Td*Kb/(echarge*Te);
			}else{ // If the grain is negative...
				Potential = DTOKSsolveOML3( Sec + Therm, TiTe, Potential);
				if( Potential < 0.0 ){ // But if it's positive
					// But! If it's now positive, our assumptions must be wrong!
					// So now we assume it's positive and calculate the potential with a well.
					Potential = DTOKSsolveOML3(0.0, TiTe, Potential)-Td*Kb/(echarge*Te);
				}
			}

//			std::cout << "\n" << Td << "\t" << TiTe << "\t" << (Sec+Therm) << "\t" << Potential;
		}
	}
*/
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
//	std::cout << "\n\n*****\n\nUnitTest 4 completed in " << elapsd_secs << "s\n";
}
