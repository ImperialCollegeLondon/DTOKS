#include "Constants.h"
#include "MathHeader.h"
#include <iostream>

// This test outputs the potential as calculated by part of the DTOKS solution to the OML equation, considering a potential well.
void FortovIonDragTest(){
	clock_t begin = clock();

	// TEST TO CALCULATE THE DTOKS FLOATING POTENTIAL FOR CONSTANT ELECTRON DENSITY AND ELECTRON TEMPERATURE

	double Radius = 1e-6; 		// Dust radius in m
	double IonDensity = 1e18; 	// Ion Density in m^-3
	double PlasmaVel = 1e4;		// PlasmaVelocity in m s^-^1

	double ConvertKelvsToeV(8.621738e-5);

	for( double u(0.01); u < 100; u *= 1.1){	// Normalised ion flow velocity
		for( double Te(10000); Te < 10001; Te *= 1.1 ){	// Electron temperature, K
			for( double Ti(5000); Ti < 50001; Ti *= 1.1 ){	// Ion temperature, K
				for( double potential(0.1); potential < 5.1; potential += 0.1 ){	// Normalised dust potential

					double Fid(0);
					if(u<2.0){ // Relative speed less than twice mach number, use Fortov et al theory with screening length 'Lambda'.
			
						double lambda = sqrt(epsilon0/(IonDensity*echarge*exp(-u*u/2)
							*(1.0/(Ti*ConvertKelvsToeV))+1.0/(Te*ConvertKelvsToeV)));
						double beta = Te*ConvertKelvsToeV*Radius
							*fabs(potential)/(lambda*Ti*ConvertKelvsToeV);
						F_Debug("\nlambda = " << lambda << "\nbeta = " << beta << "\nPot = " << potential);
	
						double Lambda = -exp(beta/2.0)*Exponential_Integral_Ei(-beta/2.0); 
						double FidS = u*(sqrt(32*PI)/3.0*epsilon0*pow(Ti*ConvertKelvsToeV,2)
								*Lambda*pow(beta,2));
			
						double FidC = u*sqrt(Kb*Ti/Mp)*4.0*PI*pow(Radius,2)*IonDensity*Mp
							*sqrt(Kb*Te/(2.0*PI*Me))*exp(-potential); 
						//for John's ion drag... I assume here and in other places in the 
						//calculation that the given potential is normalised to kTe/e
				    
						//Do The same but with the ion current instead of the electron current
				    
						//.....
				    
						Fid=FidS+FidC;
//						Fid=0.1*(FidS+FidC); // Fudge
						F_Debug("\nFidS = " << FidS << "\nFidC = " << FidC);
					}else{	// Relative speed greater than twice the mach number, use just plain collection area
						Fid = u*u*PI*Ti*ConvertKelvsToeV*pow(Radius,2)*IonDensity*echarge;
					}

					std::cout << "\n" << u << "\t" << Te << "\t" << Ti << "\t" << Ti/Te << "\t"
						<< potential << "\t" << Fid;
				}
			}
		}
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nFortovIonDrag UnitTest completed in " << elapsd_secs << "s\n";
}
