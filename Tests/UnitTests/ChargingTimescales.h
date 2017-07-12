#include "Constants.h"
#include "DTOKSsolveOML.h"
#include <iostream>

// Plot results using something
void ChargingTimescales(){
	clock_t begin = clock();

	double ePlasmaFreq, iPlasmaFreq, DebyeLength, timestep;
	for( double ne(1e18); ne < 1e20; ne *= 1.1){ // Loop over Electron density
		ePlasmaFreq = sqrt((ne*pow(echarge,2))/(epsilon0*Me));
		iPlasmaFreq = sqrt((ne*pow(echarge,2))/(epsilon0*Mp));
		for( double Te(1); Te < 100; Te *= 1.1 ){ // Loop over Electron Temperature
			DebyeLength=sqrt((epsilon0*echarge*Te)/(ne*pow(echarge,2)));
//			To just compare the plasma freq, debyelength as a function of temperature and electron density
//			std::cout << "\n" << ne << "\t" << Te << "\t" << ePlasmaFreq << "\t" << iPlasmaFreq << "\t" << DebyeLength;
			for( double Ti(1); Ti < 100; Ti *= 1.1 ){ // Loop over Ion Temperature
                                for( double a(1e-7); a < 1e-6; a *= 1.1 ){ // Loop over radius
					timestep = sqrt(2*PI)*((DebyeLength)/a)*(1/
					(ePlasmaFreq*(1+Te/Ti+fabs(DTOKSsolveOML(Ti,Te,-0.5)))));
					// Proof that ePlasmaFreq is never quicker than timestep (i.e, we expect no output)...
					if( 1/ePlasmaFreq > timestep ){
						// Proof that timestep is always smaller than ePlasmaFreq
						std::cout << "\n" << 1/ePlasmaFreq << "\t" << Te/Ti << "\t" 
							<< DTOKSsolveOML(Ti,Te,-0.5) << "\t" <<  DebyeLength/a << "\t" << timestep;
					}
				}
			}
		}
	}

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
//	std::cout << "\n\n*****\n\nUnitTest 4 completed in " << elapsd_secs << "s\n";
}
