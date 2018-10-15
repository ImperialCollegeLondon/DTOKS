#include "PlasmaData.h"
#include "Functions.h"
#include "Constants.h"
#include <iostream>
#include <math.h>

// Investigate the influence of fraction of backscattered ions and and fraction of backscattered energy
// This test prints the values of the fraction of backscattered energy and the fraction of back scattered particles
// as calculated by the backscatter function from DTOKS. This can be readily compared to the results published in
// the DTOKS papers
void BackscatterTest(){
	clock_t begin = clock();

	double RE(0), RN(0);

	PlasmaData Pdata;
	Pdata.NeutralDensity =  2e17;//10e17;	// m^-3, Neutral density
	Pdata.ElectronDensity = 2e19;//10e18; 	// m^-3, Electron density
	double Potential = 1;			// arb, assumed negative, potential normalised to dust temperature, (-e*phi)/(Kb*Td)
	double NumOfev = 1;
	Pdata.IonTemp = NumOfev*1.16e4;	 	// K, Ion Temperature
	Pdata.ElectronTemp = NumOfev*1.16e4; 	// K, Electron Temperature, convert from eV
	Pdata.NeutralTemp = NumOfev*1.16e4; 	// K, Neutral Temperature, convert from eV

	for(int i(0); i < 5; i ++){ // Loop over Electron Temperature
		for(int l(0); l < 5; l ++){ // Loop over Ion Temperatures
			for(int k(0); k < 10; k ++){ // Loop over Dust potentials
				double t1 = backscatter(Pdata.ElectronTemp*pow(10,i-2),Pdata.IonTemp*pow(10,l-2),Mp,
							Potential*((k-5)*0.2),'w',RE,RN);
				std::cout << "\n" << t1 << "\t" << Pdata.ElectronTemp*pow(10,i-2) 
						<< "\t" << Pdata.IonTemp*pow(10,l-2) << "\t" << Potential*((k-5)*0.2) 
						<< "\t" << RE << "\t" << RN;
			}
		}
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
//	std::cout << "\n\n*****\n\nUnitTest 1 completed in " << elapsd_secs << "s\n";
}
