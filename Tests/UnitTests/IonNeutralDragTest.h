#include <iostream>
#include "Constants.h"
#include "MathHeader.h" 

double NeutralDragForce(threevector Vrel, double NeutralDensity, double NeutralTemp, double radius){
	return Vrel.mag3()*Mp*sqrt(4*PI)*NeutralDensity*sqrt(Kb*NeutralTemp/(2*PI*Mp)) *PI*pow(radius,2);
}

double IonDragForce(threevector Vrel, double IonDensity, double IonTemp, double ElectronTemp, double radius, double potential){
	double ConvertKelvsToeV(8.621738e-5);
	threevector Mt = Vrel*sqrt(Mp/(Kb*IonTemp)); 
	double lambda = sqrt(epsilon0/(IonDensity*echarge*exp(-Mt.mag3()*Mt.mag3()/2)
		*(1.0/(IonTemp*ConvertKelvsToeV))+1.0/(ElectronTemp*ConvertKelvsToeV)));
	double beta = ElectronTemp*ConvertKelvsToeV*radius
		*fabs(potential)/(lambda*IonTemp*ConvertKelvsToeV);
	double Lambda = -exp(beta/2.0)*Exponential_Integral_Ei(-beta/2.0); 
	threevector FidS = Mt*(sqrt(32*PI)/3.0*epsilon0*pow(IonTemp*ConvertKelvsToeV,2)*Lambda*pow(beta,2));
	threevector FidC = Vrel*4.0*PI*pow(radius,2)*IonDensity*Mp
		*sqrt(Kb*ElectronTemp/(2.0*PI*Me))*exp(-potential);
	return ((FidS+FidC).mag3());
}

// Test orbital motion limited theory
void IonNeutralDragTest(){
	clock_t begin = clock();

	double ConvertKelvsToeV(8.621738e-5);

	threevector PlasmaVel(0.0,0.0,0.0);
	threevector DustVel(0.0,0.0,10.0);
	threevector RelativeVelocity = PlasmaVel-DustVel;
	double ElectronTemp = 10/ConvertKelvsToeV;		// 10ev converted to Kelvin
	double Radius = 1e-6;					// metres, radius of dust
	double Potential = 1;
	int i(0), j(0);
	for( double Density(0.01e18); Density < 100e18; Density *= 1.5 ){
		i = i + 1;
		for( double Temp(0.1/ConvertKelvsToeV); Temp < 1000/ConvertKelvsToeV; Temp *= 1.5 ){
			if( i == 1 ) j = j + 1;

			double Iondrag = IonDragForce(RelativeVelocity,Density,Temp,ElectronTemp,Radius,Potential);
			double Neutraldrag = NeutralDragForce(RelativeVelocity,Density,Temp,Radius);
			std::cout << RelativeVelocity.mag3() << "\t" << Density << "\t" << Temp*ConvertKelvsToeV << "\t"  << ElectronTemp << "\t" 
				<< Radius << "\t" << Potential << "\t" << Iondrag << "\t" << Neutraldrag << "\t" << Iondrag/Neutraldrag << "\n";
		}
	}
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIonNeutralDrag UnitTest completed in " << elapsd_secs << "s\n" << i << " : " << j << "\n";
}
