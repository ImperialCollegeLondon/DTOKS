#include "Functions.h"
#include "Constants.h"

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double solveSOML(double TemperatureRatio, double PlasmaFlowSpeed, double MassRatio){
	PlasmaFlowSpeed = PlasmaFlowSpeed;
	double s1 = sqrt(PI)*(1+2*pow(PlasmaFlowSpeed,2))*erf(PlasmaFlowSpeed)/(4*PlasmaFlowSpeed)+exp(-pow(PlasmaFlowSpeed,2))/2;
	double s2 = sqrt(PI)*erf(PlasmaFlowSpeed)/(2*PlasmaFlowSpeed);

	return TemperatureRatio*s1/s2-LambertW(sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio*s1/s2)/s2);
}
/* 17/11/18, got half-way through implementing this, realised it might become too complicated.
	C_Debug("\tIn ChargingModel::solveSOML(double guess)\n\n");
Will continue to test inside of DTOKS-U code
double solveSOML_Plus(double IonDensity, double IonTemperature, double ElectronTemperature, 
			double PlasmaFlowSpeed, double MassRatio, double Z, double guess){

	double IonThermalVelocity = sqrt((Kb*IonTemperature)/(MassRatio*Me));
	double Tau = IonTemperature/ElectronTemperature;
	double uz = PlasmaFlowSpeed/IonThermalVelocity;

	
	double IonFluxDiff(0.0), ElectronFluxDiff(0.0), x1(0.0), dT(0.0), dTdiff(0.0);
	int i(0);
	double guessinit = guess;
	do{
		if( i > 0 ) // If first time in loop, don't update guess. Avoids numerical instabilities.
			guess = x1;
		if( i == 1000 ){ // Too many loops! Maybe we need to flip sign?

			std::cerr << "\nError in ChargingModel::solveSOML()! Newton Rhapson Method didn't converge to accuracy: " << Accuracy;
			return 0.0;
		}

		if( guessinit >= 0.0 ){ // OML Electron flux, SOML Ion flux
			if( uz == 0.0 ){
				IonFluxDiff = IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*Z/Tau;
			}else{
				IonFluxDiff = (IonDensity*IonThermalVelocity/(2.0*sqrt(2.0)*uz))*
					Z*erf(uz)/Tau;
			}
			ElectronFluxDiff = -OMLElectronFlux(fabs(guess));
			dT = Richardson*pow(Sample->get_temperature(),2)*exp(-echarge*Sample->get_workfunction()/(Kb*Sample->get_temperature()))/echarge;
			dTdiff = 0.0;
		}else{ // OML Electron flux, SOML Ion flux
			double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
			double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
			if( uz == 0.0 ){
				IonFluxDiff = (Pdata->ElectronTemp/Pdata->IonTemp)*SOMLIonFlux(guess);
			}else{
				IonFluxDiff = Pdata->IonDensity*
					IonThermalVelocity*(1.0/(4.0*sqrt(2.0)*uz))*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm));
			}
			
			ElectronFluxDiff = -Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
			dT = Richardson*pow(Sample->get_temperature(),2)*(1.0-guess
					*(Pdata->ElectronTemp/Sample->get_temperature()))
					*exp((-echarge*Sample->get_workfunction()+guess*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge;
			dTdiff = (Pdata->ElectronTemp/Sample->get_temperature())*(dT-
					Richardson*pow(Sample->get_temperature(),2)
					*exp((-echarge*Sample->get_workfunction()+guess*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge);
		}

		x1 = guess - ( (( 1.0-DeltaSec())*OMLElectronFlux(guess) - dT - SOMLIonFlux(guess))
			/((1.0-DeltaSec())*ElectronFluxDiff - dTdiff - IonFluxDiff ) );
		i ++;
	}while( fabs(guess-x1) > Accuracy );

	return guess;
}
*/
