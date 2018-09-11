#include "Constants.h"
#include "Functions.h"
#include "math.h"

// TH Charging Test:
// This test is designed to find the floating potential for small dust grains in magnetised plasmas 
// The calculation follows the work by Drew Thomas and Josh Holgate and implements a semi-empirical model
// see D. M. Thomas and J. T. Holgate, ArXiv Prepr. (2016). for details
double solveTH(double Beta, double TemperatureRatio, double MassRatio, double Ionization, double DebyeLength){ 

	double Betai = Beta;
	double Betae = sqrt(TemperatureRatio*MassRatio)*Betai;

	double BfieldReduceFactor = 0.25;//0.75;
	double BFieldScaleFactor = 8.0;//8.0;
	double Coeff = exp(-1.5/(DebyeLength))*sqrt(TemperatureRatio/MassRatio)/(exp(-1.5/(DebyeLength))-BfieldReduceFactor*erf(Betae/BFieldScaleFactor));
	double a = Coeff*(1.0-BfieldReduceFactor*erf(Betai/BFieldScaleFactor));
	double b = Ionization*Coeff*exp(-(Betai))/TemperatureRatio;

//	double BfieldReduceFactor = 0.75;//0.75;
//	double BFieldScaleFactor = 3.0;//8.0;
//
//	
//	double Coeff =  sqrt(TemperatureRatio/MassRatio)/( exp(-Betae/BFieldScaleFactor)+ exp(-1.0/(DebyeLength))*BfieldReduceFactor*erf(Betae/BFieldScaleFactor));
//	double a = Coeff*(  exp(-1.0/(DebyeLength))*exp(-Betai/BFieldScaleFactor)+ BfieldReduceFactor*erf(Betai/BFieldScaleFactor));
//	double b = Ionization*Coeff*exp(-(Betai))/TemperatureRatio;





	if( std::isinf(exp(a/b)/b) || b == 0.0 ){
		return log( Coeff*(1.0-BfieldReduceFactor*erf(Betai/BFieldScaleFactor)) );
//		return log( Coeff*BfieldReduceFactor*erf(Betai/BFieldScaleFactor));
	}

	return a/b-LambertW(exp(a/b)/b);
}

