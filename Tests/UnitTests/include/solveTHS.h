#include "Constants.h"
#include "Functions.h"
#include "math.h"

// THS Charging Test:
// This test is designed to find the floating potential for small dust grains in magnetised plasmas 
// The calculation follows the work by Drew Thomas, Josh Holgate and Luke Simons and implements a semi-empirical model
// see D. M. Thomas, J. T. Holgate & L. Simons, ArXiv Prepr. (2016). for details
double solveTHS(double Beta, double TemperatureRatio, double MassRatio, double Ionization, double DebyeLength){ 

    double Betai = Beta;
    double Betae = sqrt(TemperatureRatio*MassRatio)*Betai;

    // Model 1: Errorfunctions, Has no initial dip and fits well
//    double BfieldReduceFactor = 0.25;//0.75;
//    double BFieldScaleFactor = 8.0;//8.0;
//    double Coeff = exp(-1.5/(DebyeLength))*sqrt(TemperatureRatio/MassRatio)/(exp(-1.5/(DebyeLength))-BfieldReduceFactor*erf(Betae/BFieldScaleFactor));
//    double a = Coeff*(1.0-BfieldReduceFactor*erf(Betai/BFieldScaleFactor));
//    double Coeff = sqrt(TemperatureRatio/MassRatio)/(1.0-exp(1.5/DebyeLength)*BfieldReduceFactor*erf(Betae/BFieldScaleFactor));
//    double a = Coeff*(1.0-exp(1.5/DebyeLength)*BfieldReduceFactor*erf(Betai/(BFieldScaleFactor*2)));
//    double b = Ionization*Coeff*exp(-(Betai))/TemperatureRatio;

    // Model 2: Exponentials, has a small initial dip that fits worse
//    double BfieldReduceFactor = 0.75;//0.75;
//    double BFieldScaleFactor = 3.0;//8.0;
//    double Coeff =  sqrt(TemperatureRatio/MassRatio)/( exp(-Betae/BFieldScaleFactor)+ exp(-1.0/(DebyeLength))*BfieldReduceFactor*erf(Betae/BFieldScaleFactor));
//    double a = Coeff*(  exp(-1.0/(DebyeLength))*exp(-Betai/BFieldScaleFactor)+ BfieldReduceFactor*erf(Betai/BFieldScaleFactor));
//    double b = Ionization*Coeff*exp(-(Betai))/TemperatureRatio;

    // Model 3: Scale factor function of debye length, Fits upper curve best
//    double BfieldReduceFactor = 0.75;
//    double BFieldScaleFactor = 2.0*(1.0+1.0/(DebyeLength*DebyeLength)); // Similar to exp(1.0/Lambda^2)
//    double Coeff =  sqrt(TemperatureRatio/MassRatio)/(  exp(-Betae/BFieldScaleFactor)+ exp(-1.0/DebyeLength)*BfieldReduceFactor*erf(Betae/BFieldScaleFactor));
//    double a = Coeff*( exp(-1.0/DebyeLength)*exp(-Betai*2.0/BFieldScaleFactor)+BfieldReduceFactor*erf(Betai*2.0/BFieldScaleFactor));
//    double b = Ionization*Coeff*exp(-(Betai)*2.0/BFieldScaleFactor)/TemperatureRatio;

    // Model 4: Scale function of debye length, fits upper curve okay
    double BfieldReduceFactor = 0.85;
    // DebyeLength = DebyeLength*DebyeLength;
    // Betae = Betae*Betae;
    // Betai = Betai*Betai;
    double BFieldScaleFactor = 2.0*exp(1.0/(DebyeLength));//(1.0+1.0/DebyeLength); // Similar to exp(1.0/Lambda)
    double Coeff =  sqrt(TemperatureRatio/MassRatio)/(exp(1.0/(DebyeLength))*exp(-Betae/BFieldScaleFactor)+ BfieldReduceFactor*erf(Betae/BFieldScaleFactor));
    double a = Coeff*( exp(1.0/(DebyeLength))*exp(-Betai*2.0/BFieldScaleFactor)+ exp(1.0/(DebyeLength))*BfieldReduceFactor*erf((Betai)*2.0/BFieldScaleFactor));
    double b = exp(1.0/(DebyeLength))*Ionization*Coeff*exp(-(Betai)*2.0/BFieldScaleFactor)/TemperatureRatio;


    if( std::isinf(exp(a/b)/b) || b == 0.0 ){
//        return log( Coeff*(1.0-BfieldReduceFactor*erf(Betai/BFieldScaleFactor)) );
//        return log( Coeff*BfieldReduceFactor*erf(Betai/BFieldScaleFactor));
//        return log( Coeff*BfieldReduceFactor );
        return log( sqrt(TemperatureRatio/MassRatio) );
    }

    return a/b-LambertW(exp(a/b)/b);
}

