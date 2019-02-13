#include "Constants.h"
#include "Functions.h"
#include <math.h>

// Solve the Modified orbital motion limited potential for large emitting dust 
// grains.
// See the paper by Minas and Nikoleta, equation (1) and (2)
// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
double solveMOMLWEM(double TemperatureRatio, double MassRatio, 
    double Ionization, double HeatCapacityRatio, double DeltaTot){

	double Delta_Phi_em = 0.5*log((2.0*PI/MassRatio)*
        (1+HeatCapacityRatio*TemperatureRatio)/pow(1.0-DeltaTot,2.0));
	// Uncomment following line to compare MOMLWEM results with figure (1) of 
    // paper:
	// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
	std::cout << DeltaTot << "\t" << Delta_Phi_em << "\n";
	return TemperatureRatio/Ionization+Delta_Phi_em/Ionization
		-LambertW((1.0-DeltaTot)*sqrt(MassRatio*TemperatureRatio)
		*exp(TemperatureRatio/Ionization+Delta_Phi_em/Ionization)/Ionization);
}

