#include "Constants.h"
#include <math.h>
#include <iostream>

double Tungsten_VapourPressure(double Temperature){
//	double AmbientPressure = 0;
	return 101325*pow(10,2.945 - 44094/Temperature + 1.3677*log10(Temperature)); // http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
	
	// Answer is in pascals
}

double EvaporationFlux(double DustTemperature, double SA, double AtomicMass){
//	H_Debug("\n\nIn HeatingModel::EvaporationFlux():");
	double AmbientPressure = 0;
	double StickCoeff = 1.0;
	return (StickCoeff*SA*AvNo*(Tungsten_VapourPressure(DustTemperature)-AmbientPressure))/
			sqrt(2*PI*AtomicMass*R*DustTemperature);
}

// This test was performed to investigate the relative influence of variation in mass and heat capacity with time
// on the change in temperature with time. This was performed below for Tungsten
// There is a major issue with calculating the change in internal energy, as the values have been given for a 'power' 
// but these should be multiplied by a time step. Essential, the equation dT = mCv dE is not well defined here when taking
// A partial derrivative with respect to time. As we are left with terms like dm/dt Cv dE where dE is not well defined.
// This analysis found that the change in heat capacity with time was of the order of 1e-9 relative to one.

int TemperatureEvolutionEquation(char Element){
	clock_t begin = clock();
	// ********************************************************** //

	int ReturnVal = 0;

	double mass = 1e-6;		//kg
	double density = 19600;		// kg/m^3
	double SA = 4*PI*pow((3*mass)/(4*PI*density),2./3.);
	double AtomicMass = 0.18384; // kg
	double HeatCapacity = 35/(1000 * AtomicMass); // Convert from kJ/mol K to J/ kJ K
	double Power = 1;		// kW
	double Temp = 5500; 		// K
	double Bondenergy = 774;	// kJ/mol

	double EvapFlux = EvaporationFlux(Temp,SA,AtomicMass);

	std::cout << "\nTerm1 = " << mass*HeatCapacity*(Power-EvapFlux*((3*Kb*Temp)/(2*1000)+Bondenergy/AvNo)) << "K/s";	
//	std::cout << "\nbit1 = " << mass*HeatCapacity*Power << "K/s";	
//	std::cout << "\nbit2 = " << -mass*HeatCapacity*EvapFlux*((3*Kb*Temp)/(2*1000)+Bondenergy/AvNo) << "K/s";
	
	std::cout << "\nTerm2 = " << EvapFlux*(AtomicMass/AvNo)*HeatCapacity*(Power-EvapFlux*((3*Kb*Temp)/(2*1000)+Bondenergy/AvNo)) << "K/s^2";
//	std::cout << "\nEvapFlux = " << EvapFlux;
//	std::cout << "\n(AtomicMass/AvNo) = " << (AtomicMass/AvNo);
	
	
	double HeatCapacityDiff = -1.551741e-7 + 2*2.915253e-8*Temp - 3*1.891725e-9*Temp+2*4.107702e-7*Temp;
//	std::cout << "\nHeatCapacityDiff = " << HeatCapacityDiff;
//	std::cout << "\nmass = " << mass;
//	std::cout << "\nPower-EvapFlux*((3*Kb*Temp)/(2*1000)+Bondenergy/AvNo) = " << Power-EvapFlux*((3*Kb*Temp)/(2*1000)+Bondenergy/AvNo);
	double term3 = 1-mass*HeatCapacityDiff*(Power-EvapFlux*((3*Kb*Temp)/(2*1000)+Bondenergy/AvNo));
	std::cout << "\nTerm3 = " << term3 << "K/s^2";
	std::cout << "\nbit2 = " << -mass*HeatCapacityDiff*(Power-EvapFlux*((3*Kb*Temp)/(2*1000)+Bondenergy/AvNo));
	std::cout << "\nEvapFlux*Ma = " << EvapFlux*AtomicMass/AvNo;

	// ********************************************************** //
	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nErrorEstimate 1 completed in " << elapsd_secs << "s\n";

	return ReturnVal;
}
