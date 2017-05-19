#ifndef __PLASMADATA_H_INCLUDED__   // if PlasmaData.h hasn't been included yet...
#define __PLASMADATA_H_INCLUDED__

// Structure containing all the information defining the plasma conditions
struct PlasmaData{
	// Parameters used to calculate particle fluxes
	double NeutralDensity;	// m^-3, Neutral Density
	double ElectronDensity;	// m^-3, Electron Density
	double Potential;	// arb, Normalised potential
	double IonTemp;		// K, Ion Temperature
	double ElectronTemp;	// K, Electron Temperature
	double NeutralTemp;	// K, Neutral Temperature
	double AmbientTemp;	// K, Ambient Temperature
};

#endif
