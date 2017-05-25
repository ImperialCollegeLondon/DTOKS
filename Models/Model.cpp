//#define PAUSE
//#define MODEL_DEBUG

#include "Model.h"
/*
	Pdata.NeutralDensity =  3e19;		// m^-3, Neutral density
	Pdata.ElectronDensity = 8e17;	 	// m^-3, Electron density
	Pdata.Potential = 1;			// arb, assumed negative, potential normalised to dust temperature, (-e*phi)/(Kb*Td)
	Pdata.IonTemp = 100*1.16e5;	 	// K, Ion Temperature
	Pdata.ElectronTemp = 100*1.16e5;	// K, Electron Temperature, convert from eV
	Pdata.NeutralTemp = 100*1.16e5; 	// K, Neutral Temperature, convert from eV
*/

const struct PlasmaData ModelDefaults = {
	1e20,		// m^-3, Neutral Density
	1e20,		// m^-3, Electron Density
	1e20,		// m^-3, Electron Density
	0,		// arb, Normalised potential
	116045.25,	// K, Ion Temperature
	116045.25,	// K, Electron Temperature
	116045.25,	// K, Neutral Temperature
	300,		// K, Ambient Temperature
	threevector(),	// m s^-1, Plasma Velocity (Should eventually be normalised to sound speed cs)
	threevector(),	// V m^-1, Electric field at dust location (Normalised later) 
	threevector(),	// T, Magnetic field at dust location (Normalised later)

};

Model::Model():Pdata(ModelDefaults){
}

Model::Model( std::shared_ptr <Matter> const& sample, PlasmaData const& pdata ):Sample(sample),Pdata(pdata){
}

/*void Model::Reset_Data( std::shared_ptr <Matter> const& sample, PlasmaData const& pdata ){
	Sample = sample;
	Pdata = pdata;
}*/
