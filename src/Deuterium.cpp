//#define ELEMENT_DEBUG
//#define ELEMENT_DEEP_DEBUG
#include "Deuterium.h"
#include "Constants.h"

const struct ElementConsts DeuteriumConsts = {
	'D',		// Specifies the element
	18.73, 		// K, Melting temperature at atmospheric pressure https://encyclopedia.airliquide.com/deuterium
	23.31,		// K, Boiling temperature at atmospheric pressure
	443.546,	// Units of kJ/mol, Bond Energy  https://labs.chem.ucsb.edu/zakarian/armen/11---bonddissociationenergy.pdf
	13.6,		// ev, Work Function, from first ionization energy of hydrogen atom
	1.0, 		// HeatTransfer coefficient (estimate!)
	0.004032,	// AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
	49.261, 	// kJ/kg, Latent Fusion Energy From Wikipedia
	322.215, 	// kJ/kg, Latent Vapour Energy From Wikipedia
	3.5*10-3,	// N/m, Surface Tension B. Keene, 1993, Review of data for the surface tension of pure metal
	171.0,		// (kg/m^3) from https://www.bnl.gov/magnets/staff/gupta/cryogenic-data-handbook/Section4.pdf
	0.0001382,	// kW/mK, invalid for liquid hydrogen
};


Deuterium::Deuterium():Matter(DeuteriumConsts){
	set_defaults();
	update();
}

Deuterium::Deuterium(double radius):Matter(radius,DeuteriumConsts){
	set_defaults();
	update();
}

Deuterium::Deuterium(double radius, double tempin):Matter(radius,tempin,DeuteriumConsts){
	E_Debug("\n\nIn Deuterium::Deuterium(double radius, double tempin)");
	set_defaults();

	update_state(0.0);		// Temperature dependant
	update_models('c','c','c','y','n');
	update();

	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Deuterium::Deuterium(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,DeuteriumConsts){
	E_Debug("\n\nIn Deuterium::Deuterium(double radius, double tempin, std::array<char,CM> &constmodels)");
	set_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	update();

	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Deuterium::Deuterium(double radius, double tempin, std::array<char,CM> &constmodels, const threevector& position, const threevector& velocity):Matter(radius,tempin,DeuteriumConsts){
	E_Debug("\n\nIn Deuterium::Deuterium(double radius, double tempin, std::array<char,CM> &constmodels, const threevector & position, const threevector& velocity)");
	set_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	update_motion(position,velocity,0.0);
	update();

	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

void Deuterium::set_defaults(){
	E_Debug("\n\nIn Deuterium::set_defaults()");
	// http://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.html
	St.HeatCapacity = 2.950;	 	// kJ/(kg-K), https://en.wikipedia.org/wiki/Deuterium
	St.Emissivity = 0.01; 			// Arb, A guess for Deuterium
	St.SuperBoilingTemp = Ec.BoilingTemp; 	// K, At any pressure
	St.Density = 162.4; 			// (kg/m^3) from Wikipedia
	update_models('c','c','c','y','n');
}

/// ************************************************************************************************* \\\

void Deuterium::update_heatcapacity(){ // Calculates the heat capacity in units of J mol^-1 K^-2
	E_Debug("\n\nIn Deuterium::update_heatcapacity()");
	static bool runOnce = true;
	WarnOnce(runOnce,"Deuterium HeatCapacity assumed constant! Variable heat capacity not possible");
	St.HeatCapacity = St.HeatCapacity;
}

void Deuterium::update_radius(){
	E_Debug("\n\nIn Deuterium::update_radius():");
	St.LinearExpansion=1.0;
	static bool runOnce = true;
	WarnOnce(runOnce,"Deuterium LinearExpansion == 1.0 assumed! Temperature dependant radius not possible");
	St.Radius=St.UnheatedRadius*St.LinearExpansion;
	assert(St.Radius>0); // Assert radius is positive	
}

void Deuterium::update_vapourpressure(){
	St.VapourPressure = probe_vapourpressure(St.Temperature);
}

double Deuterium::probe_vapourpressure(double Temperature)const{
	double VapourPressure(0.0);
	if( !St.Liquid && !St.Gas ){
		// Page 466, equation 7.15
		// H.W. Woolley, R.B. Scott, and F.G. Brickwedde, J. Res. Natl. Bur. Stand. (1934). 41, 379 (1948).
		VapourPressure = pow(10,5.1625-67.9119/Temperature+0.03102*Temperature);
	}else if( St.Liquid ){
		// Page 466, equation 7.14
		// H.W. Woolley, R.B. Scott, and F.G. Brickwedde, J. Res. Natl. Bur. Stand. (1934). 41, 379 (1948).
		VapourPressure = pow(10,4.7367-58.54440/Temperature+0.02670*Temperature);
	}else{
		VapourPressure = 0;
		std::cout << "\nWarning, Sample is assumed gas! St.VapourPressure = 0";
	}
	VapourPressure = 133.322*VapourPressure; // Convert to Pascals
	return VapourPressure;
}