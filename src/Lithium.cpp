//#define ELEMENT_DEBUG
//#define ELEMENT_DEEP_DEBUG

#include "Lithium.h"
#include "Constants.h"

const struct ElementConsts LithiumConsts = {
	'L',			// Specifies the element
	453.65, 			// K, Melting temperature at atmospheric pressure
	1603.0,			// K, Boiling temperature at atmospheric pressure
	520.0,			// Units of kJ/mol, Bond Energy, Estimated as being equal to Latent Vapour Energy
	2.9,			// ev, Work Function, Taken from DTOKS and matched with wikipedia
	1.0,	 		// kJ/(m^2 K), HeatTransfer coefficient, 
					// http://www.engineeringtoolbox.com/overall-heat-transfer-coefficients-d_284.html
	0.006941,		// AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
	432.28, 		// kJ/kg, Latent Fusion Energy From Wikipedia
	19593.71, 		// kJ/kg, Latent Vapour Energy From Wikipedia
	0.35,			// N/m, Surface Tension B. Keene, 1993, Review of data for the surface tension of pure metals
	534,			// (kg/m^3) from Wikipedia, Room temperature density
	0.0848			// kW/m K at 20 degrees celsius
};

Lithium::Lithium():Matter(LithiumConsts){
	E_Debug("\n\nIn Lithium::Lithium():Matter(&LithiumConsts)\n\n");
	lithium_defaults();
	update();
}

Lithium::Lithium(double radius):Matter(radius,LithiumConsts){
	E_Debug("\n\nIn Lithium::Lithium(double radius):Matter(radius,&LithiumConsts)\n\n");
	lithium_defaults();
	update();
}

Lithium::Lithium(double radius, double tempin):Matter(radius,tempin,LithiumConsts){
	E_Debug("\n\nIn Lithium::Lithium(double radius, double tempin):Matter(radius,tempin,&LithiumConsts)\n\n");
	lithium_defaults();
	update_state(0.0);		// Temperature dependant
	update_models('c','c','c','y','n');
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Lithium::Lithium(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,LithiumConsts){
	E_Debug("\n\nIn Lithium::Lithium(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,&LithiumConsts)\n\n\t");
	lithium_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Lithium::Lithium(double radius, double tempin, std::array<char,CM> &constmodels, const threevector& position, const threevector& velocity):Matter(radius,tempin,LithiumConsts){
	E_Debug("\n\nIn Lithium::Lithium(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,&LithiumConsts)\n\n\t");
	lithium_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	update_motion(position,velocity,0.0);
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

void Lithium::lithium_defaults(){
	E_Debug("\tIn Lithium::lithium_defaults()\n\n");

	// https://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.htm
	St.HeatCapacity = 3.55878; 		// kJ/(kg-K)         
	St.Emissivity = 0.1;			// http://www.fusion.ucla.edu/APEX/meeting15/Apex4_01-tanaka.pdf
	St.SuperBoilingTemp = Ec.BoilingTemp; 	// K, At any pressure
	St.Density = 534;			// (kg/m^3) from Wikipedia
	update_models('c','c','c','y','n');
	

//	St.ThermConduct = 0.163;		// kW/m K at 20 degrees celsius
}

/// ************************************************************************************************* \\\

void Lithium::update_heatcapacity(){ // Calculates the heat capacity in units of J mol^-1 K^-2
	E_Debug("\tIn Lithium::update_heatcapacity()\n\n");
	/*
	if( St.Temperature < 300 ){ 

		static bool runOnce = true;
		WarnOnce(runOnce,
			"In Lithium::update_heatcapacity():\nExtending heat capacity model outside temperature range! T < 300K");

		St.HeatCapacity = 24.943 - 7.72e4*pow(St.Temperature,-2) + 2.33e-3*St.Temperature + 1.18e-13*pow(St.Temperature,4);
	}else if( St.Temperature > 300 && St.Temperature < Ec.MeltingTemp ){ 
		// http://nvlpubs.nist.gov/nistpubs/jres/75a/jresv75an4p283_a1b.pdf
		St.HeatCapacity = 24.943 - 7.72e4*pow(St.Temperature,-2) + 2.33e-3*St.Temperature + 1.18e-13*pow(St.Temperature,4);
	}else{
		// http://webbook.nist.gov/cgi/inchi?ID=C7440337&Mask=2

		double t = St.Temperature/1000;
                St.HeatCapacity = 35.56404 -1.551741e-7*t + 2.915253e-8*pow(t,2) -1.891725e-9*pow(t,3)-4.107702e-7*pow(t,-2);
	}

//	E1_Debug("\n\nTemperature is : " << St.Temperature << "\nSt.Gas = " << St.Gas << "\nSt.Liquid = " << St.Liquid 
//				<< "\nCv of Solid: " << St.HeatCapacity/Ec.AtomicMass << "[kJ/(kg K)]"; );
	St.HeatCapacity = (St.HeatCapacity /(1000 * Ec.AtomicMass)); // Conversion kJ/(mol K) to kJ/( kg K ), AtomicMass [kg mol^-1]

	*/
	St.HeatCapacity = St.HeatCapacity;
}

void Lithium::update_radius(){
	E_Debug("\tIn Lithium::update_radius()\n\n");
	/*
	if( St.Temperature > 173 && St.Temperature <= 1500 ){
		static bool runOnce = true;
		WarnOnce(runOnce,
				"In Lithium::update_radius():\nExtending model outside temperature range!(from 738K to 1500K)");
		St.LinearExpansion = 1+(4.28*St.Temperature)*1e-6;
	}else if( St.Temperature > 1500 && St.Temperature < Ec.MeltingTemp ){
		St.LinearExpansion = 1+3.9003e-4+1.3896*1e-3-8.2797*1e-7*St.Temperature + 4.0557*1e-9*pow(St.Temperature,2)
					-1.2164*1e-12*pow(St.Temperature,3)+1.7034*1e-16*pow(St.Temperature,4);
		// 1.02648269488
	}else if( St.Temperature == Ec.MeltingTemp ){ // Model while it's melting
		// http://nvlpubs.nist.gov/nistpubs/jres/75a/jresv75an4p283_a1b.pdf
		St.LinearExpansion = 1.02105 + 0.03152*St.FusionEnergy/(Ec.LatentFusion*St.Mass);
		// 1.02648269488 - 1.05257238281 = 0.02608968793
	}else if( St.Temperature > Ec.MeltingTemp && St.Temperature <= St.SuperBoilingTemp){
		St.LinearExpansion = pow(1.18+6.20*1e-5*(St.Temperature-3680)+3.23*1e-8*pow((St.Temperature-3680),2),(1./3.));
	}*/
	St.LinearExpansion = 1.0;
	St.Radius=St.UnheatedRadius*St.LinearExpansion;
	E1_Debug("\nTemperature = " << St.Temperature << "\n\nSt.LinearExpansion = " << St.LinearExpansion
			<< "\nSt.Radius = " << St.Radius);
	assert(St.Radius>0); // Assert radius is positive	
}

void Lithium::update_vapourpressure(){
	E_Debug("\tIn Lithium::update_vapourpressure()\n\n");
//	double AmbientPressure = 0;
	St.VapourPressure = probe_vapourpressure(St.Temperature);
	// Answer is in pascals
}

double Lithium::probe_vapourpressure(double Temperature)const{
	double VapourPressure(0.0);
	if( Temperature < 298 ){
		static bool runOnce = true;
		WarnOnce(runOnce,
				"In Lithium::update_vapourpressure():\nExtending model outside temperature range!(from 298K< to 0K<)");
		VapourPressure = 101325*pow(10,5.667 - 8310/Temperature); // http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
	}else if( Temperature >= Ec.MeltingTemp && Temperature < 800 ){
		VapourPressure = 101325*pow(10,5.055 - 8023/Temperature); // http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
	}else if( Temperature >= 800 ){
		// https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19680018893.pdf
		VapourPressure = 101325*pow(10,10.015-8064.5/Temperature);
	}else{
		std::cerr << "\nError! Negative Temperature in Lithium::probe_vapourpressure(double Temperature)";
	}
	
	return VapourPressure;
}