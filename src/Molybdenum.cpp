//#define ELEMENT_DEBUG
//#define ELEMENT_DEEP_DEBUG

#include "Molybdenum.h"
#include "Constants.h"

const struct ElementConsts MolybdenumConsts = {
	'M',			// Specifies the element
	2896.0, 			// K, Melting temperature at atmospheric pressure
	4912.0,			// K, Boiling temperature at atmospheric pressure
	681.0,			// Units of kJ/mol, Bond Energy, Estimated as being equal to Latent Vapour Energy
	4.2,			// ev, Work Function, Taken from DTOKS and matched with wikipedia
	1.0,	 		// kJ/(m^2 K), HeatTransfer coefficient, 
					// http://www.engineeringtoolbox.com/overall-heat-transfer-coefficients-d_284.html
	0.09595,		// AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
	390.62, 		// kJ/kg, Latent Fusion Energy From Wikipedia
	6232.41, 		// kJ/kg, Latent Vapour Energy From Wikipedia
	2.24,			// N/m, Surface Tension B. Keene, 1993, Review of data for the surface tension of pure metals
	10280,			// (kg/m^3) from Wikipedia, Room temperature density
	0.138			// kW/m K at 20 degrees celsius
};

Molybdenum::Molybdenum():Matter(MolybdenumConsts){
	E_Debug("\n\nIn Molybdenum::Molybdenum():Matter(&MolybdenumConsts)\n\n");
	molybdenum_defaults();
	update();
}

Molybdenum::Molybdenum(double radius):Matter(radius,MolybdenumConsts){
	E_Debug("\n\nIn Molybdenum::Molybdenum(double radius):Matter(radius,&MolybdenumConsts)\n\n");
	molybdenum_defaults();
	update();
}

Molybdenum::Molybdenum(double radius, double tempin):Matter(radius,tempin,MolybdenumConsts){
	E_Debug("\n\nIn Molybdenum::Molybdenum(double radius, double tempin):Matter(radius,tempin,&MolybdenumConsts)\n\n");
	molybdenum_defaults();
	update_state(0.0);		// Temperature dependant
	update_models('c','c','c','y','n');
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Molybdenum::Molybdenum(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,MolybdenumConsts){
	E_Debug("\n\nIn Molybdenum::Molybdenum(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,&MolybdenumConsts)\n\n\t");
	molybdenum_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Molybdenum::Molybdenum(double radius, double tempin, std::array<char,CM> &constmodels, const threevector& position, const threevector& velocity):Matter(radius,tempin,MolybdenumConsts){
	E_Debug("\n\nIn Molybdenum::Molybdenum(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,&MolybdenumConsts)\n\n\t");
	molybdenum_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	update_motion(position,velocity,0.0);
	update();
	E1_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

void Molybdenum::molybdenum_defaults(){
	E_Debug("\tIn Molybdenum::molybdenum_defaults()\n\n");

	// http://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.html
	St.HeatCapacity = 0.27716616; 		// kJ/(kg-K)         
	St.Emissivity = 0.1; 			// Arb, http://www.engineeringtoolbox.com/emissivity-coefficients-d_447.html
	St.SuperBoilingTemp = Ec.BoilingTemp; 	// K, At any pressure
	St.Density = 10280;			// (kg/m^3) from Wikipedia
	update_models('c','c','c','y','n');
	

//	St.ThermConduct = 0.163;		// kW/m K at 20 degrees celsius
}

/// ************************************************************************************************* \\\

void Molybdenum::update_heatcapacity(){ // Calculates the heat capacity in units of J mol^-1 K^-2
	E_Debug("\tIn Molybdenum::update_heatcapacity()\n\n");
/*
	if( St.Temperature < 300 ){ 

		static bool runOnce = true;
		WarnOnce(runOnce,
			"In Molybdenum::update_heatcapacity():\nExtending heat capacity model outside temperature range! T < 300K");

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

void Molybdenum::update_radius(){
	E_Debug("\tIn Molybdenum::update_radius()\n\n");
	/*
	if( St.Temperature > 173 && St.Temperature <= 1500 ){
		static bool runOnce = true;
		WarnOnce(runOnce,
				"In Molybdenum::update_radius():\nExtending model outside temperature range!(from 738K to 1500K)");
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
	}
	*/
	St.LinearExpansion = 1.0;
	St.Radius=St.UnheatedRadius*St.LinearExpansion;

	E1_Debug("\nTemperature = " << St.Temperature << "\n\nSt.LinearExpansion = " << St.LinearExpansion
			<< "\nSt.Radius = " << St.Radius);
	assert(St.Radius>0); // Assert radius is positive	
}

void Molybdenum::update_vapourpressure(){
	E_Debug("\tIn Molybdenum::update_vapourpressure()\n\n");
//	double AmbientPressure = 0;
	St.VapourPressure = probe_vapourpressure(St.Temperature);
	// Answer is in pascals
}

double Molybdenum::probe_vapourpressure(double Temperature)const{
	double VapourPressure(0.0);
	
	if( Temperature < Ec.MeltingTemp ){
		VapourPressure = 101325*pow(10,11.529 -34626/Temperature -1.1331*log10(Temperature));
	}else{
		VapourPressure = 101325*pow(10,11.529 -34626/Temperature -1.1331*log10(Temperature));
		static bool runOnce = true;
		WarnOnce(runOnce,
				"In Molybdenum::probe_vapourpressure():\nExtending model outside temperature range! ( from St.MeltingTemp > to St.BoilingTemp > )");
		std::cerr << "\nError! Negative Temperature in Molybdenum::probe_vapourpressure(double Temperature)";
	}
	return VapourPressure;
}