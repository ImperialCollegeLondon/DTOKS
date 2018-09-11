//#define ELEMENT_DEBUG
//#define ELEMENT_DEEP_DEBUG
#include "Graphite.h"
#include "Constants.h"

const struct ElementConsts GraphiteConsts = {
	'G',			// Specifies the element
	4500, 			// K, Melting temperature at atmospheric pressure
	4000,			// K, Boiling temperature at atmospheric pressure
	345,			// Units of kJ/mol, Bond Energy, (+/- 1) 
				// Chapter 2. Carbon (Graphene/Graphite), Springer.
	4.80,			// ev, Work Function, taken from DTOKS
	100, 			// kJ/(m^2 K), HeatTransfer coefficient, (ROUGH ESTIMATE)
				// http://cr4.globalspec.com/thread/72542/Overall-heat-Transfer-Coefficient-for-a-Graphite-Condenser
	0.0120107,		// AMU in kg mol^-1, http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782425&Mask=2
	0,		 	// kJ/kg, Latent Fusion Energy From Wikipedia
	712.912/0.0120107,	// kJ/kg, Latent Vapour Energy, http://aip.scitation.org/doi/10.1063/1.1746999 (+/-0.8)
	0,			// N/m, Surface Tension, Is never a liquid
	2260,			// (kg/m^3) from Wikipedia, Room temperature density (+/- 100)
	0.140			// kW/m K, Thermal Conductivity at 293K, (+/- 0.001)
				// http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.736.9349&rep=rep1&type=pdf
};


Graphite::Graphite():Matter(GraphiteConsts){
	set_defaults();
	update();
}

Graphite::Graphite(double radius):Matter(radius,GraphiteConsts){
	set_defaults();
	update();
}

Graphite::Graphite(double radius, double tempin):Matter(radius,tempin,GraphiteConsts){
	E_Debug("\n\nIn Graphite::Graphite(double radius, double tempin)");
	set_defaults();

	update_state(0.0);		// Temperature dependent
	update_models('c','c','c','y','n');
	update();
	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Graphite::Graphite(double radius, double tempin, std::array<char,CM> &constmodels):Matter(radius,tempin,GraphiteConsts){
	E_Debug("\n\nIn Graphite::Graphite(double radius, double tempin, std::array<char,CM> &constmodels)");
	set_defaults();

	update_state(0.0);		// Temperature dependent
	update_models(constmodels);
	update();

	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Graphite::Graphite(double radius, double tempin, std::array<char,CM> &constmodels, const threevector & position, const threevector& velocity):Matter(radius,tempin,GraphiteConsts){
	E_Debug("\n\nIn Graphite::Graphite(double radius, double tempin, std::array<char,CM> &constmodels, const threevector & position, const threevector& velocity)");
	set_defaults();

	update_state(0.0);		// Temperature dependent
	update_models(constmodels);
	update_motion(position,velocity,0.0);
	update();

	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

void Graphite::set_defaults(){
	E_Debug("\n\nIn Graphite::set_defaults()");
        St.HeatCapacity = 0.846; 		// kJ/(kg-K)			(+/- 0.001)
	// http://www-eng.lbl.gov/~dw/projects/DW4229_LHC_detector_analysis/calculations/emissivity2.pdf
	St.Emissivity = 0.70; 			// Arb, Polished Graphite, 	(+/- 0.1)
        St.SuperBoilingTemp = Ec.BoilingTemp;	// K, At any pressure,
	St.Density = 2260; 			// (kg/m^3)			(+/- 100)
	update_models('c','c','c','y','n');
}


/// ************************************************************************************************* \\\

void Graphite::update_heatcapacity(){ // Calculates the heat capacity in units of J mol^-1 K^-2
	E_Debug("\n\nIn Graphite::update_heatcapacity()");
	
	if( St.Temperature > 200 && St.Temperature <= 3500){
		// Redirected from NIST: http://webbook.nist.gov/cgi/formula?ID=C7782425&Mask=2
		// http://ac.els-cdn.com/0022311573900603/1-s2.0-0022311573900603-main.pdf?_tid=3478dcf2-0338-11e7-a73c-00000aacb361&acdnat=1488892739_9595b510ccead6c7baa0f96a11b9d39d
		St.HeatCapacity = 0.54212-2.42667e-6*St.Temperature-90.2725/St.Temperature-43449.3/(pow(St.Temperature,2))
					+1.59309e7/(pow(St.Temperature,3))-1.43688e9/(pow(St.Temperature,4));

		St.HeatCapacity = St.HeatCapacity*4.184; // Convert from calorie/gram to KiloJoule / Kilogramme
	}else if( St.Temperature > 3500 && St.Temperature <= 4000){
		static bool runOnce = true;
		WarnOnce(runOnce,
				"In Graphite::update_heatcapacity():\nExtending model outside range!(from T < 3500K to T < 4000K)");

		St.HeatCapacity = 0.54212-2.42667e-6*St.Temperature-90.2725/St.Temperature-43449.3/(pow(St.Temperature,2))
					+1.59309e7/(pow(St.Temperature,3))-1.43688e9/(pow(St.Temperature,4));

		St.HeatCapacity = St.HeatCapacity*4.184; // Convert from calorie/gram to KiloJoule / Kilogramme	
	}
	E_Debug("\nTemperature is : " << St.Temperature << "\nSt.Gas = " << St.Gas << "\nSt.Liquid = " << St.Liquid 
				<< "\nCv of Solid: " << St.HeatCapacity/Ec.AtomicMass << "[kJ/(kg K)]"; );
		// http://webbook.nist.gov/cgi/cbook.cgi?ID=C7440440&Type=JANAFG&Plot=on
}

void Graphite::update_radius(){
	E_Debug("\n\nIn Graphite::update_radius():");
/*	if( St.Temperature > 273 && St.Temperature <= 1800 ){
		// Day and sosman (1912) 
		St.LinearExpansion = 1+0.5*1e-6+3.2*1e-9*(St.Temperature-273);
	}*/
	assert(St.Temperature > 273);
	if(St.Temperature > 273 && St.Temperature < 4000){ // Old limits: 1800 : 2500
		// http://aip.scitation.org/doi/pdf/10.1063/1.2945915
		// Accurate to 2.5%, Temperature in degrees Celsius.
		//St.LinearExpansion = 1-5.27e-2 + 6.98e-4*(St.Temperature-273)+7.76e-8*pow(St.Temperature-273,2);
		St.LinearExpansion = 1-5.27e-2 + 6.98e-4*St.Temperature+7.76e-8*pow(St.Temperature,2);
	}

	// Assert change in radius is less than 1% of current radius.
	assert(abs(St.Radius-St.UnheatedRadius*St.LinearExpansion)<(St.Radius*0.01));
 	St.Radius=St.UnheatedRadius*St.LinearExpansion;
	E_Debug("\nTemperature = " << St.Temperature << "\n\nSt.LinearExpansion = " << St.LinearExpansion << "\nSt.Radius = " 
			<< St.Radius);
	assert(St.Radius>0); // Assert radius is positive	
}

void Graphite::update_vapourpressure(){
//	double AmbientPressure = 0;
	// http://pubs.acs.org/doi/pdf/10.1021/ja01161a081

	// The factor alpha is the accommodation coefficient and represents the ratio of the 
	// rate at which atoms actually condense on the surface to the rate at which they strike the surface
	//double alpha = 1; 
	//St.VapourPressure = pow(10,log10((St.AtomicMass*1000*alpha)/(2*PI*83.15*1e6*St.Temperature))+0.5*log10(St.Temperature)-2.187); 
	St.VapourPressure = 0; // Carbon is never in the liquid state
}
