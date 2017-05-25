//#define ELEMENT_DEBUG
//#define ELEMENT_DEEP_DEBUG
#include "Beryllium.h"
#include "Constants.h"

const struct ElementConsts Beryllium::BerylliumConsts = {
	'B',			// Specifies the element
	1560, 			// K, Melting temperature at atmospheric pressure
	2742,			// K, Boiling temperature at atmospheric pressure
	320.3,	 		// Units of kJ/mol, Bond Energy http://www.periodni.com/be.html
	3.4,			// ev, Work Function
				// goodfellow//5.0;//from http://hyperphysics.phy-astr.gsu.edu/hbase/tables/photoelec.html
	5.7, 			// HeatTransfer coefficient 
				// http://www.engineeringtoolbox.com/overall-heat-transfer-coefficients-d_284.html
	0.009012182,		// AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
	12.2/0.009012182, 	// kJ/kg, Latent Fusion Energy From Wikipedia
	292/0.009012182, 	// kJ/kg, Latent Vapour Energy From Wikipedia
	3.093e14,		// rad s^-1, Plasma Frequency
	1.100,			// N/m, Surface Tension B. Keene, 1993, Review of data for the surface tension of pure metal
	1840,			// (kg/m^3) from Wikipedia, Room temperature density
	0.200,			// kW/mK, from Wikipedia,
};


Beryllium::Beryllium():Matter(&BerylliumConsts){
	set_defaults();
}

Beryllium::Beryllium(double radius):Matter(radius,&BerylliumConsts){
	set_defaults();
}

Beryllium::Beryllium(double radius, double tempin):Matter(radius,tempin,&BerylliumConsts){
	E_Debug("\n\nIn Beryllium::Beryllium(double radius, double tempin)");
	set_defaults();

	update_state(0.0);		// Temperature dependant
	update_models('c','c','c','y');

	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Beryllium::Beryllium(double radius, double tempin, std::array<char,4> &constmodels):Matter(radius,tempin,&BerylliumConsts){
	E_Debug("\n\nIn Beryllium::Beryllium(double radius, double tempin)");
	set_defaults();

	update_state(0.0);		// Temperature dependant
	update_models(constmodels);

	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

void Beryllium::set_defaults(){
	E_Debug("\n\nIn Beryllium::set_defaults()");
	// http://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.html
        St.HeatCapacity = 1.8254448;	 	// kJ/(kg-K)         
	St.Emissivity = 0.18; 			// Arb, Polished Beryllium
        St.SuperBoilingTemp = Ec.BoilingTemp; 	// K, At any pressure
	St.Density = 1840; 			// (kg/m^3) from Wikipedia
	update_models('c','c','c','y');
}

/// ************************************************************************************************* \\\

void Beryllium::update_heatcapacity(){ // Calculates the heat capacity in units of J mol^-1 K^-2
	E_Debug("\n\nIn Beryllium::update_heatcapacity()");

	double t = St.Temperature/1000;
	if( St.Temperature > 250 && St.Temperature <= 298 ){ 

		static bool runOnce = true;
		WarnOnce(runOnce,
				"In Beryllium::update_heatcapacity():\nExtending model outside range! from T > 298 to T > 250");	
		// http://webbook.nist.gov/cgi/inchi?ID=C7440417&Mask=2
		St.HeatCapacity = 21.20694+5.688190*t+0.968019*pow(t,2)-0.001749*pow(t,3)-0.587526/(pow(t,2));
	}else if( St.Temperature > 298 && St.Temperature <= 1527 ){ 
		// http://webbook.nist.gov/cgi/inchi?ID=C7440417&Mask=2
		St.HeatCapacity = 21.20694+5.688190*t+0.968019*pow(t,2)-0.001749*pow(t,3)-0.587526/(pow(t,2));
	}else if( St.Temperature > 1527 && St.Temperature <= Ec.MeltingTemp ){
 		// http://webbook.nist.gov/cgi/inchi?ID=C7440417&Mask=2
		St.HeatCapacity = 30.00037-0.000396*t+0.000169*pow(t,2)-0.000026*pow(t,3)-0.000105/(pow(t,2));
	}else if( St.Temperature > Ec.MeltingTemp && St.Temperature <= St.SuperBoilingTemp ){
 		// http://webbook.nist.gov/cgi/inchi?ID=C7440417&Mask=2
		St.HeatCapacity = 25.42516+2.157953*t-0.002573*pow(t,2)+0.000287*pow(t,3)+0.003958/(pow(t,2));
//		std::cout << "\nHey, HeatCapacity = " << St.HeatCapacity; std::cin.get();
	}

	E_Debug("\n\nTemperature is : " << St.Temperature << "\nSt.Gas = " << St.Gas << "\nSt.Liquid = " << St.Liquid 
				<< "\nCv of Solid: " << St.HeatCapacity/Ec.AtomicMass << "[kJ/(kg K)]"; );
	St.HeatCapacity = (St.HeatCapacity /(1000 * Ec.AtomicMass)); // Conversion J/(mol K) to kJ/( kg K ), AtomicMass [kg mol^-1]
}

void Beryllium::update_radius(){
	E_Debug("\n\nIn Beryllium::update_radius():");
	if( St.Temperature > 250 && St.Temperature <= 298 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"In Beryllium::update_radius():\nExtending model outside range!(from T<298K to T<250K)");
		// www-ferp.ucsd.edu/LIB/PROPS/PANOS/be.html
		St.LinearExpansion = 1+St.Temperature*1e-6*(8.4305+1.1464e-2*St.Temperature+2.9752e-6*pow(St.Temperature,2));
	}else if( St.Temperature > 298 && St.Temperature < Ec.MeltingTemp ){
		// www-ferp.ucsd.edu/LIB/PROPS/PANOS/be.html
		St.LinearExpansion = 1+St.Temperature*1e-6*(8.4305+1.1464e-2*St.Temperature+2.9752e-6*pow(St.Temperature,2));
	}else if( St.Temperature >= Ec.MeltingTemp && St.Temperature <= St.SuperBoilingTemp){ 
		// This could do with fixing at some point...
		static bool runOncetwo = true;
		WarnOnce(runOncetwo,"In Beryllium::update_radius():\nExtending model outside range!(from T<Tmelt to T<Tboil)");
		// www-ferp.ucsd.edu/LIB/PROPS/PANOS/be.html
		St.LinearExpansion = 1+St.Temperature*1e-6*(8.4305+1.1464e-2*St.Temperature+2.9752e-6*pow(St.Temperature,2));	
	}
	St.Radius=St.UnheatedRadius*St.LinearExpansion;
	E_Debug("\nTemperature = " << St.Temperature << "\n\nSt.LinearExpansion = " << St.LinearExpansion
			<< "\nSt.Radius = " << St.Radius);
	assert(St.Radius>0); // Assert radius is positive	
}

void Beryllium::update_vapourpressure(){

	// Also available from wikipedia in pascal: https://en.wikipedia.org/wiki/Vapor_pressures_of_the_elements_(data_page)
	// St.VapourPressure = (pow(10,13.048-17020/St.Temperature - 0.4440*log10(St.Temperature)))/101325;

	if( !St.Liquid && !St.Gas ){
		// http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
		St.VapourPressure = pow(10,8.042 -17020/St.Temperature-0.444*log10(St.Temperature));
	}else if( St.Liquid ){
		// http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf
		St.VapourPressure = pow(10,5.786 -15731/St.Temperature);
	}else{
		St.VapourPressure = 0;
		std::cout << "\nWarning, Sample is assumed gas! St.VapourPressure = 0";
	}
	St.VapourPressure = 101325*St.VapourPressure; // Convert to Pascals
}
