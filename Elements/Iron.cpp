//#define ELEMENT_DEEP_DEBUG
//#define ELEMENT_DEBUG
#include "Iron.h"
#include "Constants.h"

const struct ElementConsts Iron::IronConsts = {
	'F',			// Specifies the element
	1811, 			// K, Melting temperature at atmospheric pressure
	3134,			// K, Boiling temperature at atmospheric pressure
	13.810,			// Units of kJ/mol, Bond Energy, EQUAL TO LATENT HEAT
	4.5,			// ev, Work Function, http://hyperphysics.phy-astr.gsu.edu/hbase/tables/photoelec.html
	0.0057, 		// kJ/(m^2 K), HeatTransfer coefficient, 
					// http://www.engineeringtoolbox.com/overall-heat-transfer-coefficients-d_284.html
	0.055845,		// AMU in kg mol^-1, (+/- 0.00001) Atomic Mass
	13.810/0.055845, 	// kJ/kg, Latent Fusion Energy From Wikipedia
	340/0.055845, 		// kJ/kg, Latent Vapour Energy From Wikipedia
	3.29e16,		// rad s^-1, Plasma Frequency, Calculated from Robbies code with 4 electrons per atom
	1.862,			// N/m, Surface Tension B. Keene, 1993, Review of data for the surface tension of pure metals
	7874,			// (kg/m^3) from Wikipedia, Room temperature density
	0.0795			// kW/m K, Thermal Conductivity at 20 degrees C
				// http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/thrcn.html
};

Iron::Iron():Matter(&IronConsts){
	set_defaults();
}

Iron::Iron(double radius):Matter(radius,&IronConsts){
	set_defaults();
}

Iron::Iron(double radius, double tempin):Matter(radius,tempin,&IronConsts){
	E_Debug("\n\nIn Iron::Iron(double radius, double tempin)");
	set_defaults();

	update_state(0.0);		// Temperature dependant
	update_models('c','c','c','y');
	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

Iron::Iron(double radius, double tempin, std::array<char,4> &constmodels):Matter(radius,tempin,&IronConsts,constmodels){

	E_Debug("\n\nIn Iron::Iron(double radius, double tempin, std::array<char,4> &constmodels)");

	set_defaults();
	update_state(0.0);		// Temperature dependant
	update_models(constmodels);
	E_Debug("\nMass after = " << St.Mass << "\nRadius After = " << St.Radius << "\nSt.Density = " << St.Density 
			<< "\nSt.Volume = " << St.Volume);
}

void Iron::set_defaults(){
	E_Debug("\n\nIn Iron::set_defaults()");

        St.HeatCapacity = 0.1;	  		// J/(gK) or kJ/(kgK)
        St.Emissivity = 0.2;			// Arb, Emissivity
        St.SuperBoilingTemp =  Ec.BoilingTemp; 	// K, At any pressure
	St.Density = 7874;			// kg/m^3 From wikipedia
	update_models('c','c','c','y');


//	St.FractionalExpansion = 0.348; 	// Units of kg m^-3 K^-1: Calculated ROUGLY by measuring the gradient of Figure 1 found at https://aaltodoc.aalto.fi/bitstream/handle/123456789/14194/D4_tesfaye_fiseha_2010.pdf?sequence=1. See also https://www.jstage.jst.go.jp/article/matertrans1960/14/2/14_2_120/_pdf.  https://en.wikipedia.org/wiki/Iron#Mechanical_properties
//	Ec.BondEnergy = Ec.LatentFusion; 	// Units of kJ/mol, from http://ac.els-cdn.com/0009261490871992/1-s2.0-0009261490871992-main.pdf?_tid=15d8855a-d42a-11e6-bda3-00000aacb360&acdnat=1483718970_aad5e56f1c3219f79aec68fc1725b992, taken from figure 1b, TMBE-2000, then divided by AvNo
}

void Iron::update_heatcapacity(){ // Calculates the heat capacity in units of J mol^-1 K^-2
	E_Debug("\n\nIn Iron::update_heatcapacity()");
	double t = St.Temperature / 1000;
	E_Debug("\n St.Temperature is : " << St.Temperature << " K\n Gas=" << St.Gas << " \n Liquid = " << St.Liquid);

	if( St.Liquid == true ){
		if(St.Temperature < 2200){
			// Recommended values from: http://nist.gov/data/PDFfiles/jpcrd298.pdf
			St.HeatCapacity = 46.632;  // +/- 3
		}else{
			// From NIST: http://webbook.nist.gov/cgi/inchi?ID=C7439896&Mask=2			
			St.HeatCapacity = 46.02400 - 1.884667e-8*t + 6.094750e-9 *pow(t,2) - 6.640301e-10*pow(t,3) - (8.246121e-9/pow(t,2));
		}
	}else if( St.Gas == true ){
		
	}else{ // Must be a solid
		// All values (Except Temp < 298) from from:
		// http://webbook.nist.gov/cgi/cbook.cgi?ID=C7439896&Units=SI&Mask=2&Type=JANAFS&Plot=on#JANAFS
		if(St.Temperature <= 298){
			// HeatCapacity = ( 4.942*Temperature + (1943.75/pow(465,3))*pow(Temperature,4) ) / 1000;
			// This is an 8th order polynomial fit to the low temperature heat capacity data found at
			// http://nist.gov/data/PDFfiles/jpcrd298.pdf
			double z=(St.Temperature-120)/99;
			St.HeatCapacity = -0.71*pow(z,8) + 2.2*pow(z,7) + 0.38*pow(z,6) - 6.8*pow(z,5) 
				+ 5.3*pow(z,4) + 3.5*pow(z,3) - 8.7*pow(z,2) + 12*z + 15;
		}else if(St.Temperature > 298 && St.Temperature <= 700){ 		// 298 - 700	
			St.HeatCapacity = 18.42868 + 24.64301*t - 8.91372*pow(t,2) + 9.664706*pow(t,3) - (0.012643/pow(t,2));
		}else if(St.Temperature > 700 && St.Temperature <= 1042){	 	// 700 - 1042
			St.HeatCapacity = -57767.65 + 137919.7*t - 122773.2*pow(t,2) + 38682.42*pow(t,3) + (3993.080/pow(t,2));
		}else if(St.Temperature > 1042 && St.Temperature <= 1100){ 		// 1042 - 1100
			St.HeatCapacity = -325.8859 + 28.92876*t + (411.9629/pow(t,2));
		}else if(St.Temperature > 1100 && St.Temperature <= Ec.MeltingTemp ){	// 1100 - 1809
			St.HeatCapacity = -776.7387 + 919.4005*t - 383.7184*pow(t,2) + 57.08148*pow(t,3) + (242.1369/pow(t,2));	
		}else{
			St.HeatCapacity = 46.632;
		}
	}
//	E_Debug("\nCv of Solid: " << St.HeatCapacity << "[J/(mol K)]");
//	E_Debug("\nCv of Solid: " << St.HeatCapacity/(St.AtomicMass*1000) << "[kJ/(kg K)]");
//	Pause();
	St.HeatCapacity = (St.HeatCapacity/(Ec.AtomicMass*1000)); // Convert from J/(mol K) to kJ/( kg K ), Atomic Mass [kg mol^-1]
};

void Iron::update_radius(){
	E_Debug("\n\nIn Iron::update_radius()");

	St.LinearExpansion = 1+(10.8*St.Temperature)*1e-6;
	St.Radius=St.UnheatedRadius*St.LinearExpansion;
	if(St.Temperature == Ec.MeltingTemp){
		static bool runOnce;
		WarnOnce(runOnce,"Linear exansion Discontinuous in time");
	}
	E_Debug("\nTemperature = " << St.Temperature << "\n\nSt.LinearExpansion = " << St.LinearExpansion
			<< "\nSt.Radius = " << St.Radius);
	assert(St.Radius>0); // Assert radius is positive	
//	std::cerr << "Error! Diameter variation with temperature not specified!";
//	throw std::exception();
}

void Iron::update_vapourpressure(){
	E_Debug("\n\nIn Iron::update_vapourpressure()");
//	St.VapourPressure = pow(10,6.347 - 19574/St.Temperature); // http://mmrc.caltech.edu/PVD/manuals/Metals%20Vapor%20pressure.pdf

	// https://en.wikipedia.org/wiki/Vapor_pressures_of_the_elements_(data_page)
	if( !St.Liquid && !St.Gas ) St.VapourPressure = 0;//pow(10,12.106 - 21723/St.Temperature + 0.4536*log(St.Temperature) - 0.5846/pow(St.Temperature,3));
	if( St.Liquid ){
		static bool runOnce = true;
		WarnOnce(runOnce,"Temperature range of model extended from 2100K to 3134K");
		St.VapourPressure = pow(10,11.353 - 19574/St.Temperature);
	}
	// St.VapourPressure *= 9.86923e-6; // Conversion from pascals to atmospheres. WE WANT IT IN PASCALS IDIAT
	E_Debug("VapourPressure = " << St.VapourPressure);
}
