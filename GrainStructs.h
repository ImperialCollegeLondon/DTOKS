// Structure containing all the intrinsic and extrinsic information about a material
// Comprised of a particular element

#ifndef __GRAINSTRUCTS_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __GRAINSTRUCTS_H_INCLUDED__

// Data that varies with time and is generally a result of the simulation
struct GrainData{
	// Define if the material is a liquid or a gas
	bool Liquid, Gas;

	// Do vary with Temperature
	double SuperBoilingTemp;	// Kelvin, Super heated boiling temperature
	double UnheatedRadius;		// metres
	double Mass; 			// Kilogrammes,
	double Radius;			// metres
	double SurfaceArea;		// metres^2
	double Volume;			// metres^3
	double Density;			// Kilogrammes / (metre^3)
	double Temperature; 		// Kelvin
	double VapourPressure; 		// N/m^2 or Pa, Pascals
	double Emissivity;		// Arb, deviation of EM radiation from black body spectrum
	double LinearExpansion; 	// Units of m K^-1
	double HeatCapacity;		// kJ/(kg·K)	(Constant Pressure!)
	// Latent Heat
	double FusionEnergy;		// kJ, Energy put into the material to break bonds of solid
	double VapourEnergy;		// kJ, Energy put into the material to break bond of liquid

};

// Input parameters which are specific to the type of matter being tested
struct ElementConsts{
	char Elem;		// Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.

	// Don't vary with Temperature
	double MeltingTemp;	// Kelvin
	double BoilingTemp;	// Kelvin
	double BondEnergy; 	// kJ/mol Amount of energy required to break a single atomic bond
	double WorkFunction;	// kJ, Amount of energy required to remove an electron
	double HeatTransAir; 	// W/m^2 K 
	double AtomicMass;	// kg/mol
	double LatentFusion;	// kJ/kg, Energy required to melt 1kg of solid
	double LatentVapour; 	// kJ/kg, Energy required to vapourize 1kg of liquid
	double PlasmaFrequency;	// rad s^-1, plasma frequency of electrons in the material
	double SurfaceTension;	// N/m,
	double RTDensity;	// Kilogrammes / (metre^3), Denisty at room temperature
	double ThermConduct; 	// KiloWatts / (Metre * Kelvin) 

};

#endif
