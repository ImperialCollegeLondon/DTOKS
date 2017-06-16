//#define PAUSE
//#define MATTER_DEEP_DEBUG
//#define MATTER_DEBUG

#include "Matter.h"
struct GrainData MatterDefaults = {
	false,		// Is Liquid,		fine
	false,		// Is Gas,		fine
	1e-6,		// Unheated Radius,	should be updated
	8.21e-14,	// Mass,		should be updated
	1e-6,		// Radius,		fine
	1.26e-11,	// Surface Area		should be updated
	4.189e-18,	// Volume		should be updated
	19600,		// Density		should be updated
	3000,		// SuperBoilingTemp	should be updated
	270,		// Temperature		fine
	1e5,		// Vapour Pressure	should be updated
	1.0,		// Emissivity		should be updated
	1.0,		// Linear Expansion	should be updated
	0.5,		// Heat Capacity	should be updated
	0,		// DeltaSec		should be updated
	0,		// DeltaTherm		should be updated
	0,		// Potential		should be updated
	false,		// Is Positive		should be updated
	{0,0,0},	// Dust position	case dependant
	{0,0,0},	// Dust Velocity	case dependant
	0,		// FusionEnergy		fine
	0		// VapourEnergy		fine
};


// Constructors
Matter::Matter(const ElementConsts *elementconsts):Ec(*elementconsts),St(MatterDefaults){
	M_Debug("\n\nIn Matter::Matter(const ElementConsts *elementconsts):Ec(*elementconsts),St(MatterDefaults)\n\n");
};
Matter::Matter(double rad, const ElementConsts *elementconsts):Ec(*elementconsts),St(MatterDefaults){
	M_Debug("\n\nIn Matter::Matter(double rad, const ElementConsts *elementconsts):Ec(*elementconsts),St(MatterDefaults)\n\n");
	St.Radius = rad;			// m
	St.UnheatedRadius = St.Radius;		// m
	assert(St.Radius > 0 && St.UnheatedRadius > 0);
//	M_Debug("\nSt.Radius = " << St.Radius << "\nSt.Mass = " << St.Mass);
};

Matter::Matter(double rad, double temp, const ElementConsts *elementconsts):Ec(*elementconsts),St(MatterDefaults){
	M_Debug("\n\nIn Matter::Matter(double rad, double temp, const ElementConsts *elementconsts):Ec(*elementconsts),St(MatterDefaults)\n\n");
	ConstModels = {'c','c','c','y'};
	St.Radius = rad;					// m
	St.UnheatedRadius = St.Radius;				// m
	St.Temperature = temp;					// K
	assert(St.Radius > 0 && St.UnheatedRadius > 0 && St.Temperature > 0);

//	M_Debug("\nSt.Radius = " << St.Radius << "\nSt.Mass = " << St.Mass << "\nSt.Temperature = " << St.Temperature);
//	M_Debug("\nEc.LatentVapour = " << Ec.LatentVapour << "\nEc.AtomicMass = " << Ec.AtomicMass);
};
Matter::Matter(double rad, double temp, const ElementConsts *elementconsts, std::array<char,4> &constmodels):Ec(*elementconsts),St(MatterDefaults){
	M_Debug("\n\nIn Matter::Matter(double rad, double temp, const ElementConsts *elementconsts, std::array<char,4> &constmodels):Ec(*elementconsts),St(MatterDefaults)\n\n");
	ConstModels = constmodels;
	St.Radius = rad;					// m
	St.UnheatedRadius = St.Radius;				// m
	St.Temperature = temp;					// K
	assert(St.Radius > 0 && St.UnheatedRadius > 0 && St.Temperature > 0);

//	M_Debug("\nSt.Radius = " << St.Radius << "\nSt.Mass = " << St.Mass << "\nSt.Temperature = " << St.Temperature);
//	M_Debug("\nEc.LatentVapour = " << Ec.LatentVapour << "\nEc.AtomicMass = " << Ec.AtomicMass);
};

// Update geometric properties of matter:
// Volume, Surface area and Density. Also initialises mass
void Matter::update_dim(){ // Assuming spherical particle.		
	M_Debug("\tIn Matter::update_dim():\n\n");

	if(ConstModels[1] == 'v' || ConstModels[1] == 'V'){ // THIS WHOLE PART NEEDS CHECKING...
		update_radius();
		St.Density = Ec.RTDensity / pow(St.LinearExpansion,3); // SHOULD THIS BE : pow(2*St.LinearExpansion,3) ???
		M2_Debug("\nUnheatedRadius was = " << St.UnheatedRadius);
		// Change the radius according to the amount of mass lost
		St.UnheatedRadius = St.UnheatedRadius*(pow((3*St.Mass)/(4*PI*St.Density),1./3.)/St.Radius);
		M2_Debug("\nUnheatedRadius now = " << St.UnheatedRadius);

	}else if(ConstModels[1] == 'c' || ConstModels[1] == 'C'){
//		St.Radius = St.UnheatedRadius; // THIS IS SO WRONG. RADIUS SHOULD DECREASE AS MASS IS LOST
		St.Density = Ec.RTDensity;				// Density is RT density
	}else if(ConstModels[1] == 's' || ConstModels[1] == 's'){
		St.Density = Ec.RTDensity;				// Fix the density
	}else{
		std::cout << "\nError! In Matter::update_dim(char ConstModels[1])\n"
                                << "Invalid input for Expansion Model.";
		assert( strchr("vVcCsS",ConstModels[1]) );
	}
	St.Volume = (4*PI*pow(St.Radius,3))/3;
	St.SurfaceArea = 4*PI*pow(St.Radius,2);
	// This is a weird way of updating mass because:
	// When mass is lost, it is removed from St.Mass, then the new 'St.UnheatedRadius' is calculated.
	// This is then used to calculate the new smaller St.Volume. Finally the St.Mass is calculated again.
	static bool InitialiseMass = true;
	if( InitialiseMass ){
	        St.Mass = St.Density*St.Volume; 
	        InitialiseMass = false;
	}

//	assert( abs((St.Density - St.Mass/St.Volume)/St.Density) < 0.000001 ); // Sanity Check, this may be an issue
//	St.Mass = St.Density*St.Volume;
	
//	std::cout << "\n = " << St.Mass; 
//	std::cout << "\nSt.Mass = " << St.Mass; std::cin.get();
	
	M2_Debug("\nSt.Mass = " << St.Mass << "\nSt.Volume = " << St.Volume << "\nSt.Density = " << St.Density
			<< "\nSt.Mass/St.Volume = " << St.Mass/St.Volume);

	assert(St.Mass > 10e-25 );
}

void Matter::update_emissivity(){
	M_Debug("\tIn Matter::update_emissivity()\n\n");
	if( ConstModels[0] == 'c' || ConstModels[0] == 'C' ){ // Assume constant emissivity
		St.Emissivity = St.Emissivity;
	}else if( ConstModels[0] == 'f' || ConstModels[0] == 'F' ){ // Get emissivity from file
		if( St.Radius < 2e-8 || St.Radius > 1e-4 ){
			std::cout << "\nError! In Matter::update_emissivity().\n"
					<< "Radius = " << St.Radius << " outside limit of Emissivty model.\n";
			assert(St.Radius > 2e-8);
			assert(St.Radius < 1e-4);
		}
		assert(St.Temperature > 275);
		std::string DirName, FileName("/Temp_"), TempString, txtstring(".txt");
		if(Ec.Elem == 'w' || Ec.Elem == 'W' ) DirName = "EmissivityData/EmissivityDataTungsten";
		else if(Ec.Elem == 'f' || Ec.Elem == 'F' ) DirName = "EmissivityData/EmissivityDataIron";
		else if(Ec.Elem == 'g' || Ec.Elem == 'G' ) DirName = "EmissivityData/EmissivityDataGraphite";
		else if(Ec.Elem == 'b' || Ec.Elem == 'B' ) DirName = "EmissivityData/EmissivityDataBeryllium";
		std::ifstream TempFile;
		std::stringstream ss;

		ss << round(St.Temperature);
		ss >> TempString;
		TempFile.open((DirName+FileName+TempString+txtstring).c_str());
		bool Found = false;

		char buffer;
		double TestRad;
		while( TempFile >> TestRad >> buffer >> St.Emissivity && Found != true ){
			if(TestRad == round_to_digits(St.Radius,2)){ 
				Found = true;
				M_Debug( "\nround(St.Radius,9) = " << round_to_digits(St.Radius,2) 
					<< "\nTestRad = " << TestRad );
				//Pause();
			}
		}
		ss.clear(); ss.str("");
		TempFile.close();
	}else{ // Invalid input for Emissivity model
		std::cout << "\nError! In Matter::update_emissivity()\n"
				<< "Invalid input for Emissivity Model.";
		assert( strchr("cCvVfF",ConstModels[1]) );
	}

	if(St.Emissivity > 1.0){
		static bool runOnce = true;
		WarnOnce(runOnce,"Emissivity > 1! Emissivity being forced equal to 1");
		St.Emissivity = 1.0;
	}
}

void Matter::update_boilingtemp(){
	M_Debug("\tIn Matter::update_boilingtemp()\n\n");
	// Conver LatentVapour from kJ to J for the Gas constant
//	St.CapillaryPressure = exp((St.LatentVapour*St.AtomicMass)/(1000*R)*(1/St.Temperature() - 1/(St.MeltingTemp)));
	double CapillaryPressure = Ec.SurfaceTension*2/(St.Radius*101325);
	M2_Debug("\n\nCapillaryPressure = " << CapillaryPressure << "atmospheres.");
	if( ConstModels[3] == 's' || ConstModels[3] == 'S' ){ // For super boiling
		assert(Ec.Elem != 'g' && Ec.Elem != 'G'); // We can't do Graphite, it has no vapour pressure
		St.SuperBoilingTemp = 1/(1/Ec.BoilingTemp 
					- R*log(CapillaryPressure)/(1000*Ec.LatentVapour*Ec.AtomicMass));
		if(Ec.BoilingTemp > St.SuperBoilingTemp){
			St.SuperBoilingTemp = Ec.BoilingTemp;
		}
	}else if( ConstModels[3] == 't' || ConstModels[3] == 'T' ){ // For Thomson Boiling
		assert(Ec.Elem != 'g' && Ec.Elem != 'G'); // We can't do Graphite, it has no vapour pressure
		// Derived from Kelvin and Clausius-Clapyeron Equation, makes a difference of ~5% for
		// Grains of order 1e-8m.
		St.SuperBoilingTemp = Ec.BoilingTemp*exp(-(2*Ec.SurfaceTension)/
			(St.Radius*Ec.LatentVapour*1000*St.Density));
		// Kelvin Equation, Makes a difference of ~ 10% for Grains of size 1e-8m,
		// This is an increasingly small effect for larger grains
		St.VapourPressure = St.VapourPressure*exp((2*Ec.SurfaceTension*Ec.AtomicMass)
								/(St.Density*St.Radius*R*St.Temperature));
	}else if( ConstModels[3] == 'y' || ConstModels[3] == 'Y' ) St.SuperBoilingTemp = Ec.BoilingTemp;
	else if( ConstModels[3] == 'n' || ConstModels[3] == 'N' ) St.SuperBoilingTemp = std::numeric_limits<double>::max();
	else{
		std::cout << "\nError! In Matter::update_boilingtemp()\n"
				<< "Invalid input '" << ConstModels[3] << "' for Emissivity Model.";
		assert( strchr("yYnNsStT",ConstModels[1]) );
	}
	M2_Debug("\nR*log(CapillaryPressure)/(1000*Ec.LatentVapour*Ec.AtomicMass) = " 
		<< R*log(CapillaryPressure)/(1000*Ec.LatentVapour*Ec.AtomicMass) << "\nSt.SuperBoilingTemp = " 
		<< St.SuperBoilingTemp);
}

void Matter::update_state(double EnergyIn){
	M_Debug("\tIn Matter::update_state(double EnergyIn)\n\n");

	if( St.Temperature < St.SuperBoilingTemp && St.Temperature >= Ec.MeltingTemp){ // Melting or Liquid
		if( St.Temperature > Ec.BoilingTemp){
			static bool runOnce = true;
			WarnOnce(runOnce,"Substance is superheated! Temperature > BoilingTemp");
		}
		if( St.FusionEnergy < Ec.LatentFusion*St.Mass ){ // Must be melting
			M2_Debug("\n*Liquid is melting*, T = " << St.Temperature
				<< "\nFusionEnergy is < LatentFusion*Mass : " << St.FusionEnergy << " < " 
					<< Ec.LatentFusion*St.Mass);
			// Add energy to Latent heat and set Temperature to Melting Temperature
			St.FusionEnergy += EnergyIn; 
			St.Temperature = Ec.MeltingTemp;
			if( St.FusionEnergy > Ec.LatentFusion*St.Mass ){ // if it melts fully
				St.Liquid = true; St.Gas = false;
				St.Temperature = Ec.MeltingTemp + 
				        (St.FusionEnergy-Ec.LatentFusion*St.Mass)/St.HeatCapacity;
			}

		}else{ St.Liquid = true; St.Gas = false; }  // Else it has melted!
	}else if( St.Temperature >= St.SuperBoilingTemp ){ // Boiling or Gas
		if( St.VapourEnergy < (Ec.LatentVapour*Ec.RTDensity*PI*pow(St.UnheatedRadius,3.0)*4.0)/3.0 ){ // Must be boiling
			// Add energy to Latent heat and set Temperature to Melting Temperature
//			M1_Debug( "\nEnergyIn = " << EnergyIn << "\nEnergyIn/Ec.LatentVapour = " << EnergyIn/Ec.LatentVapour << "\nMass = " << St.Mass << "\nEc.RTDensity*PI*pow(St.UnheatedRadius,3)*4/3 = "  <<  Ec.RTDensity*PI*pow(St.UnheatedRadius,3)*4/3 << "\nEv.LatentVapour = " << Ec.LatentVapour << "\nSt.VapourEnergy = " << St.VapourEnergy);
			update_mass(EnergyIn/Ec.LatentVapour);
			St.VapourEnergy += EnergyIn; 
			St.Temperature = St.SuperBoilingTemp;

			M2_Debug("\n*Liquid is boiling*, T = " << St.Temperature 
				<< " K\nVapourEnergy is < St.LatentVapour*Mass : " << St.VapourEnergy 
				<< " < " << Ec.LatentVapour*St.Mass);
			if( St.VapourEnergy > (Ec.LatentVapour*Ec.RTDensity*PI*pow(St.UnheatedRadius,3.0)*4.0)/3.0 ){
				St.Liquid = false; St.Gas = true; 
	                        St.Temperature = St.SuperBoilingTemp + 
	                                (St.VapourEnergy-Ec.LatentVapour*St.Mass)/St.HeatCapacity;
	                        std::cout << "\n\n***** Sample has Boiled THIS STEP! *****\n";
			}
		}else{	// Else it has vapourised!
			St.Liquid = false; St.Gas = true; 
			std::cout << "\n\n***** WARNING! SAMPLE IS GASEOUS *****\n";
		}
	}
	M2_Debug("\nSt.Temperature = " << St.Temperature);
	if( St.Mass < 10e-24 ){ // Lower limit for mass of dust
		std::cout << "\nSt.Mass = " << St.Mass << " < 10e-24, vaporisation assumed.";
		St.Gas = true;
	}
}

// Change geometric and variable properties of matter; 
// volume, surface area, density, heat capacity, emissivity and vapour pressure
void Matter::update(){
	M_Debug("\tIn Matter::update()\n\n");
	// ORDER DEPENDENT UPDATES:
	M_Debug("\t"); update_dim();
	
	// ORDER INDEPENDENT UPDATES:
	if( ConstModels[2] == 'v'|| ConstModels[2] == 'V') 	update_heatcapacity();
	else if( ConstModels[2] == 'c'||ConstModels[2] == 'C') 	St.HeatCapacity = St.HeatCapacity;
	else if( ConstModels[2] == 's'|| ConstModels[2] == 'S')	St.HeatCapacity = 0.56; // This is an arbitrary fixed number
	else{
		std::cout << "\nError! In Matter::update(char EmissivModel, char ExpansionModel, char ConstModels[0])\n"
                                << "Invalid input for HeatCapacity Model.";
                assert( strchr("vVcCsS",ConstModels[1]) );	
	}

	M_Debug("\t"); update_emissivity();
	M_Debug("\t"); update_vapourpressure();
	M_Debug("\t"); update_boilingtemp();
};

// Change geometric and variable properties of matter; 
// volume, surface area, density, heat capacity, emissivity and vapour pressure

void Matter::update_models(char emissivmodel, char linexpanmodel, char heatcapacitymodel, char boilingmodel){
        M_Debug("\tIn Matter::update_models(char emissivmodel, char linexpanmodel, char heatcapacitymodel, char boilingmodel)\n\n");
	ConstModels[0] = emissivmodel;
	ConstModels[1] = linexpanmodel;
	ConstModels[2] = heatcapacitymodel;
	ConstModels[3] = boilingmodel;
};

void Matter::update_models(std::array<char,4> &constmodels){
	M_Debug("\tIn Matter::update_models(std::array<char,4> &constmodels)\n\n");
//	update(constmodels[0],constmodels[1],constmodels[2],constmodels[3],);
	ConstModels = constmodels;
}

// Mass lost in Kilogrammes
void Matter::update_mass(double LostMass){ 
	M_Debug("\tIn Matter::update_mass(double LostMass)\n\n");
	M2_Debug("\nMassLoss = " << LostMass  << "\nSt.Mass = " << St.Mass << "\nRadius = " << St.Radius 
			<< "\nDensity = " << St.Density);
	St.Mass -= LostMass;
		
	//Pause();
	if(St.Mass < 0){ 
		std::cout << "\n\nSt.Mass = " << St.Mass << "\nSample mass is Negative! Evaporation assumed\n\n";
		St.Mass = 0;
		St.Gas = true;
		St.Liquid = false;
		//throw std::exception();
	}
}

// Takes argument of amount of energy lost in Kilo Joules and changes temperature in Degrees
void Matter::update_temperature(double EnergyIn){

	M_Debug("\tIn Matter::update_temperature(double EnergyIn = " << EnergyIn << " kJ)\n\n");
//	M_Debug("\n...\nTemp = " << St.Temperature << "K \nHeatCapacity = " << St.HeatCapacity 
//		<< " kJ/(kg K)\nFusionEnergy = " << St.FusionEnergy << " kJ\nVapourEnergy = " 
//		<< St.VapourEnergy << " kJ\n\n");

	// Assert temperature changes smaller than 10% of current temperature
	if( abs(EnergyIn/(St.Mass*St.HeatCapacity)) > St.Temperature*0.1 ){
		M2_Debug("\n\nSt.Temperature = " << St.Temperature << "\nSt.Mass = " << St.Mass << "\nSt.HeatCapacity = " 
			<< St.HeatCapacity << "\nEnergyIn = " << EnergyIn);
	}

	if(abs(EnergyIn/(St.Mass*St.HeatCapacity)) > St.Temperature*0.1){
		std::cout << "\nEnergyIn = " << EnergyIn << "\nSt.Mass = " << St.Mass << "\nSt.HeatCapacity = " << St.HeatCapacity;
		std::cout << "\nError, EnergyIn/(St.Mass*St.HeatCapacity) = " << EnergyIn/(St.Mass*St.HeatCapacity); // Assert temperature change is less than 10%
		std::cout << "\nTemperature = " << St.Temperature;
	}
	assert(abs(EnergyIn/(St.Mass*St.HeatCapacity)) < St.Temperature*0.1); // Assert temperature change is less than 10%

	St.Temperature += EnergyIn/(St.Mass*St.HeatCapacity); // HeatCapacity in Units of kJ/(kg K)
	M2_Debug( "\n**** Temp change = " << EnergyIn/(St.Mass*St.HeatCapacity) << "\n**** Ein = " << EnergyIn 
			<< "\n**** Mass = " << St.Mass);

	assert(St.Temperature > 0);
	M_Debug("\t"); update_state(EnergyIn);
	M2_Debug("\nSt.Temperature = " << St.Temperature);
	if(!St.Gas)	assert(St.Temperature <= St.SuperBoilingTemp);
};

void Matter::update_motion(threevector &ChangeInPosition,threevector &ChangeInVelocity){
	M_Debug("\tIn Matter::update_motion(threevector &ChangeInPosition,threevector &ChangeInVelocity)\n\n");
	// Calculate new position
	St.DustPosition += ChangeInPosition;
	St.DustVelocity += ChangeInVelocity;
};

void Matter::update_charge(double potential, double deltat, double deltas){
	M_Debug("\tIn Matter::update_charge(double potential)\n\n");
	St.Potential = potential;
	
	St.DeltaTherm = deltat;
	St.DeltaSec = deltas;

	if ( (St.DeltaSec + St.DeltaTherm) >= 1.0 || St.Potential < 0.0 ){ // If the grain is in fact positive ...
		St.Positive = true;
	}else{
		St.Positive = false;
	}
}

