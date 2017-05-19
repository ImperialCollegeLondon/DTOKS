//#define PAUSE
//#define HEATING_DEBUG

#include "HeatingModel.h"
#include "Constants.h"
#include "Checks.h"
#include "Functions.h"

// Default Constructor, no arguments
HeatingModel::HeatingModel():Type("constant"){
	Defaults();
	Matter *Sample = new Tungsten(10.0);
}

// Basic constructor, only specify power supplied to the dust grain
HeatingModel::HeatingModel(std::string name, double power):Type(name){
	Defaults();
	Matter *Sample = new Tungsten(10.0);
	PowerIncident = power;
}

// Constructor specifying, dust radius, temperature, time step and an array of bools specifying the different heating models.
// In this case, all material constants are assumed to be true constants taking value 'c'.
HeatingModel::HeatingModel(std::string name, char element, double radius, double temp, double power, std::array<bool,9> &models, 
				double timestep ):Type(name){
        Defaults();
	UseModel = models;
	PowerIncident = power; 		// kJ, Power incident
	TimeStep = timestep;
	CheckPos(radius,"Radius");	// Check input values are positive
	CheckPos(timestep,"Time Step");	// Check input values are positive
	if	(element=='w'||element=='W') Sample = new Tungsten(radius,temp);
	else if	(element=='f'||element=='F') Sample = new Iron(radius,temp);
	else if	(element=='g'||element=='G') Sample = new Graphite(radius,temp);
	else if	(element=='b'||element=='B') Sample = new Beryllium(radius,temp);
	else{
		std::cerr << "\nError! Unrecongnised Element Character: -" << element << "-\n";
		assert( strchr("wWfFgGbB",element) );
	}
}

// Constructor specifying in addition to the previoous, the constant models being used.
HeatingModel::HeatingModel(std::string name, char element, double radius, double temp, double power, std::array<bool,9> &models, 
				std::array<char,4> &constmodels, double timestep ):Type(name){
        Defaults();
	UseModel = models;
	PowerIncident = power; 		// kJ, Power incident
	TimeStep = timestep;
	CheckPos(radius,"Radius");
	CheckPos(timestep,"Time Step");
	if	(element=='w'||element=='W') Sample = new Tungsten(radius,temp);
	else if	(element=='f'||element=='F') Sample = new Iron(radius,temp);
	else if	(element=='g'||element=='G') Sample = new Graphite(radius,temp);
	else if	(element=='b'||element=='B') Sample = new Beryllium(radius,temp);
	else{
		std::cerr << "\nError! Unrecongnised Element Character: -" << element << "-\n"; 
		assert( strchr("wWfFgGbB",element) );
	}
}

// Constructor which specifies dust radius, element and ConstModels vector by passing a pointer to a Matter reference.
HeatingModel::HeatingModel(std::string name, char element, double power, std::array<bool,9> &models, double timestep, 
				Matter *&sample, PlasmaData &pdata):Type(name){
	Defaults();
	CheckPos(timestep,"Time Step");
	UseModel 	= models;
	PowerIncident 	= power; 		// kJ, Power incident
	TimeStep 	= timestep;		// s,
	Sample 		= sample;		// Contains element data
	Pdata 		= pdata;		// Contains plasma data
}

void HeatingModel::Defaults(){
	UseModel = {true,false,false,false,false,false,false};	
	PowerIncident = 0.5;
	Pdata.AmbientTemp = 300;
	OldTemp = Sample->get_temperature(); 	// Same sign as temperature 
	TimeStep = 0.1;
	TotalTime = 0;
	ThermalEquilibrium = false;
	ForceNegative = true;

	Pdata.NeutralDensity =  3e19;		// m^-3, Neutral density
	Pdata.ElectronDensity = 8e17;	 	// m^-3, Electron density
	Pdata.Potential = 1;			// arb, assumed negative, potential normalised to dust temperature, (-e*phi)/(Kb*Td)
	Pdata.IonTemp = 100*1.16e5;	 	// K, Ion Temperature
	Pdata.ElectronTemp = 100*1.16e5;	// K, Electron Temperature, convert from eV
	Pdata.NeutralTemp = 100*1.16e5; 	// K, Neutral Temperature, convert from eV
}

void HeatingModel::Reset(std::string filename, double radius, double temp, double timestep, PlasmaData &pdata, Matter *&sample){
	HeatingFile.close();
	ThermalEquilibrium = false;
	OldTemp = 0;
	TimeStep = timestep;
	Sample = sample;
	Pdata = pdata;
	TotalTime = 0;
}

void HeatingModel::CreateFile(std::string filename, bool PrintPhaseData, std::array<char,4> &ConstModels){

	HeatingFile.open(filename);
	HeatingFile << "Time\tTemp\tMass\tDensity\tEnergyIn";
	if( PrintPhaseData ) 					HeatingFile << "\tFusionE\tVapourE";
	if( ConstModels[0] == 'v' && ConstModels[0] == 'V' ) 	HeatingFile << "\tCv";
	if( ConstModels[1] == 'v' && ConstModels[1] == 'V' ) 	HeatingFile << "\tVapourP";
	if( ConstModels[2] == 'v' && ConstModels[2] == 'V' ) 	HeatingFile << "\tLinearExpansion";
       	if( UseModel[0] )       				HeatingFile << "\tEmissLoss\tEmissiv";
       	if( UseModel[1] )       				HeatingFile << "\tEvapRate\tEvapLoss\tEvapMassLoss";
       	if( UseModel[2] )       				HeatingFile << "\tNewton";
       	if( UseModel[3] )       				HeatingFile << "\tIonFlux\tIonHeatFlux";
       	if( UseModel[4] )       				HeatingFile << "\tElectronFlux\tElectronHeatFlux";
       	if( UseModel[5] )       				HeatingFile << "\tNeutronFlux\tNeutronHeatFlux";
       	if( UseModel[6] )       				HeatingFile << "\tNeutralRecomb";
       	if( UseModel[7] )       				HeatingFile << "\tSEE";
       	if( UseModel[8] )       				HeatingFile << "\tTEE";
	HeatingFile << "\n";
}

const int HeatingModel::Vapourise(std::string filename, std::array<char,4> &ConstModels, char TimeStepType){
	H_Debug("\n\nIn HeatingModel::Vapourise(char EmissivModel, char TimeStepType)");

	CreateFile(filename,false, ConstModels);
	
	Sample->update_models(ConstModels);			// Update density, dimensions, Cp and emissivity
	while( !Sample->is_gas() && !ThermalEquilibrium ){ // If the sample is gaseous or in TE, the model ends.
		Heat(ConstModels,TimeStepType);
	}

	int rValue(0);
	if( Sample->is_gas() && Sample->get_superboilingtemp() <= Sample->get_temperature() ){
		std::cout << "\n\nSample has Boiled ";
		rValue = 1;
	}else if( Sample->is_gas() && Sample->get_superboilingtemp() > Sample->get_temperature() ){
		std::cout << "\n\nSample has Evaporated ";
		rValue = 2;
	}else if( ThermalEquilibrium ){
		std::cout << "\n\nSample has reached Thermal Equilibrium ";
		rValue = 3;
	}
	std::cout << "at T = " << Sample->get_temperature() << "K in " << TotalTime << "s!\n\n*********\n\n";
	HeatingFile.close();
	return rValue; // 0, running normally. 1; Sample boiled. 2; Sample Evaporated. 3; Thermal equilibrium.
}

void HeatingModel::Heat(std::array<char,4> &ConstModels, char TimeStepType){
	H_Debug("\n\nIn HeatingModel::Heat(char EmissivModel, char TimeStepType)");

	assert( Sample->get_mass() > 0 );

	double TotalPower = CalculatePower(Sample->get_temperature());	// Take Eularian step to get initial time step
	CheckTimeStep(TotalPower,TimeStepType);				// Check Time step length is appropriate

//	std::cout << "\nTotalPower = " << TotalPower; 
	// This is a complete Hack...
//	if( ThermalEquilibrium == true && UseModel[1] == true){ // If it's evaporating and we want to record evaporation
//		ThermalEquilibrium = false;
//		TimeStepType = 'c'; // Set time step such that it produces a change in mass of 0.001%
//		TimeStep = (0.00001*Sample->get_mass()*AvNo)/(EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass());
//	}
	double TotalEnergy = RungeKutta4();			// Calculate total energy through RungeKutta4 method

	// Account for evaporative mass loss
	if( UseModel[1] && Sample->is_liquid() )
		Sample->update_mass( (TimeStep*EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass())/AvNo );

	Sample->update();					// Update the values of constants
  	Sample->update_temperature(TotalEnergy);		// Update Temperature
  	Print(TotalEnergy/TimeStep,ConstModels);		// Print data to file
//	Print(TotalPower,ConstModels);

	// This is also a part of the Hack...
//	if( ThermalEquilibrium == true && UseModel[1] == true && Sample->get_temperature() > Sample->get_meltingtemp() ){ // If it's evaporating and we want to record evaporation
//		ThermalEquilibrium = false;
//		std::cout << "\n--- Thermal Equilibrium = " << ThermalEquilibrium << " ---\n";
//	}
}

double HeatingModel::CalculatePower(double DustTemperature)const{
	H_Debug( "\n\nIn HeatingModel::CalculatePower(double DustTemperature)" );
	double TotalPower = PowerIncident;

	if( UseModel[0] )				TotalPower -= EmissivityModel		(DustTemperature);
	if( UseModel[1] && Sample->is_liquid() )	TotalPower -= EvaporationModel		(DustTemperature); 		
	if( UseModel[2] )				TotalPower -= NewtonCooling		(DustTemperature);
	if( UseModel[3] )				TotalPower += IonHeatFlux		(DustTemperature);
	if( UseModel[4] )				TotalPower += ElectronHeatFlux		(DustTemperature);
	if( UseModel[5] )				TotalPower += NeutralHeatFlux		();
	if( UseModel[6] )				TotalPower += NeutralRecombination	(DustTemperature);
	if( UseModel[7] )				TotalPower -= SEE			(DustTemperature);
	if( UseModel[8] )				TotalPower -= TEE			(DustTemperature);	

	H_Debug("\n\nTotalPower = \t" 		<< TotalPower 				<< "kW");
	H_Debug("\nEmissivityModel() = \t" 	<< -EmissivityModel(DustTemperature) 	<< "kW");
	H_Debug("\nEvaporationModel() = \t" 	<< -EvaporationModel(DustTemperature) 	<< "kW");
	H_Debug("\nNewtonCooling() = \t" 	<< -NewtonCooling(DustTemperature) 	<< "kW");
	H_Debug("\nIonHeatFlux() = \t" 		<< IonHeatFlux(DustTemperature) 	<< "kW");
	H_Debug("\nElectronHeatFlux() = \t" 	<< ElectronHeatFlux(DustTemperature) 	<< "kW");
	H_Debug("\nNeutralHeatFlux() = \t" 	<< NeutralHeatFlux() 			<< "kW");
	H_Debug("\nNeutralRecombination() = \t"	<< NeutralRecombination(DustTemperature)<< "kW");
	H_Debug("\nSEE() = \t" 			<< -SEE(DustTemperature) 		<< "kW");
	H_Debug("\nTEE() = \t" 			<< -TEE(DustTemperature) 		<< "kW\n");
	return TotalPower;
}

double HeatingModel::RungeKutta4(){
	double k1 = CalculatePower(Sample->get_temperature()); 
	if(k1<0 && fabs(k1/2) > Sample->get_temperature()){
		ThermalEquilibrium = true;
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k1 step negative and larger than Td!";
		return 0;
	}
	double k2 = CalculatePower(Sample->get_temperature()+k1/2); 
	if( k2<0 && fabs(k2/2) > Sample->get_temperature() ){
		ThermalEquilibrium = true;
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k2 step negative and larger than Td!";
		return (TimeStep/6)*k1;
	}
	double k3 = CalculatePower(Sample->get_temperature()+k2/2);
	if( k3<0 && fabs(k3) > Sample->get_temperature() ){
		ThermalEquilibrium = true;
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k3 step negative and larger than Td!";
		return (TimeStep/6)*(k1+2*k2);
	}
	double k4 = CalculatePower(Sample->get_temperature()+k3);
	if( k4<0 && fabs(k4/2) > Sample->get_temperature() ){
		ThermalEquilibrium = true;
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k4 step negative and larger than Td!";
		return (TimeStep/6)*(k1+2*k2+2*k3);
	}
	return (TimeStep/6)*(k1+2*k2+2*k3+k4);
};

void HeatingModel::CheckTimeStep(double TotalPower, char TimeStepType){
	H_Debug( "\n\nIn HeatingModel::CheckTimeStep(double TotalPower, char TimeStepType)" );
	// Deal with case where power/time step causes large temperature change.
	assert(TimeStep > 0);
	double DeltaTempTest = TotalPower*TimeStep/(Sample->get_mass()*Sample->get_heatcapacity());
	H_Debug( "\nTotalPower = " << TotalPower << "\nTimeStep " << TimeStep << "\nSample->get_mass() = " << Sample->get_mass()
			<< "\nSample->get_heatcapacity() = " << Sample->get_heatcapacity());

	// This model allows custom time steps but will set the time step such that it is AT LEAST
	// small enough to only change the temperature by 1 degree.
	if(TimeStepType == 's' || TimeStepType == 'S'){
		if( fabs(DeltaTempTest) > 1 ){
			std::cout << "\n\nWarning! Large temperature change = " << DeltaTempTest << " K";
			std::cout << "\nChanging Time step\nFrom dt = " << TimeStep << " s ..."; 
			while( fabs(DeltaTempTest) > 0.1 && DeltaTempTest!=0){
				TimeStep = TimeStep/10;
				DeltaTempTest = TotalPower*TimeStep/(Sample->get_mass()*Sample->get_heatcapacity());
			}
			std::cout << "\nto dt = " << TimeStep << " s";
		}
	}else if(TimeStepType == 'o' || TimeStepType == 'O'){
		// This model forces the time step to be the value which produces a change in temperature or 1 degree
		TimeStep = fabs((Sample->get_mass()*Sample->get_heatcapacity())/TotalPower);
	}else if(TimeStepType == 't' || TimeStepType == 'T'){
		// This model forces the time step to be the value which produces a change in temperature or 0.01 degree
		TimeStep = fabs((Sample->get_mass()*Sample->get_heatcapacity())/(TotalPower*100));
	}else if(TimeStepType == 'f' || TimeStepType == 'F'){
		if(DeltaTempTest>2){
			static bool runOnce = true;
			WarnOnce(runOnce,"For fixed time step, DeltaTempTest > 2!"); // Check temperature change isn't too large
			std::cout << "\n\nTemp = " << Sample->get_temperature();
			std::cout << "\nDeltaTemp = " << DeltaTempTest;
		}
		
	}
	assert(TimeStep > 0 && TotalTime/TimeStep < 1e10);

	if( TotalTime/TimeStep > 1e9 ){
		static bool runOnce = true;
		std::cout << "\n\nThermal Equilibrium reached on condition (4):";
		WarnOnce(runOnce,"At least 10^8 steps taken, Thermal Equilibrium assumed.");
		ThermalEquilibrium = true;
	}

	if( (Sample->get_temperature()-OldTemp > 0 && DeltaTempTest < 0) // If temperature change sign changed between last step
		|| (Sample->get_temperature()-OldTemp < 0 && DeltaTempTest > 0) ){
		if( Sample->get_temperature() != Sample->get_superboilingtemp() ){
			std::cout << "\n\nThermal Equilibrium reached on condition (1): Sign change of Temperature change!";
			ThermalEquilibrium = true;
			TimeStep = 0;
		}
	}if(  fabs(DeltaTempTest/TimeStep) < 0.01 ){ // If Temperature gradient is less than 1%
		if( Sample->get_temperature() != Sample->get_superboilingtemp() ){
			std::cout << "\n\nThermal Equilibrium reached on condition (2): Temperature Gradient < 0.01!";
			ThermalEquilibrium = true;
			TimeStep = 0;
		}
	}

	OldTemp = Sample->get_temperature();

	H_Debug( "\nTimeStep = " << TimeStep );
	TotalTime += TimeStep;
}

void HeatingModel::Print(double TotPower, std::array<char,4> &ConstModels){
	H_Debug("\n\nIn HeatingModel::Print(double TotPower, std::array<char,4> &ConstModels)");

	HeatingFile 	<< TotalTime << "\t" << Sample->get_temperature() << "\t" << Sample->get_mass() << "\t" << Sample->get_density()
		<< "\t" << TotPower;

	bool PrintPhaseData = false;
	if( PrintPhaseData )	HeatingFile 	<< "\t" << Sample->get_fusionenergy() << "\t" << Sample->get_vapourenergy();

	if( ConstModels[0] == 'v' && ConstModels[0] == 'V' ) HeatingFile << "\t" << Sample->get_heatcapacity();
	if( ConstModels[1] == 'v' && ConstModels[1] == 'V' ) HeatingFile << "\t" << Sample->get_vapourpressure();
	if( ConstModels[2] == 'v' && ConstModels[2] == 'V' ) HeatingFile << "\t" << Sample->get_linearexpansion();


	if( UseModel[0] ) 	HeatingFile 	<< "\t" << EmissivityModel(Sample->get_temperature()) << "\t" << Sample->get_emissivity();
	if( UseModel[1] && Sample->is_liquid() )	
				HeatingFile 	<< "\t" << EvaporationFlux(Sample->get_temperature()) 
					<< "\t" << EvaporationModel(Sample->get_temperature()) 
					<< "\t" << EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass()/AvNo;
	else if( UseModel[1] )	HeatingFile 	<< "\t" << 0 << "\t" << 0 << "\t" << 0;
	if( UseModel[2] )	HeatingFile 	<< "\t" << NewtonCooling(Sample->get_temperature());
	if( UseModel[3] )	HeatingFile 	<< "\t" << IonFlux(Sample->get_temperature()) 
					<< "\t" << IonHeatFlux(Sample->get_temperature());
	if( UseModel[4] )	HeatingFile 	<< "\t" << ElectronFlux(Sample->get_temperature()) << "\t" 
					<< ElectronHeatFlux(Sample->get_temperature());
	if( UseModel[5] )	HeatingFile 	<< "\t" << NeutralFlux() << "\t" << NeutralHeatFlux();
	if( UseModel[6] )	HeatingFile 	<< "\t" << NeutralRecombination(Sample->get_temperature());
	if( UseModel[7] )	HeatingFile 	<< "\t" << SEE(Sample->get_temperature());
	if( UseModel[8] )	HeatingFile 	<< "\t" << TEE(Sample->get_temperature());
	HeatingFile << "\n";
}


// *************************************************** HEATING MODELS *************************************************** //

// http://users.wfu.edu/ucerkb/Nan242/L06-Vacuum_Evaporation.pdf, 
// https://en.wikipedia.org/wiki/Hertz%E2%80%93Knudsen_equation
// Using Hertz–Knudsen equation, returns Energy lost per second in Kila Joules
// MASS LOSS EQUATION 
// http://www.leb.eei.uni-erlangen.de/winterakademie/2006/result/content/course01/pdf/0102.pdf
const double HeatingModel::EvaporationModel(double DustTemperature)const{
	H_Debug("\n\nIn HeatingModel::EvaporationModel():");

	// Approximate emitted energy as mean maxwell
	// This used to be multiplied by 1000 which is why it was significant.
	double MaxwellEnergy = (3*Kb*DustTemperature/(2*1000)); // Converted to kJ.
// 	double MaxwellMeanVelocity = sqrt((8*Kb*DustTemperature)/(AtomicMass*1000*AMU));

	double EvapFlux = EvaporationFlux(DustTemperature);

	if( EvapFlux != EvapFlux ){ 
		std::cout << "\n\nError! EvapFlux = " << EvapFlux << "\n";
		throw std::exception(); 
	}
	H_Debug("\n\nIn HeatingModel::EvaporationModel()\nMaxwellEnergy = " << MaxwellEnergy << " kJ" << "\nEvapFlux = " 
		<< EvapFlux << " s^-1\nBondEnergy = " << Sample->get_bondenergy()/AvNo << " kJ\nreturn = " 
		<< EvapFlux*(MaxwellEnergy+Sample->get_bondenergy()/AvNo) << "\nDustTemperature = " << DustTemperature 
		<<  "\nSample->get_surfacearea() = " << Sample->get_surfacearea()
		<< "\nSample->get_vapourpressure() = " << Sample->get_vapourpressure());
	
	// See ElementData.h for more info on 'bondenergy'. Added to account for energy lost by breaking bonds.
	return EvapFlux*(MaxwellEnergy+Sample->get_bondenergy()/AvNo); 
}

const double HeatingModel::EvaporationFlux(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::EvaporationFlux():");
	double AmbientPressure = 0;
	double StickCoeff = 1.0;
	return (StickCoeff*Sample->get_surfacearea()*AvNo*(Sample->get_vapourpressure()-AmbientPressure))/
			sqrt(2*PI*Sample->get_atomicmass()*R*DustTemperature);
}

// Using Stefan-Boltzmann Law, returns Energy lost per second in Kila Joules
const double HeatingModel::EmissivityModel(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::EmissivityModel():");
//	H_Debug("\nSample->get_surfacearea() = " << Sample->get_surfacearea() << "\nSigma = " << Sigma 
//		<< "\nDustTemperature = " << DustTemperature);
	// Energy emitted from a sample converted to kJ
	return Sample->get_emissivity()*Sample->get_surfacearea()*Sigma*(pow(DustTemperature,4)-pow(Pdata.AmbientTemp,4))/1000;
}

// VERY APPROXIMATE MODEL: Atmosphere assumed to be 300 degrees always, rough heat transfer coefficient is use
// https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling
const double HeatingModel::NewtonCooling(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::NewtonCooling():");
	static bool runOnce = true;
	WarnOnce(runOnce,"In HeatingModel::NewtonCooling():\nHeatTransair Coefficient wrong for Tungsten, Beryllium and Graphite.");

	return (Sample->get_heattransair()*Sample->get_surfacearea()*(DustTemperature-Pdata.AmbientTemp))/1000; // convert to kJ
}	

// Neutral Recombination assuming Rn=0; fraction of backscattered ions/neutrals is zero.
const double HeatingModel::NeutralRecombination(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::NeutralRecombination():");
//	H_Debug("\n14.7*echarge - 2*Kb*DustTemperature*NeutralFlux()\n");
//	H_Debug("\nIonFlux = " << IonFlux());
//	H_Debug("\nReturn = " << (14.7*echarge*IonFlux()));
	double RN(0), RE(0);
//	backscatter(Pdata.ElectronTemp,Pdata.IonTemp,Mi,Pdata.Potential,Element,RE,RN);
	H_Debug( "\nRN = " << RN );
	if( RN > 0.1 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"In HeatingModel::NeutralRecombination(double DustTemperature)\nRN > 0.1. Neutral Recombination affected by backscattering by more than 10%!");
	}

	return Sample->get_surfacearea()*(1-RN)*(14.7*echarge*IonFlux(DustTemperature) - 2*Kb*DustTemperature*NeutralFlux())/1000; // Convert from J to kJ
	return Sample->get_surfacearea()*(14.7*echarge*IonFlux(DustTemperature) - 2*Kb*DustTemperature*NeutralFlux())/1000; // Convert from J to kJ
}

const double HeatingModel::SEE(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::SEE():");
//	H_Debug("\neFlux*DeltaSec()*(3*echarge+Sample->get_bondenergy())\n");
//	H_Debug("\neFlux=" << eFlux << "\nDeltaSec() = " << DeltaSec() 
//			<< "\nSample->get_bondenergy() = " << Sample->get_bondenergy());
//	H_Debug("\nReturn = " << eFlux*DeltaSec()*(3*echarge+Sample->get_bondenergy()));
	double SEE=0;
	if(ForceNegative || DeltaTot(DustTemperature) <= 1 )
		SEE = Sample->get_surfacearea()*ElectronFlux(DustTemperature)*DeltaSec()*echarge*
                                        (3+Sample->get_workfunction())/1000; // Convert
	else if( DeltaTot(DustTemperature) > 1 ) 	SEE = 0; // Electrons captured by positive grain
	return SEE;
}

const double HeatingModel::TEE(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::TEE():");
//	H_Debug("\neFlux*DeltaTherm()*(2*Kb*DustTemperature+Sample->get_bondenergy())\n");
//	H_Debug("\neFlux=" << eFlux << "\nDeltaTherm() = " << DeltaTherm() << "\n2*Kb*DustTemperature = "
//			<< 2*Kb*DustTemperature << "\nSample->get_bondenergy() = " << Sample->get_bondenergy());
//	H_Debug("\nReturn = " << eFlux*DeltaTherm()*(2*Kb*DustTemperature+Sample->get_bondenergy()));
	double TEE=0;

	if( ForceNegative || DeltaTot(DustTemperature) <= 1 )
		TEE = Sample->get_surfacearea()*DeltaTherm(DustTemperature)*
			(2*Kb*DustTemperature+echarge*Sample->get_workfunction())/1000; // Convert to kJ
	else if( DeltaTot(DustTemperature) > 1 ) 	TEE = 0; // Electrons captured by positive grain
	return TEE;
}

const double HeatingModel::IonHeatFlux(double DustTemperature)const{ // Assuming Re = 0
	H_Debug("\n\nIn HeatingModel::IonHeatFlux():");
	H_Debug("\nIonFlux() = " << IonFlux(DustTemperature) << "\nPdata.Potential = " << Pdata.Potential);

	double RN(0), RE(0);
//	backscatter(Pdata.ElectronTemp,Pdata.IonTemp,Mi,Pdata.Potential,Element,RE,RN);
	H_Debug( "\nRE = " << RE );
	if( RE > 0.1 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"In HeatingModel::IonHeatFlux(double DustTemperature)\nRE > 0.1. Ion Heat Flux affected by backscattering by more than 10%!");
	}
	return (Sample->get_surfacearea()*(1-RE)*IonFlux(DustTemperature)*Pdata.IonTemp*Kb/1000) // Convert from Joules to KJ
	*(2+2*Pdata.Potential*(Pdata.ElectronTemp/Pdata.IonTemp)+pow(Pdata.Potential*(Pdata.ElectronTemp/Pdata.IonTemp),2))/(1+Pdata.Potential*(Pdata.ElectronTemp/Pdata.IonTemp)); 

}

const double HeatingModel::ElectronHeatFlux(double DustTemperature)const{ // Only for a negative grain
//	H_Debug("\n\nIn HeatingModel::ElectronHeatFlux():");
//	H_Debug("\nSample->get_surfacearea()*2*ElectronFlux()*Pdata.ElectronTemp*Kb/1000\n");
	return Sample->get_surfacearea()*2*ElectronFlux(DustTemperature)*Pdata.ElectronTemp*Kb/1000; // Convert from Joules to KJ
}

const double HeatingModel::NeutralHeatFlux()const{
//	H_Debug("\n\nIn HeatingModel::NeutralHeatFlux():");
//	H_Debug("\nNeutralFlux()*Pdata.NeutralTemp*Kb*2\n");
	return Sample->get_surfacearea()*2*NeutralFlux()*Pdata.NeutralTemp*Kb/1000; // Convert from Joules to KJ
}

const double HeatingModel::IonFlux(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::IonFlux():");
//	H_Debug("\nElectronFlux()*(1-DeltaTherm()-DeltaSec())\n");
//	H_Debug("\nElectronFlux() =" << ElectronFlux() << "\nDeltaTherm() = " << DeltaTherm() 
//			<< "\nDeltaSec() = " << DeltaSec());
//	H_Debug("\nReturn = " << ElectronFlux()*(1-DeltaTherm()-DeltaSec()));
	double IonFlux=0;

	if(ForceNegative || DeltaTot(DustTemperature) > 1 )IonFlux = ElectronFlux(DustTemperature); //Positive grain, DeltaTot() > 1
	else	IonFlux = ElectronFlux(DustTemperature)*(1-DeltaTherm(DustTemperature)-DeltaSec());
	return IonFlux;
}

const double HeatingModel::ElectronFlux(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::ElectronFlux():");
//	H_Debug("\nPdata.NeutralDensity*exp(-Pdata.Potential)*sqrt(Kb*DustTemperature/(2*PI*Me))\n");
//	H_Debug("\nexp(-Pdata.Potential) = " << exp(-Pdata.Potential)
//			<< "\n(2*PI*Me) = " << (2*PI*Me) << "sqrt(Kb*DustTemperature/(2*PI*Me)) = " 
//			<< sqrt(Kb*DustTemperature/(2*PI*Me)));
//	H_Debug("\nReturn = " << Pdata.NeutralDensity*exp(-Pdata.Potential)*sqrt(Kb*DustTemperature/(2*PI*Me)));
//	std::cout << "\neFlux = " << Pdata.ElectronDensity*exp(-Pdata.Potential)*sqrt(Kb*Pdata.ElectronTemp/(2*PI*Me)); std::cin.get();
	return Pdata.ElectronDensity*exp(-Pdata.Potential)*sqrt(Kb*Pdata.ElectronTemp/(2*PI*Me));
}

const double HeatingModel::NeutralFlux()const{
//	H_Debug("\n\nIn HeatingModel::NeutralFlux():");
//	H_Debug("\nPdata.NeutralTemp*sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mi))/4\n");
//	H_Debug("(2*PI*Mi) = " << (2*PI*Mi) << "\nsqrt(Kb*Pdata.NeutralTemp/(2*PI*Mi)) = "<< sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mi)));
//	H_Debug("\nReturn = " << 2*Pdata.NeutralTemp*Pdata.NeutralDensity*sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mi)));
	return Pdata.NeutralDensity*sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mi));
}
const double HeatingModel::DeltaTherm(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::DeltaTherm():");
//	H_Debug("\n(Kb*DustTemperature+Sample->get_bondenergy())*Richardson*pow(DustTemperature,2)*exp(-(Sample->get_bondenergy())/(Kb*DustTemperature))\n");
//	H_Debug("\nDeltaTherm = " << (Richardson*pow(DustTemperature,2)*exp(-(Sample->get_workfunction()*echarge)/(Kb*DustTemperature)))/(echarge);
	return (Richardson*pow(DustTemperature,2)*exp(-(Sample->get_workfunction()*echarge)/(Kb*DustTemperature)))
		/echarge;
}

const double HeatingModel::DeltaSec()const{
//	H_Debug("\n\nIn HeatingModel::DeltaSec():");
//	H_Debug("\nPdata.ElectronTemp = " << Pdata.ElectronTemp);
//	H_Debug("\nDeltaSec = " << sec(10,Element));
	return sec(Pdata.ElectronTemp/1.16e5,Sample->get_elem()); 	// Convert from K to ev
}

const double HeatingModel::DeltaTot(double DustTemperature)const{
	return (DeltaSec() + DeltaTherm(DustTemperature));
}

// ************************************* \\
