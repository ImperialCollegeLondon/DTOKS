//#define PAUSE
//#define HEATING_DEBUG
//#define HEATING_DEEP_DEBUG

#include "HeatingModel.h"
#include "Constants.h"
#include "Checks.h"
#include "Functions.h"

HeatingModel::HeatingModel():Type("constant"),Model(){
	H_Debug("\n\nIn HeatingModel::HeatingModel():Type(constant),Model()\n\n");
	Defaults();
	CreateFile("Default_Heating_filename.txt",false);
}


// Constructor which specifies dust radius, by passing a pointer to a Matter reference.
HeatingModel::HeatingModel(std::string filename, double accuracy, std::array<bool,9> &models,
				Matter *& sample, PlasmaData &pdata) : Model(sample,pdata,accuracy){
	H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, double accuracy, std::array<bool,9> &models, Matter *& sample, PlasmaData const &pdata) : Model(sample,pdata,accuracy)\n\n");
	Defaults();
	UseModel 		= models;
	CreateFile(filename,false);
}

HeatingModel::HeatingModel(std::string filename, double accuracy, std::array<bool,9> &models,
				Matter *& sample, PlasmaGrid &pgrid) : Model(sample,pgrid,accuracy){
	H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename,double accuracy, std::array<bool,9> &models, Matter *& sample, PlasmaGrid const &pgrid) : Model(sample,pgrid,accuracy)\n\n");
	Defaults();
	UseModel 	= models;
	CreateFile(filename,false);
}

void HeatingModel::Defaults(){
	H_Debug("\tIn HeatingModel::Defaults()\n\n");
	UseModel = {false,false,false,false,false,false,false};
	PowerIncident = 0;			// kW, Power Incident
	TimeStep = 0;				// s, Time step
	TotalTime = 0;				// s, Total Time
	ForceNegative = true;			// Force the sample to behave as negatively charged
}

void HeatingModel::CreateFile(std::string filename, bool PrintPhaseData){
	H_Debug("\tIn HeatingModel::CreateFile(std::string filename, bool PrintPhaseData)\n\n");
	ModelDataFile.open(filename);
	ModelDataFile << "Time\tTemp\tMass\tDensity";
	if( PrintPhaseData ) 						ModelDataFile << "\tFusionE\tVapourE";
	if( Sample->get_c(0) == 'v' || Sample->get_c(0) == 'V' ) 	ModelDataFile << "\tCv";
	if( Sample->get_c(1) == 'v' || Sample->get_c(1) == 'V' ) 	ModelDataFile << "\tVapourP";
	if( Sample->get_c(2) == 'v' || Sample->get_c(2) == 'V' ) 	ModelDataFile << "\tLinearExpansion";
       	if( UseModel[0] )       					ModelDataFile << "\tEmissLoss";
       	if( UseModel[0] && Sample->get_c(3) == 'v' )   			ModelDataFile << "\tEmissiv";
       	if( UseModel[1] )       					ModelDataFile << "\tEvapRate\tEvapLoss\tEvapMassLoss";
       	if( UseModel[2] )       					ModelDataFile << "\tNewton";
       	if( UseModel[3] )       					ModelDataFile << "\tIonFlux\tIonHeatFlux";
       	if( UseModel[4] )       					ModelDataFile << "\tElectronFlux\tElectronHeatFlux";
       	if( UseModel[5] )       					ModelDataFile << "\tNeutronFlux\tNeutronHeatFlux";
       	if( UseModel[6] )       					ModelDataFile << "\tNeutralRecomb";
       	if( UseModel[7] )       					ModelDataFile << "\tSEE";
       	if( UseModel[8] )       					ModelDataFile << "\tTEE";
	ModelDataFile << "\n";
}

const int HeatingModel::Vapourise(){
	H_Debug("\tIn HeatingModel::Vapourise()\n\n");
	// If the sample is gaseous or in TE (Given that the plasma is continuous), the model ends.
	while( !Sample->is_gas() ){ 
		Heat();
		ThermalEquilibrium = false;
		if( ContinuousPlasma && ThermalEquilibrium )
			break;
	}

	int rValue(0);
	if( Sample->is_gas() && Sample->get_superboilingtemp() <= Sample->get_temperature() ){
		std::cout << "\n\nSample has Boiled ";
		rValue = 1;
	}else if( Sample->is_gas() && Sample->get_superboilingtemp() > Sample->get_temperature() ){
		std::cout << "\n\nSample has Evaporated ";
		rValue = 2;
	}else if( ThermalEquilibrium && ContinuousPlasma ){
		std::cout << "\n\nSample has reached Thermal Equilibrium in Continuous Plasma.";
		rValue = 3;
	}
	std::cout << "\nat T = " << Sample->get_temperature() << "K in " << TotalTime << "s!\n\n*********\n\n";
	ModelDataFile.close();
	return rValue; // 0, running normally. 1; Sample boiled. 2; Sample Evaporated. 3; Thermal equilibrium.
}

void HeatingModel::Heat(double timestep){
	H_Debug("\tIn HeatingModel::Heat(double timestep)\n\n");
	
	// Make sure timestep input time is valid. Shouldn't exceed the timescale of the process.
	assert(timestep > 0 && timestep <= TimeStep );
	TimeStep = timestep;
	std::cout << "\nTimeStep = " << TimeStep << "\n";
	assert( Sample->get_mass() > 0 );

	double TotalEnergy = RungeKutta4();                     // Calculate total energy through RungeKutta4 method
	H1_Debug( "\tTotalEnergy = " << TotalEnergy << "\n");
        Sample->update_temperature(TotalEnergy);                // Update Temperature

	// Account for evaporative mass loss
	if( UseModel[1] && Sample->is_liquid() )
		Sample->update_mass( (TimeStep*EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass())/AvNo );
	
	Sample->update();
        H_Debug("\t"); Print();                // Print data to file
	std::cout << "\nTemperature = " << Sample->get_temperature();
	TotalTime += TimeStep;
}

void HeatingModel::Heat(){
	H_Debug("\tIn HeatingModel::Heat()\n\n");
	

	assert( Sample->get_mass() > 0 );

	double TotalEnergy = RungeKutta4();                     // Calculate total energy through RungeKutta4 method
	H1_Debug( "\tTotalEnergy = " << TotalEnergy << "\n");
        Sample->update_temperature(TotalEnergy);                // Update Temperature

	// Account for evaporative mass loss
	if( UseModel[1] && Sample->is_liquid() )
		Sample->update_mass( (TimeStep*EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass())/AvNo );
	
	Sample->update();
        H_Debug("\t"); Print();                // Print data to file
	std::cout << "\nTemperature = " << Sample->get_temperature();
	TotalTime += TimeStep;
}

double HeatingModel::CalculatePower(double DustTemperature)const{
	H_Debug( "\tIn HeatingModel::CalculatePower(double DustTemperature)\n\n");
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

	H1_Debug("\n\nPowerIncident = \t" 	<< PowerIncident			<< "kW");
	if( UseModel[0] )	H1_Debug("\nEmissivityModel() = \t" 	<< -EmissivityModel(DustTemperature) 	<< "kW");
	if( UseModel[1] && Sample->is_liquid() )	
				H1_Debug("\nEvaporationModel() = \t" 	<< -EvaporationModel(DustTemperature) 	<< "kW");
	if( UseModel[2] )	H1_Debug("\nNewtonCooling() = \t" 	<< -NewtonCooling(DustTemperature) 	<< "kW");
	if( UseModel[3] )	H1_Debug("\nIonHeatFlux() = \t" 	<< IonHeatFlux(DustTemperature) 	<< "kW");
	if( UseModel[4] )	H1_Debug("\nElectronHeatFlux() = \t" 	<< ElectronHeatFlux(DustTemperature) 	<< "kW");
	if( UseModel[5] )	H1_Debug("\nNeutralHeatFlux() = \t" 	<< NeutralHeatFlux() 			<< "kW");
	if( UseModel[6] )	H1_Debug("\nNeutralRecombination() = \t"<< NeutralRecombination(DustTemperature)<< "kW");
	if( UseModel[7] )	H1_Debug("\nSEE() = \t" 		<< -SEE(DustTemperature) 		<< "kW");
	if( UseModel[8] )	H1_Debug("\nTEE() = \t" 		<< -TEE(DustTemperature) 		<< "kW\n");


	return TotalPower;
}

double HeatingModel::RungeKutta4(){
	H_Debug( "\tIn HeatingModel::RungeKutta4()\n\n");
	double k1 = CalculatePower(Sample->get_temperature()); 
	if(k1<0 && fabs(k1/2) > Sample->get_temperature()){
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k1 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return 0;
	}
	double k2 = CalculatePower(Sample->get_temperature()+k1/2); 
	if( k2<0 && fabs(k2/2) > Sample->get_temperature() ){
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k2 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return (TimeStep/6)*k1;
	}
	double k3 = CalculatePower(Sample->get_temperature()+k2/2);
	if( k3<0 && fabs(k3) > Sample->get_temperature() ){
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k3 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return (TimeStep/6)*(k1+2*k2);
	}
	double k4 = CalculatePower(Sample->get_temperature()+k3);
	if( k4<0 && fabs(k4/2) > Sample->get_temperature() ){
		std::cout << "\n\nThermal Equilibrium reached on condition (3): k4 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return (TimeStep/6)*(k1+2*k2+2*k3);
	}
	H1_Debug( "\nTimeStep = " << TimeStep << "\nk1 = " << k1 << "\nk2 =" << k2 << "\nk3 = " << k3 << "\nk4 = " << k4);
	return (TimeStep/6)*(k1+2*k2+2*k3+k4);
};


// This model forces the time step to be the value which produces a change in temperature or 1*accuracy degree
double HeatingModel::UpdateTimeStep(){
	H_Debug( "\tIn HeatingModel::UpdateTimeStep()\n\n" );

	// Take Eularian step to get initial time step
	H_Debug("\t"); double TotalPower = CalculatePower(Sample->get_temperature());

	// COULD BE: If first time step or change in temp is greater than 1 degree, set time step to be equal to 1 degree step
	// Note, if the time step changes conditionally, then the time scale of the process may vary independantly of the time step.
	// This causes issues when deciding on the time ordering of processes.
//	if( TimeStep == 0 || fabs(TimeStep*TotalPower/(Sample->get_mass()*Sample->get_heatcapacity())) > 1 ) 
	TimeStep = fabs((Sample->get_mass()*Sample->get_heatcapacity())/(TotalPower*Accuracy));

	// Check Thermal Equilibrium hasn't been reached
	double DeltaTempTest = TotalPower*TimeStep/(Sample->get_mass()*Sample->get_heatcapacity());

	// May be that we could remove this condition now...
	if( TotalPower == 0 ){	
		static bool runOnce = true;
		WarnOnce(runOnce,"\nWarning! TotalPower = 0\nThermalEquilibrium Assumed, TimeStep being set to unity.");
		ThermalEquilibrium = true;
		TimeStep = 1;
	}

	if( Sample->get_temperature() != Sample->get_superboilingtemp() ){
		if( (Sample->get_temperature()-OldTemp > 0 && DeltaTempTest < 0) // If temperature changed sign changed this step
			|| (Sample->get_temperature()-OldTemp < 0 && DeltaTempTest > 0) ){
			std::cout << "\n\nThermal Equilibrium reached on condition (1): Sign change of Temperature change!";
			ThermalEquilibrium = true;
			TimeStep = 1; 
		}if(  fabs(DeltaTempTest/TimeStep) < 0.01 ){ // If Temperature gradient is less than 1%
			std::cout << "\n\nThermal Equilibrium reached on condition (2): Temperature Gradient < 0.01!";
			ThermalEquilibrium = true;
			TimeStep = 1;
		}
	}

	OldTemp = Sample->get_temperature();

	H1_Debug("\nSample->get_mass() = " << Sample->get_mass() << "\nSample->get_heatcapacity() = " << 
		Sample->get_heatcapacity() << "\nTotalPower = " << TotalPower << "\nAccuracy = " << Accuracy);
	assert(TimeStep > 0 && TimeStep != INFINITY && TimeStep == TimeStep);

	return TimeStep;
}

void HeatingModel::Print(){
	H_Debug("\tIn HeatingModel::Print()\n\n");

	ModelDataFile 	<< TotalTime << "\t" << Sample->get_temperature() << "\t" << Sample->get_mass() 
		<< "\t" << Sample->get_density();

	bool PrintPhaseData = false; // Lol
	if( PrintPhaseData )	ModelDataFile 	<< "\t" << Sample->get_fusionenergy() << "\t" << Sample->get_vapourenergy();

	if( Sample->get_c(0) == 'v' || Sample->get_c(0) == 'V' ) 	ModelDataFile << "\tCv";
	if( Sample->get_c(1) == 'v' || Sample->get_c(1) == 'V' ) 	ModelDataFile << "\tVapourP";
	if( Sample->get_c(2) == 'v' || Sample->get_c(2) == 'V' ) 	ModelDataFile << "\tLinearExpansion";

	if( UseModel[0] ) 	ModelDataFile 	<< "\t" << EmissivityModel(Sample->get_temperature());
	if( UseModel[0] && Sample->get_c(3) ) 	
				ModelDataFile 	<< "\t" << Sample->get_emissivity();
	if( UseModel[1] && Sample->is_liquid() )	
				ModelDataFile 	<< "\t" << EvaporationFlux(Sample->get_temperature()) 
					<< "\t" << EvaporationModel(Sample->get_temperature()) 
					<< "\t" << EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass()/AvNo;
	else if( UseModel[1] )	ModelDataFile 	<< "\t" << 0 << "\t" << 0 << "\t" << 0;
	if( UseModel[2] )	ModelDataFile 	<< "\t" << NewtonCooling(Sample->get_temperature());
	if( UseModel[3] )	ModelDataFile 	<< "\t" << IonFlux(Sample->get_temperature()) 
					<< "\t" << IonHeatFlux(Sample->get_temperature());
	if( UseModel[4] )	ModelDataFile 	<< "\t" << ElectronFlux(Sample->get_temperature()) << "\t" 
					<< ElectronHeatFlux(Sample->get_temperature());
	if( UseModel[5] )	ModelDataFile 	<< "\t" << NeutralFlux() << "\t" << NeutralHeatFlux();
	if( UseModel[6] )	ModelDataFile 	<< "\t" << NeutralRecombination(Sample->get_temperature());
	if( UseModel[7] )	ModelDataFile 	<< "\t" << SEE(Sample->get_temperature());
	if( UseModel[8] )	ModelDataFile 	<< "\t" << TEE(Sample->get_temperature());
	ModelDataFile << "\n";
}


// *************************************************** HEATING MODELS *************************************************** //

// http://users.wfu.edu/ucerkb/Nan242/L06-Vacuum_Evaporation.pdf, 
// https://en.wikipedia.org/wiki/Hertz%E2%80%93Knudsen_equation
// Using Hertz–Knudsen equation, returns Energy lost per second in Kila Joules
// MASS LOSS EQUATION 
// http://www.leb.eei.uni-erlangen.de/winterakademie/2006/result/content/course01/pdf/0102.pdf
const double HeatingModel::EvaporationModel(double DustTemperature)const{
	H_Debug("\tIn HeatingModel::EvaporationModel():\n\n");

	// Approximate emitted energy as mean maxwell
	// This used to be multiplied by 1000 which is why it was significant.
	double MaxwellEnergy = (3*Kb*DustTemperature/(2*1000)); // Converted to kJ.
// 	double MaxwellMeanVelocity = sqrt((8*Kb*DustTemperature)/(AtomicMass*1000*AMU));

	double EvapFlux = EvaporationFlux(DustTemperature);

	if( EvapFlux != EvapFlux ){ 
		std::cout << "\n\nError! EvapFlux = " << EvapFlux << "\n";
		throw std::exception(); 
	}
//	H1_Debug("\n\nIn HeatingModel::EvaporationModel()\nMaxwellEnergy = " << MaxwellEnergy << " kJ" << "\nEvapFlux = " 
//		<< EvapFlux << " s^-1\nBondEnergy = " << Sample->get_bondenergy()/AvNo << " kJ\nreturn = " 
//		<< EvapFlux*(MaxwellEnergy+Sample->get_bondenergy()/AvNo) << "\nDustTemperature = " << DustTemperature 
//		<<  "\nSample->get_surfacearea() = " << Sample->get_surfacearea()
//		<< "\nSample->get_vapourpressure() = " << Sample->get_vapourpressure());
	
	// See ElementData.h for more info on 'bondenergy'. Added to account for energy lost by breaking bonds.
	return EvapFlux*(MaxwellEnergy+Sample->get_bondenergy()/AvNo); 
}

const double HeatingModel::EvaporationFlux(double DustTemperature)const{
	H_Debug("\tIn HeatingModel::EvaporationFlux():\n\n");
	double AmbientPressure = 0;
	double StickCoeff = 1.0;
	return (StickCoeff*Sample->get_surfacearea()*AvNo*(Sample->get_vapourpressure()-AmbientPressure))/
			sqrt(2*PI*Sample->get_atomicmass()*R*DustTemperature);
}

// Using Stefan-Boltzmann Law, returns Energy lost per second in Kila Joules
const double HeatingModel::EmissivityModel(double DustTemperature)const{
//	H1_Debug("\n\nIn HeatingModel::EmissivityModel():");
//	H1_Debug("\nSample->get_surfacearea() = " << Sample->get_surfacearea() << "\nSigma = " << Sigma 
//		<< "\nDustTemperature = " << DustTemperature);
	// Energy emitted from a sample converted to kJ
	return Sample->get_emissivity()*Sample->get_surfacearea()*Sigma*(pow(DustTemperature,4)-pow(Pdata.AmbientTemp,4))/1000;
}

// VERY APPROXIMATE MODEL: Atmosphere assumed to be 300 degrees always, rough heat transfer coefficient is use
// https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling
const double HeatingModel::NewtonCooling(double DustTemperature)const{
	H_Debug("\tIn HeatingModel::NewtonCooling():\n\n");
	static bool runOnce = true;
	WarnOnce(runOnce,"In HeatingModel::NewtonCooling():\nHeatTransair Coefficient wrong for Tungsten, Beryllium and Graphite.");

	return (Sample->get_heattransair()*Sample->get_surfacearea()*(DustTemperature-Pdata.AmbientTemp))/1000; // convert to kJ
}	

// Neutral Recombination assuming Rn=0; fraction of backscattered ions/neutrals is zero.
const double HeatingModel::NeutralRecombination(double DustTemperature)const{
	H_Debug("\tIn HeatingModel::NeutralRecombination():\n\n");
//	H1_Debug("\n14.7*echarge - 2*Kb*DustTemperature*NeutralFlux()\n");
//	H1_Debug("\nIonFlux = " << IonFlux());
//	H1_Debug("\nReturn = " << (14.7*echarge*IonFlux()));
	double RN(0), RE(0);
//	backscatter(Pdata.ElectronTemp,Pdata.IonTemp,Mp,Sample->get_potential(),Element,RE,RN);
//	H1_Debug( "\nRN = " << RN );
	if( RN > 0.1 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"In HeatingModel::NeutralRecombination(double DustTemperature)\nRN > 0.1. Neutral Recombination affected by backscattering by more than 10%!");
	}

	return Sample->get_surfacearea()*(1-RN)*(14.7*echarge*IonFlux(DustTemperature) - 2*Kb*DustTemperature*NeutralFlux())/1000; // Convert from J to kJ
	return Sample->get_surfacearea()*(14.7*echarge*IonFlux(DustTemperature) - 2*Kb*DustTemperature*NeutralFlux())/1000; // Convert from J to kJ
}

const double HeatingModel::SEE(double DustTemperature)const{
	H_Debug("\tIn HeatingModel::SEE():\n\n");
//	H1_Debug("\neFlux*Sample->get_deltasec()*(3*echarge+Sample->get_bondenergy())\n");
//	H1_Debug("\neFlux=" << eFlux << "\nSample->get_deltasec() = " << Sample->get_deltasec() 
//			<< "\nSample->get_bondenergy() = " << Sample->get_bondenergy());
//	H1_Debug("\nReturn = " << eFlux*Sample->get_deltasec()*(3*echarge+Sample->get_bondenergy()));
	double SEE=0;
	if(ForceNegative || Sample->get_deltatot() <= 1 )
		SEE = Sample->get_surfacearea()*ElectronFlux(DustTemperature)*Sample->get_deltasec()*echarge*
                                        (3+Sample->get_workfunction())/1000; // Convert
	else if( Sample->get_deltatot() > 1 ) 	SEE = 0; // Electrons captured by positive grain
	return SEE;
}

const double HeatingModel::TEE(double DustTemperature)const{
	H_Debug("\tIn HeatingModel::TEE():\n\n");
//	H1_Debug("\neFlux*Sample->get_deltatherm()*(2*Kb*DustTemperature+Sample->get_bondenergy())\n");
//	H1_Debug("\neFlux=" << eFlux << "\nSample->get_deltatherm() = " << Sample->get_deltatherm() << "\n2*Kb*DustTemperature = "
//			<< 2*Kb*DustTemperature << "\nSample->get_bondenergy() = " << Sample->get_bondenergy());
//	H1_Debug("\nReturn = " << eFlux*Sample->get_deltatherm()*(2*Kb*DustTemperature+Sample->get_bondenergy()));
	double TEE=0;
	if( ForceNegative || Sample->get_deltatot() <= 1 )
		TEE = Sample->get_surfacearea()*Sample->get_deltatherm()*
			(2*Kb*DustTemperature+echarge*Sample->get_workfunction())/1000; // Convert to kJ
	else if( Sample->get_deltatot() > 1 ) 	TEE = 0; // Electrons captured by positive grain

	return TEE;
}

const double HeatingModel::IonHeatFlux(double DustTemperature)const{ // Assuming Re = 0
	H_Debug("\tIn HeatingModel::IonHeatFlux(double DustTemperature):\n\n");
//	H1_Debug("\nIonFlux() = " << IonFlux(DustTemperature) << "\nSample->get_potential() = " << Sample->get_potential());

	double RN(0), RE(0);
//	backscatter(Pdata.ElectronTemp,Pdata.IonTemp,Mp,Sample->get_potential(),Element,RE,RN);
	H1_Debug( "\nRE = " << RE );
	if( RE > 0.1 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"In HeatingModel::IonHeatFlux(double DustTemperature)\nRE > 0.1. Ion Heat Flux affected by backscattering by more than 10%!");
	}
//	H1_Debug("\nPdata.IonTemp = " << Pdata.IonTemp << "\nSample->get_potential() = " << Sample->get_potential() << "\nPdata.ElectronTemp = " << Pdata.ElectronTemp << "\n\n");
//	std::cout << "\n\n***** Ions = " << (Sample->get_surfacearea()*(1-RE)*IonFlux(DustTemperature)*Pdata.IonTemp*Kb/1000)
//	*(2+2*Sample->get_potential()*(Pdata.ElectronTemp/Pdata.IonTemp)+pow(Sample->get_potential()*(Pdata.ElectronTemp/Pdata.IonTemp),2))/(1+Sample->get_potential()*(Pdata.ElectronTemp/Pdata.IonTemp));

	return (Sample->get_surfacearea()*(1-RE)*IonFlux(DustTemperature)*Pdata.IonTemp*Kb/1000) // Convert from Joules to KJ
	*(2+2*Sample->get_potential()*(Pdata.ElectronTemp/Pdata.IonTemp)+pow(Sample->get_potential()*(Pdata.ElectronTemp/Pdata.IonTemp),2))/(1+Sample->get_potential()*(Pdata.ElectronTemp/Pdata.IonTemp)); 

}

const double HeatingModel::ElectronHeatFlux(double DustTemperature)const{ // Only for a negative grain
	H_Debug("\tIn HeatingModel::ElectronHeatFlux():\n\n");
//	H1_Debug("\nSample->get_surfacearea()*2*ElectronFlux()*Pdata.ElectronTemp*Kb/1000\n");
//	std::cout << "\n\n***** Electrons = " 
//		<< Sample->get_surfacearea()*2*ElectronFlux(DustTemperature)*Pdata.ElectronTemp*Kb/1000;
	return Sample->get_surfacearea()*2*ElectronFlux(DustTemperature)*Pdata.ElectronTemp*Kb/1000; // Convert from Joules to KJ
}

const double HeatingModel::NeutralHeatFlux()const{
	H_Debug("\tIn HeatingModel::NeutralHeatFlux():\n\n");
//	H1_Debug("\nNeutralFlux()*Pdata.NeutralTemp*Kb*2\n");
	return Sample->get_surfacearea()*2*NeutralFlux()*Pdata.NeutralTemp*Kb/1000; // Convert from Joules to KJ
}

const double HeatingModel::IonFlux(double DustTemperature)const{
//	H_Debug("\n\nIn HeatingModel::IonFlux():");
//	H1_Debug("\nElectronFlux()*(1-Sample->get_deltatherm()-Sample->get_deltasec())\n");
//	H1_Debug("\nElectronFlux() =" << ElectronFlux() << "\nSample->get_deltatherm() = " << Sample->get_deltatherm() 
//			<< "\nSample->get_deltasec() = " << Sample->get_deltasec());
//	H1_Debug("\nReturn = " << ElectronFlux()*(1-Sample->get_deltatherm()-Sample->get_deltasec()));
	double IonFlux=0;

	if(ForceNegative || Sample->get_deltatot() > 1 )IonFlux = ElectronFlux(DustTemperature); //Positive grain, DeltaTot() > 1
	else	IonFlux = ElectronFlux(DustTemperature)*(1-Sample->get_deltatot());
	return IonFlux;
}

const double HeatingModel::ElectronFlux(double DustTemperature)const{
	H_Debug("\tIn HeatingModel::ElectronFlux():\n\n");
//	H1_Debug("\nPdata.NeutralDensity*exp(-Sample->get_potential())*sqrt(Kb*DustTemperature/(2*PI*Me))\n");
//	H1_Debug("\nexp(-Sample->get_potential()) = " << exp(-Sample->get_potential())
//			<< "\n(2*PI*Me) = " << (2*PI*Me) << "sqrt(Kb*DustTemperature/(2*PI*Me)) = " 
//			<< sqrt(Kb*DustTemperature/(2*PI*Me)));
//	H1_Debug("\nReturn = " << Pdata.NeutralDensity*exp(-Sample->get_potential())*sqrt(Kb*DustTemperature/(2*PI*Me)));
//	std::cout << "\neFlux = " << Pdata.ElectronDensity*exp(-Sample->get_potential())*sqrt(Kb*Pdata.ElectronTemp/(2*PI*Me)); std::cin.get();
	return Pdata.ElectronDensity*exp(-Sample->get_potential())*sqrt(Kb*Pdata.ElectronTemp/(2*PI*Me));
}

const double HeatingModel::NeutralFlux()const{
	H_Debug("\tIn HeatingModel::NeutralFlux():\n\n");
//	H1_Debug("\nPdata.NeutralTemp*sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mp))/4\n");
//	H1_Debug("(2*PI*Mp) = " << (2*PI*Mp) << "\nsqrt(Kb*Pdata.NeutralTemp/(2*PI*Mp)) = "<< sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mp)));
//	H1_Debug("\nReturn = " << 2*Pdata.NeutralTemp*Pdata.NeutralDensity*sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mp)));
	return Pdata.NeutralDensity*sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mp));
}

// ************************************* //
