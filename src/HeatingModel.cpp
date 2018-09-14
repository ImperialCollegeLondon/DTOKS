//#define PAUSE
//#define HEATING_DEBUG
//#define HEATING_DEEP_DEBUG

#include "HeatingModel.h"
#include "Constants.h"
#include "Functions.h"

HeatingModel::HeatingModel():Model(){
	H_Debug("\n\nIn HeatingModel::HeatingModel():Model()\n\n");
	Defaults();
	CreateFile("Default_Heating_filename.txt",false);
}


// Constructor which specifies dust radius, by passing a pointer to a Matter reference.
HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models,
				Matter *& sample, PlasmaData &pdata) : Model(sample,pdata,accuracy){
	H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models, Matter *& sample, PlasmaData const *&pdata) : Model(sample,pdata,accuracy)\n\n");
	Defaults();
	UseModel 		= models;
	CreateFile(filename,false);
}

// Constructor which specifies dust radius, by passing a pointer to a Matter reference.
HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models,
				Matter *& sample, PlasmaData *pdata) : Model(sample,*pdata,accuracy){
	H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models, Matter *& sample, PlasmaData const *&pdata) : Model(sample,pdata,accuracy)\n\n");
	Defaults();
	UseModel 		= models;
	CreateFile(filename,false);
}

HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models,
				Matter *& sample, PlasmaGrid_Data &pgrid) : Model(sample,pgrid,accuracy){
	H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models, Matter *& sample, PlasmaGrid const &pgrid) : Model(sample,pgrid,accuracy)\n\n");
	Defaults();
	UseModel 	= models;
	CreateFile(filename,false);
}

HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models,
				Matter *& sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata) 
				: Model(sample,pgrid,pdata,accuracy){
	H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, float accuracy, std::array<bool,HMN> &models, Matter *& sample, PlasmaGrid const &pgrid) : Model(sample,pgrid,accuracy)\n\n");
	Defaults();
	UseModel 	= models;
	CreateFile(filename,false);
}

void HeatingModel::Defaults(){
	H_Debug("\tIn HeatingModel::Defaults()\n\n");
	UseModel = {false,false,false,false,false,false,false};
	PowerIncident = 0;			// kW, Power Incident
	OldTemp = Sample->get_temperature();	// Set default OldTemp
	ThermalEquilibrium = false;
	RE = 0.0;
	RN = 0.0;
}

void HeatingModel::CreateFile(std::string filename, bool PrintPhaseData){
	H_Debug("\tIn HeatingModel::CreateFile(std::string filename, bool PrintPhaseData)\n\n");
	FileName=filename;
	ModelDataFile.open(FileName);
	ModelDataFile << std::fixed << std::setprecision(16) << std::endl;
	ModelDataFile << "Time\tTemp\tMass\tDensity";
	if( PrintPhaseData ) 				ModelDataFile << "\tFusionE\tVapourE";
	if( Sample->get_c(0) == 'v' || Sample->get_c(0) == 'V' ) 	ModelDataFile << "\tCv";
	if( Sample->get_c(1) == 'v' || Sample->get_c(1) == 'V' ) 	ModelDataFile << "\tVapourP";
	if( Sample->get_c(2) == 'v' || Sample->get_c(2) == 'V' ) 	ModelDataFile << "\tLinearExpansion";
	if( UseModel[0] )					ModelDataFile << "\tEmissLoss";
	if( UseModel[0] && Sample->get_c(3) == 'v' )   			ModelDataFile << "\tEmissiv";
	if( UseModel[1] )					ModelDataFile << "\tEvapRate\tEvapLoss\tEvapMassLoss";
	if( UseModel[2] )					ModelDataFile << "\tNewton";
	if( UseModel[3] )					ModelDataFile << "\tIonFlux\tIonHeatFlux";
	if( UseModel[4] )					ModelDataFile << "\tElectronFlux\tElectronHeatFlux";
	if( UseModel[5] )					ModelDataFile << "\tNeutralFlux\tNeutralHeatFlux";
	if( UseModel[6] )					ModelDataFile << "\tNeutralRecomb";
	if( UseModel[7] )					ModelDataFile << "\tSEE";
	if( UseModel[8] )					ModelDataFile << "\tTEE";
	ModelDataFile << "\n";
	Print();
	ModelDataFile.close();
	ModelDataFile.clear();
}

// This model forces the time step to be the value which produces a change in temperature or 1*accuracy degree
double HeatingModel::ProbeTimeStep()const{
	H_Debug( "\tIn HeatingModel::ProbeTimeStep()\n\n" );

	// Take Eularian step to get initial time step
	H_Debug("\t"); double TotalPower = CalculatePower(Sample->get_temperature());
	double timestep = fabs((Sample->get_mass()*Sample->get_heatcapacity()*Accuracy)/TotalPower);

	// Calculate timestep that produces mass change of less than 0.01% of current mass.
	// If this timestep is quicker than current step, change timestep
	if( TotalPower != 0 ){
		if( UseModel[1] && Sample->is_liquid() ){
			double MassTimeStep = (0.01*Sample->get_mass()*AvNo)
						/(EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass());
			if( MassTimeStep < timestep ){
				H_Debug("\nMass is limiting time step\nMassTimeStep = " << MassTimeStep << "\ntimestep = " << timestep);
				timestep = MassTimeStep;
			}
		}
	}else{
		// Check thermal equilibrium hasn't been explicitly reached somehow.
		if( ContinuousPlasma ){	
			static bool runOnce = true;
			WarnOnce(runOnce,"\nWarning! TotalPower = 0");
			std::cout << "\nThermalEquilibrium reached on condition (1): TotalPower = 0.";
			timestep = 1;
		}else if( !ContinuousPlasma ){
			std::cout << "\nNo Net Power Region...";
			timestep = 10;
		}
	}
	// Check Thermal Equilibrium hasn't been reached for continuous plasma
	double DeltaTempTest = TotalPower*timestep/(Sample->get_mass()*Sample->get_heatcapacity());
	// If we're not boiling	and in a continuous Plasma
	if( Sample->get_temperature() != Sample->get_superboilingtemp() && 
		Sample->get_temperature() != Sample->get_meltingtemp() && ContinuousPlasma ){
// 			&& (!UseModel[1] && Sample->is_liquid()) ){
		static bool runOnce = true;
		WarnOnce(runOnce,"THERMAL EQUILIBRIUM IN CONTINUOUS PLASMA MAY NOT WORK!\nOldTemp is only redefined inside UpdateTimeStep, not ProbeTimeStep! CAUTION!");
		if( ((Sample->get_temperature()-OldTemp > 0 && DeltaTempTest < 0) // If temperature changed sign this step
			|| (Sample->get_temperature()-OldTemp < 0 && DeltaTempTest > 0)) ){
			std::cout << "\n\nThermal Equilibrium reached on condition (2): Sign change of Temperature change!";
			timestep = 1;	
		}if( (fabs(DeltaTempTest/timestep) < 0.01) ){ // If Temperature gradient is less than 1%
			std::cout << "\n\nThermal Equilibrium reached on condition (3): Temperature Gradient < 0.01!";
			timestep = 1;  
		}
	}

	H1_Debug("\nSample->get_mass() = " << Sample->get_mass() << "\nSample->get_heatcapacity() = " << 
		Sample->get_heatcapacity() << "\nTotalPower = " << TotalPower << "\nAccuracy = " << Accuracy);
	assert(timestep > 0 && timestep != INFINITY && timestep == timestep);

	return timestep;
}

// This model forces the time step to be the value which produces a change in temperature or 1*accuracy degree
double HeatingModel::UpdateTimeStep(){
	H_Debug( "\tIn HeatingModel::UpdateTimeStep()\n\n" );

	TimeStep = ProbeTimeStep();
	if( TimeStep == 1 ) 	ThermalEquilibrium = true;
	OldTemp = Sample->get_temperature();

	return TimeStep;
}

void HeatingModel::Print(){
	H_Debug("\tIn HeatingModel::Print()\n\n");
	ModelDataFile.open(FileName,std::ofstream::app);
	ModelDataFile 	<< TotalTime << "\t" << Sample->get_temperature() << "\t" << Sample->get_mass() 
		<< "\t" << Sample->get_density();

//	bool PrintPhaseData = false; // Lol
//	if( PrintPhaseData )	ModelDataFile 	<< "\t" << Sample->get_fusionenergy() << "\t" << Sample->get_vapourenergy();

	if( UseModel[0] ) 	ModelDataFile 	<< "\t" << EmissivityModel(Sample->get_temperature());
	if( UseModel[0] && Sample->get_c(3) == 'v' )
				ModelDataFile	<< "\t" << Sample->get_emissivity();
	if( UseModel[1] && Sample->is_liquid() )	
				ModelDataFile 	<< "\t" << EvaporationFlux(Sample->get_temperature()) 
					<< "\t" << EvaporationModel(Sample->get_temperature())*1000
					<< "\t" << EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass()/AvNo;
	else if( UseModel[1] )	ModelDataFile 	<< "\t" << 0 << "\t" << 0 << "\t" << 0;
	if( UseModel[2] )	ModelDataFile 	<< "\t" << NewtonCooling(Sample->get_temperature());
	if( UseModel[3] )	ModelDataFile 	<< "\t" << IonFlux(Sample->get_temperature()) 
					<< "\t" << DTOKSIonHeatFlux(Sample->get_temperature());
	if( UseModel[4] )	ModelDataFile 	<< "\t" << ElectronFlux(Sample->get_temperature()) << "\t" 
					<< ElectronHeatFlux(Sample->get_temperature());
	if( UseModel[5] )	ModelDataFile 	<< "\t" << NeutralFlux() << "\t" << NeutralHeatFlux(Sample->get_temperature());
	if( UseModel[6] )	ModelDataFile 	<< "\t" << NeutralRecombination(Sample->get_temperature());
	if( UseModel[7] )	ModelDataFile 	<< "\t" << SEE(Sample->get_temperature());
	if( UseModel[8] )	ModelDataFile 	<< "\t" << TEE(Sample->get_temperature());
	ModelDataFile << "\n";
	ModelDataFile.close();
	ModelDataFile.clear();
}



const int HeatingModel::Vapourise(){
	H_Debug("\tIn HeatingModel::Vapourise()\n\n");
	// If the sample is gaseous or in TE (Given that the plasma is continuous), the model ends.
	while( !Sample->is_gas() ){
		Heat();
		UpdateTimeStep();
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

	assert( Sample->get_mass() > 0 );

	//backscatter(Pdata->ElectronTemp,Pdata->IonTemp,Pdata->mi,Sample->get_potential(),Sample->get_elem(),RE,RN);
	double TotalEnergy = RungeKutta4(timestep);// Calculate total energy through RungeKutta4 method
	H1_Debug( "\tTotalEnergy = " << TotalEnergy << "\n");
	Sample->update_temperature(TotalEnergy);  // Update Temperature

//	std::cout << "\nTemperature = " << Sample->get_temperature();
	// Account for evaporative mass loss, if model is turned on, if it's a liquid and not boiling!
	if( UseModel[1] && Sample->is_liquid() && (Sample->get_temperature() != Sample->get_boilingtemp()) )
		Sample->update_mass( (timestep*EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass())/AvNo );
	H1_Debug("\t\tMass Loss = " << (timestep*EvaporationFlux(Sample->get_temperature())*Sample->get_atomicmass())/AvNo << "\n");
	if( !Sample->is_gas() )
		Sample->update();
	Print();  // Print data to file
	H_Debug("\t"); 

	TotalTime += timestep;
}

void HeatingModel::Heat(){
	H_Debug("\tIn HeatingModel::Heat()\n\n");
	Heat(TimeStep);
}

double HeatingModel::CalculatePower(double DustTemperature)const{
	H_Debug( "\tIn HeatingModel::CalculatePower(double DustTemperature = " << DustTemperature << ")\n\n");
	double TotalPower = PowerIncident*1000; // This looks weird, but doing operations like this reduces the number of divisions
	
	if( UseModel[0] )				TotalPower -= EmissivityModel		(DustTemperature);
	if( UseModel[1] && Sample->is_liquid() )	TotalPower -= EvaporationModel		(DustTemperature)*1000;
	if( UseModel[2] )				TotalPower -= NewtonCooling			(DustTemperature);
	if( UseModel[3] )				TotalPower += DTOKSIonHeatFlux		(DustTemperature);
	if( UseModel[4] )				TotalPower += ElectronHeatFlux		(DustTemperature);
	if( UseModel[5] )				TotalPower += NeutralHeatFlux		(DustTemperature);
	if( UseModel[6] )				TotalPower += NeutralRecombination	(DustTemperature);
	if( UseModel[7] )				TotalPower -= SEE					(DustTemperature);
	if( UseModel[8] )				TotalPower -= TEE					(DustTemperature);	
	if( UseModel[9] )				TotalPower += DUSTTIonHeatFlux		(DustTemperature);	
	TotalPower = TotalPower/1000;

	H1_Debug("\n\n\tPowerIncident = \t" 	<< PowerIncident*1000				<< "W");
	H1_Debug("\n\tEmissivityModel() = \t" 	<< -EmissivityModel(DustTemperature) 		<< "W");
	H1_Debug("\n\tEvaporationModel() = \t" 	<< -EvaporationModel(DustTemperature)*1000 	<< "W");
	H1_Debug("\n\tNewtonCooling() = \t" 	<< -NewtonCooling(DustTemperature) 		<< "W");
	H1_Debug("\n\tDTOKSIonHeatFlux() = \t" 	<< DTOKSIonHeatFlux(DustTemperature) 		<< "W");
	H1_Debug("\n\tElectronHeatFlux() = \t" 	<< ElectronHeatFlux(DustTemperature) 		<< "W");
	H1_Debug("\n\tNeutralHeatFlux(DustTemperature) = \t" 	<< NeutralHeatFlux(DustTemperature) 				<< "W");
	H1_Debug("\n\tNeutralRecombination() = \t"<< NeutralRecombination(DustTemperature)	<< "W");
	H1_Debug("\n\tSEE() = \t" 		<< -SEE(DustTemperature) 			<< "W");
	H1_Debug("\n\tTEE() = \t" 		<< -TEE(DustTemperature) 			<< "W\n");
	//std::cin.get();
	H1_Debug("\nTotalPower = \t" 		<< TotalPower*1000		<< "W\n");
	return TotalPower;
}

double HeatingModel::RungeKutta4(double timestep){
	H_Debug( "\tIn HeatingModel::RungeKutta4()\n\n");
	double k1 = CalculatePower(Sample->get_temperature()); 
	if(k1<0 && fabs(k1/2) > Sample->get_temperature()){
		std::cout << "\n\nThermal Equilibrium reached on condition (4): k1 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return 0;
	}
	double k2 = CalculatePower(Sample->get_temperature()+k1/2); 
	if( k2<0 && fabs(k2/2) > Sample->get_temperature() ){
		std::cout << "\n\nThermal Equilibrium reached on condition (4): k2 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return (timestep/6)*k1;
	}
	double k3 = CalculatePower(Sample->get_temperature()+k2/2);
	if( k3<0 && fabs(k3) > Sample->get_temperature() ){
		std::cout << "\n\nThermal Equilibrium reached on condition (4): k3 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return (timestep/6)*(k1+2*k2);
	}
	double k4 = CalculatePower(Sample->get_temperature()+k3);
	if( k4<0 && fabs(k4/2) > Sample->get_temperature() ){
		std::cout << "\n\nThermal Equilibrium reached on condition (4): k4 step negative and larger than Td!";
		ThermalEquilibrium = true;
		return (timestep/6)*(k1+2*k2+2*k3);
	}
	H1_Debug( "\ntimestep = " << timestep << "\nk1 = " << k1 << "\nk2 =" << k2 << "\nk3 = " << k3 << "\nk4 = " << k4);
	return (timestep/6)*(k1+2*k2+2*k3+k4);
};


// *************************************************** HEATING MODELS *************************************************** //

// Using Stefan-Boltzmann Law, returns Energy lost per second in Kila Joules
const double HeatingModel::EmissivityModel(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::EmissivityModel(double DustTemperature):");
	// Energy emitted from a sample converted to kJ
	return Sample->get_emissivity()*Sample->get_surfacearea()*Sigma
			*(pow(DustTemperature,4)-pow(Pdata->AmbientTemp,4));
}

// http://users.wfu.edu/ucerkb/Nan242/L06-Vacuum_Evaporation.pdf, 
// https://en.wikipedia.org/wiki/Hertz%E2%80%93Knudsen_equation
// Using Hertz–Knudsen equation, returns Energy lost per second in Kila Joules
// MASS LOSS EQUATION 
// http://www.leb.eei.uni-erlangen.de/winterakademie/2006/result/content/course01/pdf/0102.pdf
const double HeatingModel::EvaporationModel(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::EvaporationModel():\n\n");

	// Approximate emitted energy as mean maxwell
	// This used to be multiplied by 1000 which is why it was significant.
	double MaxwellEnergy = (3*Kb*DustTemperature/(2*1000)); // Converted to kJ.
// 	double MaxwellMeanVelocity = sqrt((8*Kb*DustTemperature)/(AtomicMass*1000*AMU));

	double EvapFlux = EvaporationFlux(DustTemperature);

	if( EvapFlux != EvapFlux ){ 
		std::cout << "\n\nError! EvapFlux = " << EvapFlux << "\n";
		throw std::exception(); 
	}

	// See GrainStructs.h for more info on 'bondenergy'. 
	// Added to account for energy lost by breaking bonds.
	return EvapFlux*(MaxwellEnergy+Sample->get_bondenergy()/AvNo); 
}

const double HeatingModel::EvaporationFlux(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::EvaporationFlux():\n\n");
	double AmbientPressure = 0;
	double StickCoeff = 1.0;
	return (StickCoeff*Sample->get_surfacearea()*AvNo*(Sample->get_vapourpressure()-AmbientPressure))/
			sqrt(2*PI*Sample->get_atomicmass()*R*DustTemperature);
}

// VERY APPROXIMATE MODEL: Atmosphere assumed to be 300 degrees always, rough heat transfer coefficient is use
// https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling
const double HeatingModel::NewtonCooling(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::NewtonCooling():\n\n");
	static bool runOnce = true;
	WarnOnce(runOnce,"In HeatingModel::NewtonCooling():\nHeatTransair Coefficient wrong for Tungsten, Beryllium and Graphite.");

	return (Sample->get_heattransair()*Sample->get_surfacearea()*(DustTemperature-Pdata->AmbientTemp)); 
}	

// Neutral Recombination assuming Rn=0; fraction of backscattered ions/neutrals is zero.
const double HeatingModel::NeutralRecombination(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::NeutralRecombination():\n\n");

//	H1_Debug( "\nRN = " << RN );
//	if( RN > 0.1 ){ // Uncomment when RN is calculated
//		static bool runOnce = true;
//		WarnOnce(runOnce,"In HeatingModel::NeutralRecombination(double DustTemperature)\nRN > 0.1. Neutral Recombination affected by backscattering by more than 10%!");
//	}
	if( Sample->is_positive() ) 	
		return Sample->get_surfacearea()*(14.7*echarge - 2.0*Kb*DustTemperature)*IonFlux(DustTemperature); 
	else				
		return Sample->get_surfacearea()*(1.0-RN)*(14.7*echarge - 2.0*Kb*DustTemperature)*IonFlux(DustTemperature); 
}

const double HeatingModel::SEE(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::SEE():\n\n");

	double SEE=0; 
	if( !Sample->is_positive() ) // Dust grain is negative
		SEE = Sample->get_surfacearea()*ElectronFlux(DustTemperature)*Sample->get_deltasec()*echarge*
     (3.0+Sample->get_workfunction()); 
	// else, electrons captured by positive grain, return zero
	return SEE;
}

const double HeatingModel::TEE(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::TEE():\n\n");

	double TEE=0;
	if( !Sample->is_positive() ) // Dust grain is negative
		TEE = Sample->get_surfacearea()*Sample->get_deltatherm()*ElectronFlux(Sample->get_temperature())*
			(2*Kb*DustTemperature+echarge*Sample->get_workfunction()); 
	// else, electrons captured by positive grain, return zero

	return TEE;
}

const double HeatingModel::DTOKSIonHeatFlux(double DustTemperature)const{ // Assuming Re = 0
	H_Debug("\n\tIn HeatingModel::DTOKSIonHeatFlux(double DustTemperature):\n\n");

//	H1_Debug( "\nRE = " << RE );
//	if( RE > 0.1 ){ // Uncomment when RE is calculated
//		static bool runOnce = true;
//		WarnOnce(runOnce,"In HeatingModel::DTOKSIonHeatFlux(double DustTemperature)\nRE > 0.1. Ion Heat Flux affected by backscattering by more than 10%!");
//	}

	return (Sample->get_surfacearea()*(1-RE)*IonFlux(DustTemperature)*Pdata->IonTemp*Kb) 
				*(2+2*Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp)
				+pow(Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp),2))
				/(1+Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp)); 
				// Convert from Joules to KJ
}

// This is the heat flux of ions incident on the dust grain surface ONLY
// Calculated from equation (31) & (32), page 29, from the following reference:
// S.I. Krasheninnikov, R.D. Smirnov, and D.L. Rudakov,
// Plasma Phys. Control. Fusion 53, 083001 (2011).
const double HeatingModel::DUSTTIonHeatFlux(double DustTemperature)const{ // Assuming Re = 0
	H_Debug("\n\tIn HeatingModel::DUSTTIonHeatFlux(double DustTemperature):\n\n");

//	H1_Debug( "\nRE = " << RE );
//	if( RE > 0.1 ){ // Uncomment when RE is calculated
//		static bool runOnce = true;
//		WarnOnce(runOnce,"In HeatingModel::DUSTTIonHeatFlux(double DustTemperature)\nRE > 0.1. Ion Heat Flux affected by backscattering by more than 10%!");
//	}
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double RelativeVelocity = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	double RenormalisedPotential = (Pdata->Z*Sample->get_potential())/TemperatureRatio;
	double Term1(0.0);
	double Term2(0.0);
	double Term2Coeff = (1.0/(4.0*RelativeVelocity))*(3.0+12.0*RelativeVelocity*RelativeVelocity
					+4.0*pow(RelativeVelocity,4.0)-2.0*RenormalisedPotential*(1+2.0*RelativeVelocity*RelativeVelocity));
	if( !Sample->is_positive() ){ // Negative dust grain
		Term1 =(1.0/(2.0*sqrt(PI)))*(5.0+2.0*RelativeVelocity*RelativeVelocity
				-2.0*RenormalisedPotential)*exp(-RelativeVelocity*RelativeVelocity); 
		Term2 = erf(RelativeVelocity)*Term2Coeff; 
	}else{ // Positive dust grain
		assert( Sample->get_potential() <= 0.0 );
		double Coeff1 = 5.0+2.0*RelativeVelocity*RelativeVelocity
				-((3.0+2.0*RelativeVelocity*RelativeVelocity)/RelativeVelocity)*sqrt(-RenormalisedPotential);
		double Coeff2 = 5.0+2.0*RelativeVelocity*RelativeVelocity
				+((3.0+2.0*RelativeVelocity*RelativeVelocity)/RelativeVelocity)*sqrt(-RenormalisedPotential);
		double uip = RelativeVelocity+sqrt(-RenormalisedPotential);
		double uim = RelativeVelocity-sqrt(-RenormalisedPotential);
		Term1 = (1.0/(4.0*sqrt(PI)))*(Coeff1*exp(-uip*uip) + Coeff2*exp(-uim*uim)); 
		Term2 = (Term2Coeff/2.0)*(erf(uip)+erf(uim));
	}
	
	assert((Term1 + Term2) > 0 && (Term1 + Term2) != INFINITY && (Term1 + Term2) == (Term1 + Term2));
	return Sample->get_surfacearea()*(1.0-RE)*Pdata->IonDensity
				*Pdata->IonTemp*Kb*IonThermalVelocity*(Term1 + Term2);
				// Convert from Joules to KJ
}


const double HeatingModel::ElectronHeatFlux(double DustTemperature)const{ // Only for a negative grain
	H_Debug("\n\tIn HeatingModel::ElectronHeatFlux():\n\n");
	if( !Sample->is_positive() ){ // Negative dust grain
		return 2.0*sqrt(2.0*PI)*Sample->get_radius()*Sample->get_radius()
				*ElectronFlux(DustTemperature)*Kb*(Pdata->ElectronTemp-DustTemperature);
	}else{ // Positive dust grain
		return 2.0*sqrt(2.0*PI)*Sample->get_radius()*Sample->get_radius()
				*ElectronFlux(DustTemperature)*Kb*(Pdata->ElectronTemp-DustTemperature);
	}
}

const double HeatingModel::NeutralHeatFlux(double DustTemperature)const{
	H_Debug("\n\tIn HeatingModel::NeutralHeatFlux(double DustTemperature):\n\n");
	return Sample->get_surfacearea()*2*NeutralFlux()*(Pdata->NeutralTemp-DustTemperature)*Kb;
}

// ************************************* //
