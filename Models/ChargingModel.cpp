//#define PAUSE
//#define CHARGING_DEBUG

#include "ChargingModel.h"

ChargingModel::ChargingModel():Model(){
	C_Debug("\n\nIn ChargingModel::ChargingModel():Model()\n\n");
	UseModel[0] = true;				// Charging Models turned on of possibly 9
	CreateFile("Default_Charging_Filename.txt");
}

ChargingModel::ChargingModel(std::string filename, double accuracy, std::array<bool,1> models,
				Matter *& sample, PlasmaData *&pdata) : Model(sample,pdata,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename,double accuracy,std::array<bool,1> models,Matter *& sample, PlasmaData const& pdata) : Model(sample,pdata,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, double accuracy, std::array<bool,1> models,
				Matter *& sample, PlasmaGrid &pgrid) : Model(sample,pgrid,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename,double accuracy,std::array<bool,1> models,Matter *& sample, PlasmaGrid const& pgrid) : Model(sample,pgrid,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}

void ChargingModel::CreateFile(std::string filename){
	C_Debug("\tIn ChargingModel::CreateFile(std::string filename)\n\n");
	ModelDataFile.open(filename);
	if( UseModel[0] ) ModelDataFile << "Positive\tPotential";


	ModelDataFile << "\n";
}

void ChargingModel::Print(){
	C_Debug("\tIn ChargingModel::Print()\n\n");
	if( Sample->is_positive() )  ModelDataFile << "Pos\t";
	if( !Sample->is_positive() ) ModelDataFile << "Neg\t";
	if( UseModel[0] ) ModelDataFile << Sample->get_potential();
	ModelDataFile << "\n";
}

double ChargingModel::ProbeTimeStep()const{
	C_Debug( "\tIn ChargingModel::ProbeTimeStep()\n\n" );
	double timestep(0);	

	if( Pdata->ElectronDensity != 0 && Pdata->IonTemp != 0 ){
		// Calcualte the time scale of the behaviour from Krashinnenikovs equation
		double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));
		double PlasmaFreq = sqrt((Pdata->ElectronDensity*pow(echarge,2))/(epsilon0*Me));
		timestep = sqrt(2*PI) * ((DebyeLength)/Sample->get_radius()) 
				* (1/(PlasmaFreq*(1+Pdata->ElectronTemp/Pdata->IonTemp+fabs(Sample->get_potential()))))*Accuracy;
		C_Debug("\n\t\tDebyeLength = " << DebyeLength << "\n\t\tPlasmaFreq = " << PlasmaFreq 
			<< "\n\t\ttimestep = " << timestep << "\n\n");

	}else{	timestep = 1; } // In region of no plasma


	assert(timestep == timestep);
	assert(timestep > 0);

	return timestep;
}

double ChargingModel::UpdateTimeStep(){
	C_Debug( "\tIn ChargingModel::UpdateTimeStep()\n\n" );
	

	if( Pdata->ElectronDensity != 0 && Pdata->IonTemp != 0 ){
		// Calcualte the time scale of the behaviour from Krashinnenikovs equation

		double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));
		double PlasmaFreq = sqrt((Pdata->ElectronDensity*pow(echarge,2))/(epsilon0*Me));
		TimeStep = sqrt(2*PI) * ((DebyeLength)/Sample->get_radius()) 
				* (1/(PlasmaFreq*(1+Pdata->ElectronTemp/Pdata->IonTemp+fabs(Sample->get_potential()))))*Accuracy;
		C_Debug("\n\t\tDebyeLength = " << DebyeLength << "\n\t\tPlasmaFreq = " << PlasmaFreq 
			<< "\n\t\tElectronTemp = " << Pdata->ElectronTemp << "\n\t\tTimeStep = " << TimeStep << "\n\n");

	}else{	TimeStep = 1; } // In region of no plasma


	assert(TimeStep == TimeStep);
	assert(TimeStep > 0);

	return TimeStep;
}

void ChargingModel::Charge(double timestep){
	C_Debug("\tIn ChargingModel::Charge(double timestep)\n\n");

	// Make sure timestep input time is valid. Shouldn't exceed the timescale of the process.
	assert(timestep > 0);// && timestep <= TimeStep );
	double DSec = DeltaSec();
	double DTherm = DeltaTherm();
	// Assume the grain is negative and calculate potential
	if( UseModel[0] ){
		double Potential;
		if( (DSec + DTherm) < 1.0 ){ // solveOML only defined for deltatot < 1.0
			Potential = solveOML( DSec + DTherm,Sample->get_potential()); 
		}else{ // If the grain is in fact positive ...
			Potential = solveOML( DSec + DTherm,Sample->get_potential());
			if( Potential < 0.0 ){
				Potential = solveOML(0.0,Sample->get_potential())-Kb*Sample->get_temperature()
						/(echarge*Pdata->ElectronTemp);
			}
		}
		double ChargeOfGrain = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*Pdata->ElectronTemp)/echarge;
		Sample->update_charge(ChargeOfGrain,Potential,DSec,DTherm);
//		std::cout << "\nPotential = " << Potential << "\nDeltaSec = " << Sample->get_deltasec() << "\nDeltatherm = " 
//			<< Sample->get_deltatherm() << "\nCharge = " << ChargeOfGrain; std::cin.get();
	}
	TotalTime += timestep;

	C_Debug("\t"); Print();
}

void ChargingModel::Charge(){
	C_Debug("\tIn ChargingModel::Charge()\n\n");

	double DSec = DeltaSec();
	double DTherm = DeltaTherm();
	// Assume the grain is negative and calculate potential
	if( UseModel[0] ){
		double Potential;
		if( DSec+DTherm < 1.0 ){ // solveOML only defined for deltatot < 1.0
			Potential = solveOML(DSec+DTherm,Sample->get_potential()); 
		}else{ // If the grain is in fact positive ...
			Potential = solveOML(DSec+DTherm,Sample->get_potential());
			if( Potential < 0.0 ){
				Potential = solveOML(0.0,Sample->get_potential())-Kb*Sample->get_temperature()
						/(echarge*Pdata->ElectronTemp);
			}
		}
		double ChargeOfGrain = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*Pdata->ElectronTemp)/echarge;
		Sample->update_charge(ChargeOfGrain,Potential,DSec,DTherm);
//		std::cout << "\nPotential = " << Potential << "\nDeltaSec = " << Sample->get_deltasec() << "\nDeltatherm = " 
//			<< Sample->get_deltatherm() << "\nCharge = " << ChargeOfGrain; std::cin.get();
	}
	TotalTime += TimeStep;

	C_Debug("\t"); Print();
}

double ChargingModel::solveOML(double a, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
        if( a >= 1.0 ){
		static bool runOnce;
		WarnOnce(runOnce,"DeltaTot >= 1.0. DeltaTot being set equal to unity.");
		a = 1.0;
	}
	double b(0);
	if( Pdata->ElectronTemp == 0 ) 	b = 0;
	else 				b = Pdata->IonTemp/Pdata->ElectronTemp;
	double C = Me/Mp;

	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/((a-1)*exp(-guess) - sqrt(C/b) ) );

	while(fabs(guess-x1)>1e-2){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}

double ChargingModel::DeltaSec()const{
	C_Debug("\tIn ChargingModel::DeltaSec()\n\n");
	return (Richardson*pow(Sample->get_temperature(),2)*exp(-(Sample->get_workfunction()*echarge)
				/(Kb*Sample->get_temperature())))/echarge;
}

double ChargingModel::DeltaTherm()const{
	C_Debug("\tIn ChargingModel::DeltaTherm()\n\n");
	return sec(Pdata->ElectronTemp/1.16e5,Sample->get_elem());
}
