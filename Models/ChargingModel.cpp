//#define PAUSE
//#define CHARGING_DEBUG

#include "ChargingModel.h"

ChargingModel::ChargingModel():Model(){
	C_Debug("\n\nIn ChargingModel::ChargingModel():Model()\n\n");
	UseModel[0] = true;				// Charging Models turned on of possibly 9
	CreateFile("Default_Charging_Filename.txt");
	TimeStep = 0;
}

ChargingModel::ChargingModel(std::string filename, double accuracy, std::array<bool,1> models,
				Matter *& sample, PlasmaData &pdata) : Model(sample,pdata,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename,double accuracy,std::array<bool,1> models,Matter *& sample, PlasmaData const& pdata) : Model(sample,pdata,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
	TimeStep = 0;
}

ChargingModel::ChargingModel(std::string filename, double accuracy, std::array<bool,1> models,
				Matter *& sample, PlasmaGrid &pgrid) : Model(sample,pgrid,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename,double accuracy,std::array<bool,1> models,Matter *& sample, PlasmaGrid const& pgrid) : Model(sample,pgrid,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
	TimeStep = 0;
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

double ChargingModel::CheckTimeStep(){
	C_Debug( "\tIn ChargingModel::CheckTimeStep()\n\n" );
	// Deal with case where power/time step causes large temperature change.
	
	C_Debug("\nPdata.ElectronTemp = " << Pdata.ElectronTemp << "\nPdata.ElectronDensity = " << Pdata.ElectronDensity);
	
	double DebyeLength=sqrt((epsilon0*Kb*Pdata.ElectronTemp)/(Pdata.ElectronDensity*pow(echarge,2)));
	double PlasmaFreq = sqrt((Pdata.ElectronDensity*pow(echarge,2))/(epsilon0*Me));
	TimeStep = sqrt(2*PI) * ((DebyeLength)/Sample->get_radius()) 
			* (1/(PlasmaFreq*(1+Pdata.ElectronTemp/Pdata.IonTemp+Sample->get_potential())));

//	if(TimeStep != TimeStep)
//		TimeStep = 1e-8;	// (s), An estimate for regions of low plasma density
	C_Debug("\n\t\tDebyeLength = " << DebyeLength << "\n\t\tPlasmaFreq = " << PlasmaFreq 
			<< "\n\t\tTimeStep = " << TimeStep << "\n\n");

	assert(TimeStep == TimeStep);
	assert(TimeStep > 0);

	return TimeStep;
}

void ChargingModel::Charge(){
	C_Debug("\tIn ChargingModel::Charge()\n\n");

//	std::cout << "\n(4)Temp = " << Sample->get_temperature() << "\nVapPressure = " << Sample->get_vapourpressure() << "\nCv = " << Sample->get_heatcapacity() << "\n";

	// Assume the grain is negative and calculate potential
	if( UseModel[0] ){
		double Potential;
		if( Sample->get_deltatot() < 1.0 ){ // solveOML only defined for deltatot < 1.0
			Potential = solveOML(Sample->get_deltatot(),Sample->get_potential()); 
		}else{ // If the grain is in fact positive ...
			Potential = solveOML(Sample->get_deltatot(),Sample->get_potential());
			if( Potential < 0.0 ){
				Potential = solveOML(0.0,Sample->get_potential())-Kb*Sample->get_temperature()
						/(echarge*Pdata.ElectronTemp);
			}
		}
		Sample->update_charge(Potential,DeltaSec(),DeltaTherm());
//		std::cout << "\nPotential = " << Potential << "\nDeltaSec = " << Sample->get_deltasec() << "\nDeltatherm = " 
//			<< Sample->get_deltatherm() << "\n"; std::cin.get();

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
	double b = Pdata.IonTemp/Pdata.ElectronTemp;
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
	return sec(Pdata.ElectronTemp/1.16e5,Sample->get_elem());
}
