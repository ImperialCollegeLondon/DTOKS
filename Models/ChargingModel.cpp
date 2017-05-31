//#define PAUSE
#define CHARGING_DEBUG

#include "ChargingModel.h"

ChargingModel::ChargingModel():Model(){
	C_Debug("\n\nIn ChargingModel::ChargingModel():Model()\n\n");
	CreateFile("Default_Charging_Filename.txt");
	UseModel[0] = true;				// Charging Models turned on of possibly 9
	TimeStep = 0;
}

ChargingModel::ChargingModel(std::string filename, double accuracy, std::array<bool,1> models,
				std::shared_ptr <Matter> const& sample, PlasmaData const& pdata) : Model(sample,pdata,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename,std::array<bool,1> models,std::shared_ptr <Matter> const& sample, PlasmaData const& pdata) : Model(sample,pdata)\n\n");
	CreateFile(filename);
	UseModel = models;
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
/*
	// Take Eularian step to get initial time step
	C_Debug("\t"); double TotalPower = CalculatePower(Sample->get_temperature());
	// This model forces the time step to be the value which produces a change in temperature or 1*accuracy degree
	TimeStep = fabs((Sample->get_mass()*Sample->get_heatcapacity())/(TotalPower*accuracy));

	assert(TimeStep > 0);
*/
	TotalTime += TimeStep;
}

void ChargingModel::Charge(){
	C_Debug("\tIn ChargingModel::Charge()\n\n");
/*
	// Assume the grain is negative and calculate potential
	double Potential = solveOML(Sample->get_deltatot(),Sample->get_potential());
	if( Sample->get_deltatot() >= 1.0 || Potential < 0.0 ){ // If the grain is in fact positive ...
		Potential = solveOML(0.0,Sample->get_potential())-Kb*Sample->get_temperature()/(echarge*Pdata.ElectronTemp);
	}
	Sample->update_charge(Potential);
*/
	C_Debug("\t"); Print();
}

double ChargingModel::solveOML(double a, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
	double b = Pdata.IonTemp/Pdata.ElectronTemp;
	double C = Me/Mi;
	
	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/((a-1)*exp(-guess) - sqrt(C/b) ) );

	while(fabs(guess-x1)>1e-2){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}
