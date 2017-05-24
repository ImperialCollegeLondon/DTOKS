//#define PAUSE
//#define CHARGING_DEBUG

#include "ChargingModel.h"

// Default Constructor, no arguments
ChargingModel::ChargingModel(){
	CreateFile("Default_Charging_Filename.txt");
}

void ChargingModel::CreateFile(std::string filename){

	ChargingFile.open(filename);
	if( UseModel[0] ) ChargingFile << "Positive\tPotential";


	ChargingFile << "\n";
}

void ChargingModel::Print(){
	C_Debug("\n\nIn ChargingModel::Print(double TotPower, std::array<char,4> &ConstModels)");

	if( Sample->is_positive() )  ChargingFile << "Pos\t";
	if( !Sample->is_positive() ) ChargingFile << "Neg\t";
	if( UseModel[0] ) ChargingFile << Sample->get_potential();
	ChargingFile << "\n";
}

void ChargingModel::Charge(){
	// Assume the grain is negative and calculate potential
	double Potential = solveOML(Sample->get_deltatot(),Sample->get_potential());
	if( Sample->get_deltatot() >= 1.0 || Potential < 0.0 ){ // If the grain is in fact positive ...
		Potential = solveOML(0.0,Sample->get_potential())-Kb*Sample->get_temperature()/(echarge*Pdata.ElectronTemp);
	}
	Sample->update_charge(Potential);
	Print();
}

double ChargingModel::solveOML(double a, double guess){
	double b = Pdata.IonTemp/Pdata.ElectronTemp;
	double C = Me/Mi;
	
	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/((a-1)*exp(-guess) - sqrt(C/b) ) );

	while(fabs(guess-x1)>1e-2){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}
