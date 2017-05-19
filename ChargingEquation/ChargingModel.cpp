//#define PAUSE
//#define CHARGING_DEBUG

#include "ChargingModel.h"

// Default Constructor, no arguments
ChargingModel::ChargingModel(){
}

void ChargingModel::CreateFile(std::string filename){

	ChargingFile.open(filename);
	ChargingFile << "\n";
}

void ChargingModel::Print(){
	F_Debug("\n\nIn ChargingModel::Print(double TotPower, std::array<char,4> &ConstModels)");

	ChargingFile << "\n";
}


