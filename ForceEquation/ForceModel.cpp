//#define PAUSE
//#define FORCE_DEBUG

#include "ForceModel.h"

// Default Constructor, no arguments
ForceModel::ForceModel(){
}

void ForceModel::CreateFile(std::string filename){

	ForceFile.open(filename);
	ForceFile << "\n";
}

void ForceModel::Print(){
	F_Debug("\n\nIn ForceModel::Print(double TotPower, std::array<char,4> &ConstModels)");

	ForceFile << "\n";
}


