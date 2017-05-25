//#define PAUSE
#define DTOKSU_DEBUG
#include "DTOKSU.h"

// Default Constructor, no arguments
DTOKSU::DTOKSU(){
	D_Debug("\n\nIn DTOKSU::DTOKSU()");
}


void DTOKSU::CreateFile(std::string filename){
	D_Debug("\n\nIn DTOKSU::CreateFile(std::string filename)");
	MyFile.open(filename);
	MyFile << "Yo";
	MyFile << "\n";
}

void DTOKSU::CheckTimeStep(){
	D_Debug( "\n\nIn DTOKSU::CheckTimeStep()" );
	assert(TimeStep > 0);
	D_Debug( "\nTimeStep = " << TimeStep );
	TotalTime += TimeStep;
}

void DTOKSU::Print(){
	D_Debug("\n\nIn DTOKSU::Print()");

	MyFile 	<< "\nSomeStuff";

	MyFile << "\n";
}

int DTOKSU::Run(){
	D_Debug("\n\nIn DTOKSU::Run()");
	return 0;
}


// ************************************* \\

const double DTOKSU::DeltaTherm(double DustTemperature)const{
//	D_Debug("\n\nIn HeatingModel::DeltaTherm():");
	return (Richardson*pow(DustTemperature,2)*exp(-(Sample->get_workfunction()*echarge)/(Kb*DustTemperature)))
		/echarge;
}

const double DTOKSU::DeltaSec()const{
//	D_Debug("\n\nIn HeatingModel::DeltaSec():");
	return sec(Pdata.ElectronTemp/1.16e5,Sample->get_elem()); 	// Convert from K to ev
}

const double DTOKSU::DeltaTot(double DustTemperature)const{
//	D_Debug("\n\nIn HeatingModel::DeltaTot():");
	return (DeltaSec() + DeltaTherm(DustTemperature));
}


