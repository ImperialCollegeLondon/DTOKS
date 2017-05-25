//#define PAUSE
#define DTOKSU_DEBUG
#include "DTOKSU.h"

// Default Constructor, no arguments
DTOKSU::DTOKSU(){
	D_Debug("\n\nIn DTOKSU::DTOKSU()\n\n");
	D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
}

DTOKSU::DTOKSU( double timestep, std::shared_ptr<Matter> const& sample, PlasmaData const &pdata,
				std::array<bool,9> &heatmodels, std::array<bool,3> &forcemodels, std::array<bool,1> &chargemodels)
				:CM("cf.txt",chargemodels,sample,pdata),
				 HM("hf.txt",timestep,heatmodels,sample,pdata),
				 FM("ff.txt",forcemodels,sample,pdata){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
}

void DTOKSU::CreateFile( std::string filename ){

	D_Debug("\n\nIn DTOKSU::CreateFile(std::string filename)\n\n");
	MyFile.open(filename);
	MyFile << "\n";
}

void DTOKSU::CheckTimeStep(){
	D_Debug( "\tIn DTOKSU::CheckTimeStep()\n\n");
	assert(TimeStep > 0);
	D_Debug( "\nTimeStep = " << TimeStep );
	TotalTime += TimeStep;
}

void DTOKSU::Print(){
	D_Debug("\tIn DTOKSU::Print()\n\n");

	MyFile 	<< "\nSomeStuff";

	MyFile << "\n";
}

int DTOKSU::Run(){
	D_Debug("- In DTOKSU::Run()\n\n");
	HM.Heat('c');
	CM.Charge();
	FM.Force();
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


