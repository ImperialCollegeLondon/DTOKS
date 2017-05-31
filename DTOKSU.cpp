//#define PAUSE
#define DTOKSU_DEBUG
#include "DTOKSU.h"

// CONSIDER DEFINING A DEFAULT SAMPLE
std::array<bool, 9> DefaultHeatModels = {false,false,false,false,false,false, false,false,false};
std::array<bool,4> DefaultForceModels = {false,false,false,false};
std::array<bool,1> DefaultChargeModels = {false};
std::array<char,4> DefaultConstModels = { 'c','c','c','c'};
std::shared_ptr<Matter>	DefaultSample = std::make_shared<Tungsten>();

// Default Constructor, no arguments
DTOKSU::DTOKSU():
			CM("cf.txt",1.0,DefaultChargeModels,DefaultSample,PlasmaDefaults),
			HM("hf.txt",1e-9,1.0,DefaultHeatModels,DefaultSample,PlasmaDefaults),
			FM("ff.txt",1.0,DefaultForceModels,DefaultSample,PlasmaDefaults){
	D_Debug("\n\nIn DTOKSU::DTOKSU()\n\n");
	D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
	TimeStep = 0;
	TotalTime = 0;
	Pdata = PlasmaDefaults;
	Sample = DefaultSample;
}

DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, std::shared_ptr<Matter> const& sample, PlasmaData const &pdata,
				std::array<bool,9> &heatmodels, std::array<bool,4> &forcemodels, std::array<bool,1> &chargemodels)
				:CM("cf.txt",acclvls[0],chargemodels,sample,pdata),
				 HM("hf.txt",acclvls[1],timestep,heatmodels,sample,pdata),
				 FM("ff.txt",acclvls[2],forcemodels,sample,pdata){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
	TimeStep = 0;
	TotalTime = 0;
	Pdata = PlasmaDefaults;
	Sample = DefaultSample;
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

void DTOKSU::UpdatePData(){

}

int DTOKSU::Run(){
	D_Debug("- In DTOKSU::Run()\n\n");

	CM.CheckTimeStep();		// Check Time step length is appropriate
	CM.Charge();

	FM.CheckTimeStep();		// Check Time step length is appropriate
	FM.Force();

	for(int i = 0; i < 5; i ++){	
		HM.CheckTimeStep();		// Check Time step length is appropriate
		HM.Heat();
		std::cout << "\n(1)VapPressure = " << Sample->get_vapourpressure() << "Cv = " << Sample->get_heatcapacity();
//		Sample->update();               // CHECK THAT UPDATING HERE UPDATES EVERYWHERE!

//		std::cout << "\n(2)VapPressure = " << Sample->get_vapourpressure() << "Cv = " << Sample->get_heatcapacity();

	}

	return 0;
}


// ************************************* //

/*
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
*/

