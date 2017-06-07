//#define PAUSE
#define DTOKSU_DEBUG
#include "DTOKSU.h"

// CONSIDER DEFINING A DEFAULT SAMPLE
std::array<bool, 9> DefaultHeatModels = {false,false,false,false,false,false, false,false,false};
std::array<bool,4> DefaultForceModels = {false,false,false,false};
std::array<bool,1> DefaultChargeModels = {false};
std::array<char,4> DefaultConstModels = { 'c','c','c','c'};

/*
struct PlasmaData PlasmaDefaults1 = {
	1e20,		// m^-3, Neutral Density
	1e20,		// m^-3, Electron Density
	1e20,		// m^-3, Electron Density
	116045.25,	// K, Ion Temperature
	116045.25,	// K, Electron Temperature
	116045.25,	// K, Neutral Temperature
	300,		// K, Ambient Temperature
	threevector(),	// m s^-1, Plasma Velocity (Should eventually be normalised to sound speed cs)
	threevector(),	// m s^-2, Acceleration due to gravity
	threevector(),	// V m^-1, Electric field at dust location (Normalised later) 
	threevector(),	// T, Magnetic field at dust location (Normalised later)
};

// Default Constructor, no arguments. This is not a well defined constructor so has been commented out

Matter *DefaultSample = new Tungsten();
DTOKSU::DTOKSU():
			Pgrid('h','m',0.01),	// Default configuration for MAST
			CM("cf.txt",1.0,DefaultChargeModels,DefaultSample,PlasmaDefaults1),
			HM("hf.txt",1e-9,1.0,DefaultHeatModels,DefaultSample,PlasmaDefaults1),
			FM("ff.txt",1.0,DefaultForceModels,DefaultSample,PlasmaDefaults1){
	D_Debug("\n\nIn DTOKSU::DTOKSU()\n\n");
	D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
	TimeStep = 0;
	TotalTime = 0;

}
*/
DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, PlasmaData &pdata,
				std::array<bool,9> &heatmodels, std::array<bool,4> &forcemodels, std::array<bool,1> &chargemodels)
				: Sample(sample),
				CM("cf.txt",acclvls[0],chargemodels,sample,pdata),
				HM("hf.txt",acclvls[1],heatmodels,sample,pdata),
				FM("ff.txt",acclvls[2],forcemodels,sample,pdata){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
//	std::cout << "\nacclvls[0] = " << acclvls[0];
	TimeStep = timestep;
	TotalTime = 0;
}

DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, PlasmaGrid &pgrid,
				std::array<bool,9> &heatmodels, std::array<bool,4> &forcemodels, std::array<bool,1> &chargemodels)
				: Sample(sample),
				CM("cf.txt",acclvls[0],chargemodels,sample,pgrid),
				HM("hf.txt",acclvls[1],heatmodels,sample,pgrid),
				FM("ff.txt",acclvls[2],forcemodels,sample,pgrid){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");

//	D_Debug("\nHeatModels = " << heatmodels[0]);
//	std::cout << "\nacclvls[0] = " << acclvls[0];
	TimeStep = timestep;
	TotalTime = 0;
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

	for(size_t i = 0; i < 10000000; i ++){

		double ChargeTime = CM.CheckTimeStep();		// Check Time step length is appropriate
		double ForceTime = FM.CheckTimeStep();		// Check Time step length is appropriate
		double HeatTime = HM.CheckTimeStep();		// Check Time step length is appropriate

		TimeStep = std::max(ChargeTime,std::max(ForceTime,HeatTime));

//		CM.Charge();
//		FM.Force();
		HM.Heat();

//		HM.update();

//		Pgrid.readdata();	// Read data

//		std::cout << "\n\tTimeStep = " << TimeStep << "\n\tChargeTime = " << ChargeTime 
//			<< "\n\tForceTime = " << ForceTime << "\n\tHeatTime = " <<HeatTime << "\n";

		Sample->update();               // CHECK THAT UPDATING HERE UPDATES EVERYWHERE!

//		std::cout << "\n(1)Temperature = " << Sample->get_temperature();


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

