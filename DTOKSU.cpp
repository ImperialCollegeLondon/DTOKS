//#define PAUSE
#define DTOKSU_DEBUG
#include "DTOKSU.h"

// CONSIDER DEFINING A DEFAULT SAMPLE
std::array<bool, 9> DefaultHeatModels = {false,false,false,false,false,false, false,false,false};
std::array<bool,4> DefaultForceModels = {false,false,false,false};
std::array<bool,1> DefaultChargeModels = {false};
std::array<char,4> DefaultConstModels = { 'c','c','c','c'};


// Default Constructor, no arguments. This is not a well defined constructor so has been commented out
/*
Matter *DefaultSample = new Tungsten();
DTOKSU::DTOKSU():
			Pgrid('h','m',0.01),	// Default configuration for MAST
			CM("cf.txt",1.0,DefaultChargeModels,DefaultSample,PlasmaDefaults),
			HM("hf.txt",1e-9,1.0,DefaultHeatModels,DefaultSample,PlasmaDefaults),
			FM("ff.txt",1.0,DefaultForceModels,DefaultSample,PlasmaDefaults){
	D_Debug("\n\nIn DTOKSU::DTOKSU()\n\n");
	D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
	TimeStep = 0;
	TotalTime = 0;

}
*/
DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, PlasmaData const &pdata,
				std::array<bool,9> &heatmodels, std::array<bool,4> &forcemodels, std::array<bool,1> &chargemodels)
				: Pgrid('h','m',0.01),
				CM("cf.txt",acclvls[0],chargemodels,sample,pdata),
				HM("hf.txt",timestep,acclvls[1],heatmodels,sample,pdata),
				FM("ff.txt",acclvls[2],forcemodels,sample,pdata){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
//	std::cout << "\nacclvls[0] = " << acclvls[0];
	TimeStep = 0;
	TotalTime = 0;
	int p(0), k(0);
	Pgrid.locate(p,k,sample->get_position()); // This will cause an assertion to be raised if default position is (0,0,0)
	Pgrid.setfields(p,k);
}

DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, char plasma, char machine, double spacing,
				std::array<bool,9> &heatmodels, std::array<bool,4> &forcemodels, std::array<bool,1> &chargemodels)
				: Pgrid(plasma,machine,spacing),
				CM("cf.txt",acclvls[0],chargemodels,sample,Pgrid.get_plasmadata(sample->get_position())),
				HM("hf.txt",timestep,acclvls[1],heatmodels,sample,Pgrid.get_plasmadata(sample->get_position())),
				FM("ff.txt",acclvls[2],forcemodels,sample,Pgrid.get_plasmadata(sample->get_position())){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
//	std::cout << "\nacclvls[0] = " << acclvls[0];
	TimeStep = 0;
	TotalTime = 0;
	int p(0), k(0);
	Pgrid.locate(p,k,sample->get_position());
	Pgrid.setfields(p,k);
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

	int i(0), k(0);
	for(int i = 0; i < 5; i ++){
		
		// if time step is larger than this, it will cross a cell in one time step
//		double dt = 0.5*Pgrid.getdl()/((dust.getvd()).mag3()); 
		double ChargeTime = CM.CheckTimeStep();		// Check Time step length is appropriate
		double ForceTime = FM.CheckTimeStep();		// Check Time step length is appropriate
		double HeatTime = HM.CheckTimeStep();		// Check Time step length is appropriate

		CM.Charge();
		FM.Force();
		HM.Heat();

//		HM.update();

//		Pgrid.readdata();	// Read data
//		Pgrid.locate(i, k, Sample.get_position() ); // Locate the grid point nearest the dust
//		Pgrid.setfields(i,k); // Calculate fields

		std::cout << "\nChargeTime = " << ChargeTime << "\nForceTime = " << ForceTime << "\nHeatTime = " <<HeatTime;
//		std::cout << "\n(1)Temp = " << Sample->get_temperature() << "\nVapPressure = " << Sample->get_vapourpressure() << "\nCv = " << Sample->get_heatcapacity() << "\n";
		//Sample->update();               // CHECK THAT UPDATING HERE UPDATES EVERYWHERE!

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

