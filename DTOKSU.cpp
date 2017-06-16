//#define PAUSE
#define DTOKSU_DEBUG
#include "DTOKSU.h"

// CONSIDER DEFINING A DEFAULT SAMPLE
std::array<bool, 9> DefaultHeatModels = {false,false,false,false,false,false, false,false,false};
std::array<bool,4> DefaultForceModels = {false,false,false,false};
std::array<bool,1> DefaultChargeModels = {false};
std::array<char,4> DefaultConstModels = { 'c','c','c','c'};

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
	CreateFile("df.txt");
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
	CreateFile("df.txt");
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

	std::vector<double> Times(3,0);
	double HeatTime(0),ForceTime(0),ChargeTime(0);
	for(size_t i = 0; i < 100000; i ++){

		CM.Charge();
		Sample->update();			// Update data in GrainStructs

		ChargeTime = CM.CheckTimeStep();		// Check Time step length is appropriate
		ForceTime = FM.CheckTimeStep();			// Check Time step length is appropriate
		HeatTime = HM.CheckTimeStep();			// Check Time step length is appropriate
		if( HeatTime == 1) break;

		// We will assume Charging Time scale is much faster than either heating or moving.
		double MaxTimeStep = std::max(ForceTime,HeatTime);
		TimeStep = std::min(ForceTime,HeatTime);

		if( ChargeTime > TimeStep ){
			static bool runOnce = true;
			WarnOnce(runOnce,"*** WARNING! Charging Time scale is not the shortest timescale!! ***\n");
		}

		for(double IntermediateStep(0); IntermediateStep < MaxTimeStep; IntermediateStep += TimeStep){
			D_Debug( "\nIntermediateStep/MaxTimeStep = " << IntermediateStep << "/" << MaxTimeStep);
			double TestHeatTime = HM.CheckTimeStep();
			double TestForceTime = FM.CheckTimeStep();
			if( TimeStep == HeatTime){
				HM.Heat(); 
				TotalTime += TimeStep;
				if( TestForceTime < TimeStep || TestForceTime < HeatTime*0.1 ){
					std::cout << "\nWARNING! Smallest time scale changeover! Heat -> Force";
					TimeStep = TestForceTime;
				}
				D_Debug("\nHeat Step Taken.");
			}else if( TimeStep == ForceTime){
				FM.Force();
				TotalTime += TimeStep;
				if( TestHeatTime < TimeStep || TestHeatTime < HeatTime*0.1 ){
					std::cout << "\nWARNING! Smallest time scale changeover! Force -> Heat";
					TimeStep = TestHeatTime;
				}
				D_Debug("\nForce Step Taken.");
			}
			CM.Charge();
			Sample->update();
		}
		
		std::cout << "\ni : " << (int)i << "\nTemperature = " << Sample->get_temperature() << "\n\n"; std::cin.get();
		
		if( MaxTimeStep == HeatTime ){		HM.Heat(); D_Debug("\nHeat Step Taken."); }
		else if( MaxTimeStep == ForceTime ){	FM.Force(); D_Debug("\nForce Step Taken."); }	

		if( Sample->is_gas() && Sample->get_superboilingtemp() <= Sample->get_temperature() ){
			std::cout << "\n\nSample has Boiled ";
			break;
		}else if( Sample->is_gas() && Sample->get_superboilingtemp() > Sample->get_temperature() ){
			std::cout << "\n\nSample has Evaporated ";
			break;
		}

//		std::cout << "\n\tTimeStep = " << TimeStep << "\n\tChargeTime = " << CM.CheckTimeStep()
//			<< "\n\tForceTime = " << FM.CheckTimeStep() << "\n\tHeatTime = " << HM.CheckTimeStep() << "\n";
	}
	if( HeatTime == 1 ) std::cout << "\nThermal Equilibrium was achieved!\n\n";
	std::cout << "\nScript Complete.";
	return 0;
}
