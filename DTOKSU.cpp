//#define PAUSE
#define DTOKSU_DEBUG
//#define DTOKSU_DEEP_DEBUG
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
	MinTimeStep = timestep;
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
	MinTimeStep = timestep;
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
	
	assert(MinTimeStep > 0);

	D_Debug( "\nMinTimeStep = " << MinTimeStep );
	TotalTime += MinTimeStep;
}

void DTOKSU::Print(){
	D_Debug("\tIn DTOKSU::Print()\n\n");

	MyFile 	<< "\nSomeStuff";

	MyFile << "\n";
}

int DTOKSU::Run(){
	D_Debug("- In DTOKSU::Run()\n\n");

	double HeatTime(0),ForceTime(0),ChargeTime(0);
	for(size_t i = 0; i < 100000; i ++){

		CM.Charge();
		Sample->update();			// Update data in GrainStructs

		ChargeTime 	= CM.UpdateTimeStep();		// Check Time step length is appropriate
		ForceTime 	= FM.UpdateTimeStep();			// Check Time step length is appropriate
		HeatTime 	= HM.UpdateTimeStep();			// Check Time step length is appropriate
		if( HeatTime == 1) break;			// Thermal Equilibrium Reached

		// We will assume Charging Time scale is much faster than either heating or moving, but check for the other case.
		double MaxTimeStep = std::max(ForceTime,HeatTime);
		MinTimeStep = std::min(ForceTime,HeatTime);

		// Check Charging timescale isn't the fastest timescale.
		if( ChargeTime > MinTimeStep ){
			static bool runOnce = true;
			WarnOnce(runOnce,"*** WARNING! Charging Time scale is not the shortest timescale!! ***\n");
		}


		// Resolve region of rapid charge variation
		if( Sample->get_deltatot() > 0.5 && Sample->get_deltatot() < 1.0 ){
			D1_Debug("\n\nPotential Focus Region, steps taken at 0.01*MinTimeStep\n");
			D1_Debug("Potential = " << Sample->get_potential() << "\nDeltaTot = " << Sample->get_deltatot() << "\n\n");
			CM.Charge(MinTimeStep*0.01);
			FM.Force(MinTimeStep*0.01);
			HM.Heat(MinTimeStep*0.01);
		// Else If the timescales of the processes are comparable, step through each at the faster timescale
		}else if( MinTimeStep*2.0 > MaxTimeStep){
			D1_Debug("\nComparable TimeStep lengths, taking time steps at shorter time scale");
			HM.Heat(MinTimeStep);
			FM.Force(MinTimeStep);
		}else{ // Else, we can take a steps through the smaller one til the sum of the steps is the larger.
			unsigned int j(1);
			for( j =1; (j*MinTimeStep) < MaxTimeStep; j ++){
				D1_Debug( "\nIntermediateStep/MaxTimeStep = " << j*MinTimeStep << "/" << MaxTimeStep);

				// Take the time step in the faster time process
				if( MinTimeStep == HeatTime ){
					HM.Heat(MinTimeStep);
					D1_Debug("\nHeat Step Taken.");
				}else if( MinTimeStep == ForceTime ){
					FM.Force(MinTimeStep);	
					D1_Debug("\nForce Step Taken.");
				}else{
					std::cerr << "\nUnexpected Timescale Behaviour (1)!";
				}
				CM.Charge();
				Sample->update();
			}
			// Take a time step in the slower time process
			D1_Debug("\nStep = " << (j-1)*MinTimeStep << "\nMaxTimeStep = " << MaxTimeStep << "\nj = " << j << "\n");
			if( MaxTimeStep == HeatTime ){		HM.Heat((j-1)*MinTimeStep); D_Debug("\nHeat Step Taken."); }
			else if( MaxTimeStep == ForceTime ){	FM.Force((j-1)*MinTimeStep); D_Debug("\nForce Step Taken."); }
			else{					std::cerr << "\nUnexpected Timescale Behaviour! (2)";	}
		}
		D1_Debug("\ni : " << (int)i << "\nTemperature = " << Sample->get_temperature() << "\n\n"); 
//		if( i > 240 && i < 300 ){ std::cin.get(); }
		D_Debug("\n\tMinTimeStep = " << MinTimeStep << "\n\tChargeTime = " << ChargeTime
			<< "\n\tForceTime = " << ForceTime << "\n\tHeatTime = " << HeatTime << "\n");

		if( Sample->is_gas() && Sample->get_superboilingtemp() <= Sample->get_temperature() ){
			std::cout << "\n\nSample has Boiled ";
			break;
		}else if( Sample->is_gas() && Sample->get_superboilingtemp() > Sample->get_temperature() ){
			std::cout << "\n\nSample has Evaporated ";
			break;
		}

		TotalTime += MinTimeStep;
	}
	if( HeatTime == 1 ) std::cout << "\nThermal Equilibrium was achieved!\n\n";
	std::cout << "\nScript Complete.";
	return 0;
}

//				if( TestForceTime < MinTimeStep || TestForceTime < HeatTime*0.1 ){
				// Check if the timescale hasn't changed substantially this step... (May not be necessary)
//				double TestHeatTime = HM.CheckTimeStep();
//				double TestForceTime = FM.CheckTimeStep();
//				if( TestForceTime < MinTimeStep || TestForceTime < HeatTime*0.5 ){
//					std::cout << "\nWARNING! Smallest time scale changeover! Heat -> Force";
//					MinTimeStep = TestForceTime;
//				}
//				if( TestHeatTime < MinTimeStep || TestHeatTime < ForceTime*0.5 ){
//					std::cout << "\nWARNING! Smallest time scale changeover! Force -> Heat";
//					MinTimeStep = TestHeatTime;
//				}
//				if( MinTimeStep == HeatTime || MinTimeStep == TestHeatTime){
