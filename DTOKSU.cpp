//#define PAUSE
//#define DTOKSU_DEBUG
//#define DTOKSU_DEEP_DEBUG
#include "DTOKSU.h"

// CONSIDER DEFINING A DEFAULT SAMPLE
std::array<bool, 9> DefaultHeatModels = {false,false,false,false,false,false, false,false,false};
std::array<bool,5> DefaultForceModels = {false,false,false,false,false};
std::array<bool,1> DefaultChargeModels = {false};
std::array<char,4> DefaultConstModels = { 'c','c','c','c'};

DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, PlasmaData *&pdata,
				std::array<bool,9> &heatmodels, std::array<bool,5> &forcemodels, std::array<bool,1> &chargemodels)
			: Sample(sample),
				CM("Data/cf.txt",acclvls[0],chargemodels,sample,pdata),
				HM("Data/hf.txt",acclvls[1],heatmodels,sample,pdata),
				FM("Data/ff.txt",acclvls[2],forcemodels,sample,pdata){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
//	std::cout << "\nacclvls[0] = " << acclvls[0];
	MinTimeStep = timestep;
	TotalTime = 0;
	CreateFile("Data/df.txt");
}

DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, PlasmaGrid &pgrid,
				std::array<bool,9> &heatmodels, std::array<bool,5> &forcemodels, std::array<bool,1> &chargemodels)
				: Sample(sample),
				CM("Data/cf.txt",acclvls[0],chargemodels,sample,pgrid),
				HM("Data/hf.txt",acclvls[1],heatmodels,sample,pgrid),
				FM("Data/ff.txt",acclvls[2],forcemodels,sample,pgrid){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");

//	D_Debug("\nHeatModels = " << heatmodels[0]);
//	std::cout << "\nacclvls[0] = " << acclvls[0];
	MinTimeStep = timestep;
	TotalTime = 0;
	CreateFile("Data/df.txt");
}

void DTOKSU::CreateFile( std::string filename ){
	D_Debug("\n\nIn DTOKSU::CreateFile(std::string filename)\n\n");
	MyFile.open(filename);
	MyFile << "TotalTime\n";
}

void DTOKSU::OpenFiles( std::string filename, unsigned int i ){
	D_Debug("\n\nIn DTOKSU::OpenFiles(std::string filename, unsigned int i)\n\n");
	CreateFile(filename + "_df_" + std::to_string(i) + ".txt");
	HM.CreateFile(filename + "_hm_" + std::to_string(i) + ".txt",false);
	FM.CreateFile(filename + "_fm_" + std::to_string(i) + ".txt");
	CM.CreateFile(filename + "_cm_" + std::to_string(i) + ".txt");
}

void DTOKSU::CloseFiles(){
	D_Debug("\n\nIn DTOKSU::CloseFiles()\n\n");
	MyFile.close();
	HM.CloseFile();
	FM.CloseFile();
	CM.CloseFile();
}

void DTOKSU::CheckTimeStep(){
	D_Debug( "\tIn DTOKSU::CheckTimeStep()\n\n");
	
	assert(MinTimeStep > 0);

	D_Debug( "\nMinTimeStep = " << MinTimeStep );
	TotalTime += MinTimeStep;
}

void DTOKSU::Print(){
	D_Debug("\tIn DTOKSU::Print()\n\n");

	MyFile 	<< TotalTime;

	MyFile << "\n";
}

int DTOKSU::Run(){
	D_Debug("- In DTOKSU::Run()\n\n");

	double HeatTime(0),ForceTime(0),ChargeTime(0);

	bool InGrid = FM.update_plasmadata(Sample->get_position());
	CM.Charge();
	Sample->update();			// Need to manually update the first time as first step is not necessarily heating
	while( InGrid && !Sample->is_split() ){

		// ***** START OF : DETERMINE TIMESCALES OF PROCESSES ***** //	

		ChargeTime 	= CM.UpdateTimeStep();		// Check Time step length is appropriate
		ForceTime 	= FM.UpdateTimeStep();		// Check Time step length is appropriate
		HeatTime 	= HM.UpdateTimeStep();		// Check Time step length is appropriate
		if( HeatTime == 1) break;			// Thermal Equilibrium Reached

		// We will assume Charging Time scale is much faster than either heating or moving, but check for the other case.
		double MaxTimeStep = std::max(ForceTime,HeatTime);
		MinTimeStep = std::min(ForceTime,HeatTime);

		// Check Charging timescale isn't the fastest timescale.
		if( ChargeTime > MinTimeStep && ChargeTime != 1){
			static bool runOnce = true;
			WarnOnce(runOnce,"*** Charging Time scale is not the shortest timescale!! ***\n");
		}

		// ***** END OF : DETERMINE TIMESCALES OF PROCESSES ***** //	
		// ***** START OF : NUMERICAL METHOD BASED ON TIME SCALES ***** //	
	
		// Resolve region where the Total Power is zero for a Plasma Grid.
		// This typically occurs when plasma parameters are zero in a cell and other models are off (or zero)... 
		// Even in No Plasma Region, cooling processes can still occur, but this is specifically for zero power
		
		if( HeatTime == 10 ){
			D1_Debug("\nNo Net Power Region...");
			FM.Force();
			HM.AddTime(ForceTime);
			CM.Charge(ForceTime);
			TotalTime += ForceTime;
		}else if( Sample->get_deltatot() > 0.5 && Sample->get_deltatot() < 1.0 ){ // Region of rapid charge variation
			// WARNING, CURRENT TESTING SHOWS DTOKSU DOESN'T ENTER HERE AT ALL WITH PGRID.
			// I'm pretty sure that the grain does become positive at some point so this is likely because
			// it is going through that region in the for loop below...
			D1_Debug("\n\nPotential Focus Region, steps taken at 0.01*MinTimeStep\n");
			D1_Debug("Potential = " << Sample->get_potential() << "\nDeltaTot = " << Sample->get_deltatot() << "\n\n");
			HM.Heat  (ChargeTime);
			FM.Force (ChargeTime);
			CM.Charge(ChargeTime);
			TotalTime += ChargeTime;
		// Else If the timescales of the processes are comparable, step through each at the faster timescale
		}else if( MinTimeStep*2.0 > MaxTimeStep){
			D1_Debug("\nComparable Timescales, taking time steps through both processes at shorter time scale");
			HM.Heat(MinTimeStep);
			FM.Force(MinTimeStep);
			CM.Charge(MinTimeStep);
			TotalTime += MinTimeStep;
		}else{ // Else, we can take steps through the smaller one til the sum of the steps is the larger.
			D1_Debug("\nDifferent Timescales, taking many time steps through quicker process at shorter time scale");
			unsigned int j(1);
			for( j =1; (j*MinTimeStep) < MaxTimeStep; j ++){
				D1_Debug( "\nIntermediateStep/MaxTimeStep = " << j*MinTimeStep << "/" << MaxTimeStep);

				// Take the time step in the faster time process
				if( MinTimeStep == HeatTime ){
					HM.Heat(MinTimeStep);
					D1_Debug("\nHeat Step Taken.");
					// Check that time scales haven't changed significantly whilst looping... 
				}else if( MinTimeStep == ForceTime ){
					FM.Force(MinTimeStep);	
					D1_Debug("\nForce Step Taken.");
					// Check that time scales haven't changed significantly whilst looping... 
				}else{
					std::cerr << "\nUnexpected Timescale Behaviour (1)!";
				}
				CM.Charge(MinTimeStep);


				// Check that time scales haven't changed significantly whilst looping...
				if( ForceTime/FM.ProbeTimeStep() > 2 ){
					D1_Debug("\nForce TimeStep Has Changed Significantly whilst taking small steps...");
					j ++;
					break; // Can't do this: MaxTimeStep = j*MinTimeStep; as we change MaxTimeStep...
				}
				if( HeatTime/HM.ProbeTimeStep() > 2 ){
					D1_Debug("\nHeat TimeStep Has Changed Significantly whilst taking small steps...");
					j ++;
					break; // Can't do this: MaxTimeStep = j*MinTimeStep; as we change MaxTimeStep...
				}
			}

			// Take a time step in the slower time process
			D1_Debug("\n*STEP* = " << (j-1)*MinTimeStep << "\nMaxTimeStep = " << MaxTimeStep << "\nj = " << j << "\n");
				
			if( MaxTimeStep == HeatTime ){
				HM.Heat((j-1)*MinTimeStep); 
				D_Debug("\nHeat Step Taken."); 
			}else if( MaxTimeStep == ForceTime ){	
				FM.Force((j-1)*MinTimeStep); 
				D_Debug("\nForce Step Taken."); 
			}else{	std::cerr << "\nUnexpected Timescale Behaviour! (2)";	}
			TotalTime += (j-1)*MinTimeStep;
			// This is effectively 'an extra step' which is necessary because we need the last step to be charging
			// So we set the time step to be arbitrarily small... Not the best practice ever but okay.
			CM.Charge(1e-100);	
		}
		// ***** END OF : NUMERICAL METHOD BASED ON TIME SCALES ***** //	
		D_Debug("\nTemperature = " << Sample->get_temperature() << "\n\n"); 
		D_Debug("\n\tMinTimeStep = " << MinTimeStep << "\n\tChargeTime = " << ChargeTime
			<< "\n\tForceTime = " << ForceTime << "\n\tHeatTime = " << HeatTime << "\n");

		// Update the plasma data from the plasma grid for all models...
		InGrid = FM.update_plasmadata(Sample->get_position());
		Print();
		// ***** START OF : DETERMINE IF END CONDITION HAS BEEN REACHED ***** //
		if( Sample->is_gas() && Sample->get_superboilingtemp() <= Sample->get_temperature() ){
			std::cout << "\n\nSample has Boiled ";
			break;
		}else if( Sample->is_gas() && Sample->get_superboilingtemp() > Sample->get_temperature() ){
			std::cout << "\n\nSample has Evaporated ";
			break;
		}else if( Sample->is_split() && Sample->is_gas() ){
			std::cout << "\n\nSample has broken up into a gas electrostatically";
			break; // Could replace with return 3; and leave off the end check?...
		}else if( Sample->is_split() ){
			std::cout << "\n\nSample has broken up into two large parts";
			break; // Could replace with return 3; and leave off the end check?...
		}else if( Sample->is_gas() ){
			std::cout << "\n\nSample has vapourised";
			break;
		
		}
		// ***** END OF : DETERMINE IF END CONDITION HAS BEEN REACHED ***** //
	}

	if( fabs(1 - (HM.get_totaltime()/FM.get_totaltime())) > 0.001  
		|| fabs(1 - (FM.get_totaltime()/CM.get_totaltime())) > 0.001 ){
		std::cout << "\nWarning! Total Times recorded by processes don't match!";
		std::cout << "\nHM.get_totaltime() = " << HM.get_totaltime() << "\nFM.get_totaltime() = " 
			<< FM.get_totaltime() << "\nCM.get_totaltime() = " << CM.get_totaltime() << "\n\nTotalTime = " << TotalTime;

	}
	if( !InGrid ){
		std::cout << "\nSample has left simulation domain";
		return 1;
	}else if( HeatTime == 1 ){
		std::cout << "\nEnd of Run, exiting due to Thermal Equilibrium being reached!\n\n";
		return 2;
	}else if( Sample->is_split() && !Sample->is_gas() ){
		return 3;	// Return status for continue simulating Breakup condition
	}

	std::cout << "\nFinished DTOKS-U run.";
	return 0;
}

