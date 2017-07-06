//#define PAUSE
//#define DTOKSU_DEBUG
//#define DTOKSU_DEEP_DEBUG
#include "DTOKSU.h"

// CONSIDER DEFINING A DEFAULT SAMPLE
std::array<bool,9> DefaultHeatModels   = {false,false,false,false,false,false, false,false,false};
std::array<bool,5> DefaultForceModels  = {false,false,false,false,false};
std::array<bool,1> DefaultChargeModels = {false};
std::array<char,4> DefaultConstModels  = { 'c','c','c','c'};

DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, PlasmaData *&pdata,
				std::array<bool,9> &heatmodels, std::array<bool,5> &forcemodels, std::array<bool,1> &chargemodels){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
	Sample.push_back(sample);
	HM.push_back(HeatingModel("hf.txt",acclvls[1],heatmodels,sample,pdata));
	CM.push_back(ChargingModel("cf.txt",acclvls[0],chargemodels,sample,pdata));
	FM.push_back(ForceModel("ff.txt",acclvls[2],forcemodels,sample,pdata));
	MinTimeStep = timestep;
	TotalTime = 0;
	CreateFile("df.txt");
	D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
}

DTOKSU::DTOKSU( double timestep, std::array<double,3> acclvls, Matter *& sample, PlasmaGrid &pgrid,
				std::array<bool,9> &heatmodels, std::array<bool,5> &forcemodels, std::array<bool,1> &chargemodels){
        D_Debug("\n\nIn DTOKSU::DTOKSU( ... )\n\n");
	Sample.push_back(sample);
	HM.push_back(HeatingModel("hf.txt",acclvls[1],heatmodels,sample,pgrid));
	CM.push_back(ChargingModel("cf.txt",acclvls[0],chargemodels,sample,pgrid));
	FM.push_back(ForceModel("ff.txt",acclvls[2],forcemodels,sample,pgrid));
	MinTimeStep = timestep;
	TotalTime = 0;
	CreateFile("df.txt");
        D_Debug("\n\n************************************* SETUP FINISHED ************************************* \n\n");
}

void DTOKSU::CreateFile( std::string filename ){
	D_Debug("\n\nIn DTOKSU::CreateFile(std::string filename)\n\n");
	MyFile.open(filename);
	MyFile << "TotalTime\n";
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

//	std::cout << "\nHM.size() = " << HM.size(); std::cin.get();
	if( HM.size() == 0 || FM.size() == 0 || CM.size() == 0 ){
		std::cout << "\nSomething's gone terribly wrong...";
		return -1000;
	}

	HM[0].update_sample(Sample[0]);
	HM[0].update_plasmadata(Sample[0]->get_position());
	HM[0].CreateFile("hf.txt",false);
	PlasmaData *pdatatemp = HM[0].get_plasmadata();
	CM[0].update_sample(Sample[0]);
	CM[0].update_plasmadata(pdatatemp);
	CM[0].CreateFile("cf.txt");
	FM[0].update_sample(Sample[0]);
	FM[0].update_plasmadata(pdatatemp);
	FM[0].CreateFile("ff.txt");

	for(unsigned int i(0); i < HM.size(); i ++){
	
		CM[i].Charge();
		Sample[i]->update();		// Need to manually update the first time as first step is not necessarily heating
		bool InGrid = FM[i].update_plasmadata(Sample[i]->get_position());
		while( InGrid ){
	
			// ***** START OF : DETERMINE TIMESCALES OF PROCESSES ***** //	
			std::cout << "\ni = " << i;
			std::cout << "\nSample[i]->get_position() = " << Sample[i]->get_position();
			ChargeTime 	= CM[i].UpdateTimeStep();		// Check Time step length is appropriate
			ForceTime 	= FM[i].UpdateTimeStep();		// Check Time step length is appropriate
			HeatTime 	= HM[i].UpdateTimeStep();		// Check Time step length is appropriate
			if( HeatTime == 1) break;			// Thermal Equilibrium Reached
	
			// We will assume Charging Time scale is much faster than either heating or moving, 
			// but check for the other case.
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
				FM[i].Force();
				HM[i].AddTime(ForceTime);
				CM[i].Charge(ForceTime);
				TotalTime += ForceTime;
			}else if( Sample[i]->get_deltatot() > 0.5 && Sample[i]->get_deltatot() < 1.0 ){ // For rapid charge variation
				// WARNING, CURRENT TESTING SHOWS DTOKSU DOESN'T ENTER HERE AT ALL WITH PGRID.
				// I'm pretty sure that the grain does become positive at some point so this is likely because
				// it is going through that region in the for loop below...
				D1_Debug("\n\nPotential Focus Region, steps taken at 0.01*MinTimeStep\n");
				D1_Debug("Potential = " << Sample[i]->get_potential());
				D1_Debug("\nDeltaTot = " << Sample[i]->get_deltatot() << "\n\n");
				HM[i].Heat  (ChargeTime);
				FM[i].Force (ChargeTime);
				CM[i].Charge(ChargeTime);
				TotalTime += ChargeTime;
			// Else If the timescales of the processes are comparable, step through each at the faster timescale
			}else if( MinTimeStep*2.0 > MaxTimeStep){
				D1_Debug("\nComparable Timescales, taking time steps through both processes at shorter time scale");
				HM[i].Heat(MinTimeStep);
				FM[i].Force(MinTimeStep);
				CM[i].Charge(MinTimeStep);
				TotalTime += MinTimeStep;
			}else{ // Else, we can take steps through the smaller one til the sum of the steps is the larger.
				D1_Debug("\nDifferent Timescales,");
				D1_Debug("taking many time steps through quicker process at shorter time scale");
				unsigned int j(1);
				for( j =1; (j*MinTimeStep) < MaxTimeStep; j ++){
					D1_Debug( "\nIntermediateStep/MaxTimeStep = " << j*MinTimeStep << "/" << MaxTimeStep);
	
					// Take the time step in the faster time process
					if( MinTimeStep == HeatTime ){
						HM[i].Heat(MinTimeStep);
						D1_Debug("\nHeat Step Taken.");
						// Check that time scales haven't changed significantly whilst looping... 
					}else if( MinTimeStep == ForceTime ){
						FM[i].Force(MinTimeStep);	
						D1_Debug("\nForce Step Taken.");
						// Check that time scales haven't changed significantly whilst looping... 
					}else{
						std::cerr << "\nUnexpected Timescale Behaviour (1)!";
					}
					CM[i].Charge(MinTimeStep);
	
	
					// Check that time scales haven't changed significantly whilst looping...
					if( ForceTime/FM[i].ProbeTimeStep() > 2 ){
						D1_Debug("\nForce TimeStep Has Changed Significantly whilst taking small steps...");
						j ++;
						break; // Can't do this: MaxTimeStep = j*MinTimeStep; as we change MaxTimeStep...
					}
					if( HeatTime/HM[i].ProbeTimeStep() > 2 ){
						D1_Debug("\nHeat TimeStep Has Changed Significantly whilst taking small steps...");
						j ++;
						break; // Can't do this: MaxTimeStep = j*MinTimeStep; as we change MaxTimeStep...
					}
				}
	
				// Take a time step in the slower time process
				D1_Debug("\n*STEP* = " << (j-1)*MinTimeStep << "\nMaxTimeStep = ")
				D1_Debug( MaxTimeStep << "\nj = " << j << "\n");
					
				if( MaxTimeStep == HeatTime ){
					HM[i].Heat((j-1)*MinTimeStep); 
					D_Debug("\nHeat Step Taken."); 
				}else if( MaxTimeStep == ForceTime ){	
					FM[i].Force((j-1)*MinTimeStep); 
					D_Debug("\nForce Step Taken."); 
				}else{	std::cerr << "\nUnexpected Timescale Behaviour! (2)";	}
				TotalTime += (j-1)*MinTimeStep;
				CM[i].Charge(MinTimeStep);
			}
			// ***** END OF : NUMERICAL METHOD BASED ON TIME SCALES ***** //	
			D_Debug("\nTemperature = " << Sample[i]->get_temperature() << "\n\n"); 
			D_Debug("\n\tMinTimeStep = " << MinTimeStep << "\n\tChargeTime = " << ChargeTime
				<< "\n\tForceTime = " << ForceTime << "\n\tHeatTime = " << HeatTime << "\n");
	
			// Update the plasma data from the plasma grid for all models...
			InGrid = FM[i].update_plasmadata(Sample[i]->get_position());
//			InGrid = HM[i].update_plasmadata(Sample[i]->get_position());
//			InGrid = CM[i].update_plasmadata(Sample[i]->get_position());
			Print();
			// ***** START OF : DETERMINE IF END CONDITION HAS BEEN REACHED ***** //
//			std::cout << "\nSample[i]->is_breakingup() = " << Sample[i]->is_breakingup();
			if( Sample[i]->is_gas() && Sample[i]->get_superboilingtemp() <= Sample[i]->get_temperature() ){
				std::cout << "\n\nSample[i] has Boiled";
				break;
			}else if( Sample[i]->is_gas() && Sample[i]->get_superboilingtemp() > Sample[i]->get_temperature() ){
				std::cout << "\n\nSample[i] has Evaporated";
				break;
			}else if( Sample[i]->is_breakingup() ){
				std::cout << "\n\nElectrostatic breakup has occured for Sample[i]";
				std::cin.get();
				HM.push_back(HM[0]);
				HM[i+1].CreateFile("hf" + std::to_string(i+1) + ".txt",false);
				FM.push_back(FM[0]);
				FM[i+1].CreateFile("ff" + std::to_string(i+1) + ".txt");
				CM.push_back(CM[0]);
				CM[i+1].CreateFile("cf" + std::to_string(i+1) + ".txt");
				Matter *NewSample;
				if 	(Sample[i]->get_elem() == 'W') NewSample = new Tungsten	(Sample[i]->get_radius(),Sample[i]->get_temperature(),Sample[i]->get_models());
				else if (Sample[i]->get_elem() == 'B') NewSample = new Beryllium(Sample[i]->get_radius(),Sample[i]->get_temperature(),Sample[i]->get_models());
				else if (Sample[i]->get_elem() == 'F') NewSample = new Iron	(Sample[i]->get_radius(),Sample[i]->get_temperature(),Sample[i]->get_models());
				else if (Sample[i]->get_elem() == 'G') NewSample = new Graphite	(Sample[i]->get_radius(),Sample[i]->get_temperature(),Sample[i]->get_models());
				else{
					std::cerr << "\nInvalid Option Enetered!";
				}
				threevector temp = Sample[i]->get_position();
				threevector temp2(0.0,0.0,-100.0);
				threevector temp3(0.0,0.0,100.0);
				std::cout << "\ntemp = " << temp; std::cin.get();
				Sample.push_back(NewSample);
//				Sample[i]->update_motion(temp,temp3);
				Sample[i+1]->update_motion(temp,temp2);
				std::cout << "\nSample[i]->get_position() = " << Sample[i]->get_position();
				std::cout << "\nSample[i+1]->get_position() = " << Sample[i+1]->get_position();
				Sample[i]->reset_breakup();
				HM[i+1].update_sample(Sample[i+1]);
				HM[i].update_sample(Sample[i]);
				HM[i+1].update_plasmadata(NewSample->get_position());
				PlasmaData *pdatatemp2 = HM[i+1].get_plasmadata();
				CM[i+1].update_sample(Sample[i+1]);
				CM[i].update_sample(Sample[i]);
				CM[i+1].update_plasmadata(pdatatemp2);
				FM[i+1].update_sample(Sample[i+1]);
				FM[i].update_sample(Sample[i]);
				FM[i+1].update_plasmadata(pdatatemp2);

			}else if( Sample[i]->is_gas() ){
				std::cout << "\nSample[i] has vapourised";
				break;	
			}
			// ***** END OF : DETERMINE IF END CONDITION HAS BEEN REACHED ***** //
		}
		if( !InGrid ){
			std::cout << "\nSample has left simulation domain";
		}

	}
	if( fabs(1 - (HM[0].get_totaltime()/FM[0].get_totaltime())) > 0.001  
		|| fabs(1 - (FM[0].get_totaltime()/CM[0].get_totaltime())) > 0.001 ){
		std::cout << "\nWarning! Total Times recorded by processes don't match!";
		std::cout << "\nHM[0].get_totaltime() = " << HM[0].get_totaltime() << "\nFM[0].get_totaltime() = " 
				<< FM[0].get_totaltime() << "\nCM[0].get_totaltime() = " << CM[0].get_totaltime() 
				<< "\n\nTotalTime = " << TotalTime;
	}
	if( HeatTime == 1 ) std::cout << "\nEnd of Run, exiting due to Thermal Equilibrium being reached!\n\n";
	std::cout << "\nFinished DTOKS-U run.";
	return 0;
}

