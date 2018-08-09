//#define PAUSE
//#define BREAKUP_DEBUG
#include "Breakup.h"

// Parameterised constructor
Breakup::Breakup( DTOKSU *&dtoksu, Matter *&sample ) : Sample(sample),Sim(dtoksu){
}

// Run DTOKS many times with breakup turned on
int Breakup::Run(){
	threevector Zeros(0.0,0.0,0.0);

	unsigned int p(1);
	unsigned int i(1);
	
	double seed=std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937 randnumber(seed);
	std::uniform_real_distribution<double> rad(0.0, 1.0); // Uniformly Randomly Distributed Variable between 0.0 and 1.0
	// Loop over the number of split paths.
	for( unsigned int j(0);	j < i; j ++){
		B_Debug("\nSimulating NEGATIVE Branch "); B_Debug(i); B_Debug(" : "); B_Debug(j);
		B_Debug("\nStart Pos = "); B_Debug(Sample->get_position()); 
		B_Debug("\nVelocity = "); B_Debug(Sample->get_velocity());
		B_Debug("\nMass = "); B_Debug(Sample->get_mass()); 
		B_Debug("\nTemperature = "); B_Debug(Sample->get_temperature());

		// When breakup occurs and a path forks, track it. If it breaks up, track the subsequent particle
		// Repeat until the end condition is no-longer breakup, i.e while return of DTOKSU object isn't 3.
		while( Sim->Run() == 3 ){ // Breakup has occured...

			// Close data files and open new ones, with new names based off index
			Sim->CloseFiles();
			Sim->OpenFiles("Data/breakup",p);
			
			// Reset breakup so that it's recorded with breakup turned off
			Sample->reset_breakup();

			// Reset the end point data with the same position, no rotation and heading off in negative direction
			// Rotation occurs in random direction perpendicular to magnetic field and velocity as per theory
			double VelocityMag = 2*PI*(Sample->get_radius())*Sample->get_rotationalfreq();
			threevector Unit(2.0*rad(randnumber)-1.0,2.0*rad(randnumber)-1.0,2.0*rad(randnumber)-1.0);
			threevector VelocityUnitVec = (Unit.getunit()^Sim->get_bfielddir()).getunit();
			threevector dvMinus = VelocityMag*VelocityUnitVec;
			Sample->update_motion(Zeros,dvMinus,-Sample->get_rotationalfreq());

			// Record dust end conditions for later when we simulate other half of fork
			GDvector.push_back(Sample->get_graindata());
			HMTime.push_back(Sim->get_HMTime());
			FMTime.push_back(Sim->get_FMTime());
			CMTime.push_back(Sim->get_CMTime());
	
			// Change the dust velocity, the mass has already been halved in Matter. 
			// Add the velocity twice over as we took it away in one direction
			threevector dvPlus = -2.0*dvMinus;
			Sample->update_motion(Zeros,dvPlus,0.0);

			B_Debug("\nSimulating POSITIVE Branch "); B_Debug(i); B_Debug(" : "); B_Debug(j);
			B_Debug("\nStart Pos = "); B_Debug(Sample->get_position());
			B_Debug("\nVelocity = "); B_Debug(Sample->get_velocity());
			B_Debug("\nMass = "); B_Debug(Sample->get_mass()); 
			B_Debug("\nTemperature = "); B_Debug(Sample->get_temperature()); 
			// increment counters of number of forks and number of positive forks.
			// p is used for recording index of file
			i = i + 1;
			p = p + 1;
		} // Run the simulation again if breakup occured!
		
		B_Debug("\n***** START OF : DUST DIDN'T BREAKUP *****\n!");
		B_Debug(i); B_Debug(" : "); B_Debug(j);
		// Close data files and open new ones
		Sim->CloseFiles();
		if( GDvector.size() > 0 ){ // If we have at least one breakup event, i.e stored end point data
			// Re-initialise simulation with new Velocity... 
			// But we've already done this now when we saved the data previously...
			// So we're good to go, just copy over the previous end conditions
			Sample->update_graindata(GDvector[j]);


			// Reset Model Times, the only reason for this is to make the plotting work correctly...
			// Ensure that the time of the models is global and not local to the track
			Sim->ResetModelTime(HMTime[j]-Sim->get_HMTime(),FMTime[j]-Sim->get_FMTime(),CMTime[j]-Sim->get_CMTime());

			// Open files 
			Sim->OpenFiles("Data/breakup",p);
			p = p + 1;
			B_Debug("\nSimulating POSITIVE Branch "); B_Debug(i); B_Debug(" : "); B_Debug(j);
                        B_Debug("\nStart Pos = "); B_Debug(Sample->get_position()); 
                        B_Debug("\nVelocity = "); B_Debug(Sample->get_velocity());
                        B_Debug("\nMass = "); B_Debug(Sample->get_mass());
                        B_Debug("\nTemperature = "); B_Debug(Sample->get_temperature());
			Pause();
		}
		B_Debug("\n***** END OF : DUST DIDN'T BREAKUP *****\n!");
	}
	
}
