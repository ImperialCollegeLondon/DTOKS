//#define PAUSE
//#define BREAKUP_DEBUG
//#define BREAKUP_DEEP_DEBUG
#include "Breakup.h"

// Parameterised constructor
Breakup::Breakup( DTOKSU *&dtoksu, Matter *&sample ) : Sample(sample),Sim(dtoksu){
}

//
int Breakup::Run(){
	threevector Zeros(0.0,0.0,0.0);

	unsigned int p(1);
	unsigned int i(1);

	// Loop over the number of broken paths.
	for( unsigned int j(0);	j < i; j ++){
//		std::cout << "\nSimulating NEGATIVE Branch " << i << " : " << j << "\nStart Pos = " << Sample->get_position()
//			<< "\nVelocity = " << Sample->get_velocity() << "\nMass = " << Sample->get_mass(); std::cin.get();

		// When breakup occurs and a path forks, track it. If it breaks up, track the subsequent particle
		// Repeat until the end condition is no-longer breakup
		while( Sim->Run() == 3 ){ // Breakup has occured...

			// Close data files and open new ones, with new names based off index
			Sim->CloseFiles();
			Sim->OpenFiles("Data/breakup",p);
			
			// Reset breakup so that it's recorded with breakup turned off
			Sample->reset_breakup();

			// Reset the end point data with the same position, no rotation and heading off in negative direction
			// Rotation occurs in same direction as Ion Gyromotion
			double VelocityMag = 2*PI*(Sample->get_radius())*Sample->get_rotationalfreq();

			threevector VelocityUnitVec = (Sample->get_velocity().getunit()^Sim->get_bfielddir());
			threevector dvMinus = -VelocityMag*VelocityUnitVec;
			Sample->update_motion(Zeros,dvMinus,-Sample->get_rotationalfreq());

			// Record dust end conditions
			GDvector.push_back(Sample->get_graindata());
			HMTime.push_back(Sim->get_HMTime());
			FMTime.push_back(Sim->get_FMTime());
			CMTime.push_back(Sim->get_CMTime());
	
			// Change the dust velocity, the mass has already been halved in Matter. 
			// Add the velocity twice over as we took it away in one direction
			threevector dvPlus = 2*VelocityMag*VelocityUnitVec;
			Sample->update_motion(Zeros,dvPlus,0.0);



//			std::cout << "\nSimulating POSITIVE Branch " << i << " : " << j << "\nStart Pos = " << Sample->get_position()
//				<< "\nVelocity = " << Sample->get_velocity() << "\nMass = " << Sample->get_mass(); std::cin.get();
			i = i + 1;
			p = p + 1;
			// Run the simulation again!
		}
		
//		std::cout << "\n***** START OF : DUST DIDN'T BREAKUP *****\n!";
//		std::cout << "\ni : " << i << "\tj : " << j;
		// Close data files and open new ones
		Sim->CloseFiles();
		if( GDvector.size() > 0 ){
			// Re-initialise simulation with new Velocity... But we've already done this now when we saved the data previously...
			// So we're good to go, just copy over the previous end conditions
			Sample->update_graindata(GDvector[j]);


			// Reset Model Times, the only reason for this is to make the plotting work correctly...
			Sim->ResetModelTime(HMTime[j]-Sim->get_HMTime(),FMTime[j]-Sim->get_FMTime(),CMTime[j]-Sim->get_CMTime());

			Sim->OpenFiles("Data/breakup",p);
			p = p + 1;
//			std::cout << "\nSimulating Branch " << i  << " : " << j << "\nStart Pos = " << GDvector[j].DustPosition
//				<< "\nVelocity = " << Sample->get_velocity() << "\nMass = " << Sample->get_mass(); std::cin.get();

		}
//		std::cout << "\n***** END OF : DUST DIDN'T BREAKUP *****\n!";
	}
	
}
