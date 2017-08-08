//#define PAUSE
//#define BREAKUP_DEBUG
//#define BREAKUP_DEEP_DEBUG
#include "Breakup.h"

// Parameterised constructor
Breakup::Breakup( DTOKSU *&dtoksu, Matter *&sample ) : Sample(sample){
	Sim = dtoksu;
}

//
int Breakup::Run(){
	threevector Zeros(0.0,0.0,0.0);
	threevector dvPlus(0.0,0.0,10.0);
	threevector dvMinus(0.0,0.0,-10.0);


	unsigned int p(1);
	unsigned int i(1);

	for( unsigned int j(0);	j < i; j ++){
		while( Sim->Run() == 3 ){ // Breakup has occured...
			std::cout << "\ni : " << i << "\tj : " << j;
			// Close data files and open new ones
			Sim->CloseFiles();
			Sim->OpenFiles("Data/breakup",p);
			
			// Reset breakup so that it's recorded with breakup turned off
			Sample->reset_breakup();

			// Record dust end conditions
			GDvector.push_back(Sample->get_graindata());			
	
			// Change the dust velocity, the mass has already been halved in Matter. 
			Sample->update_motion(Zeros,dvPlus,0.0);


//			std::cout << "\nSimulating Branch " << i << "\nStart Pos = " << Sample->get_position()
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
			Sim->OpenFiles("Data/breakup",p);

			// Re-initialise simulation with new Velocity...
			GDvector[j].DustVelocity = GDvector[j].DustVelocity + dvMinus;

			Sample->update_graindata(GDvector[j]);
			p = p + 1;
		}
		std::cout << "\nSimulating Branch " << j+1 << "\nStart Pos = " << GDvector[j].DustPosition
			<< "\nVelocity = " << Sample->get_velocity() << "\nMass = " << Sample->get_mass(); std::cin.get();
//		std::cout << "\n***** END OF : DUST DIDN'T BREAKUP *****\n!";
	}
	
}
