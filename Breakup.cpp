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
	threevector dvPlus(0.0,0.0,10.0);
	threevector dvMinus(0.0,0.0,-10.0);

	unsigned int i(1);

	for( unsigned int j(0);	j < i; j ++){
//		int ReturnStat = Sim->Run();
//		if( ReturnStat == 3 ){ i = i+1; }
		while( Sim->Run() == 3 ){ // Breakup has occured...
			std::cout << "\ni : " << i << "\tj : " << j;
			// Close data files and open new ones
			Sim->CloseFiles();
			Sim->OpenFiles("Data/breakup",i);
			
			// Record dust end conditions
			EndPositions.push_back(Sample->get_position());
			EndVelocities.push_back(Sample->get_velocity());
			EndMasses.push_back(Sample->get_mass());
	
			// Re-initialise simulation with new 
//			threevector StartPos = EndPositions[i-1]-Sample->get_position();
			threevector StartPos = EndPositions[i-1];
	//		threevector StartVel = EndVelocities[i]+dvPlus;
			Sample->update_motion(StartPos,dvPlus);
			Sample->set_mass(EndMasses[i-1]);
			Sample->reset_breakup();

			std::cout << "\nSimulating Branch " << i << "\nStart Pos = " << StartPos 
				<< "\nVelocity = " << Sample->get_velocity() << "\nMass = " << Sample->get_mass(); std::cin.get();

			i = i + 1;



//			std::cout << "\ni : " << i << "\tj : " << j;

			// Run the simulation again!
		}
		std::cout << "\nIN LOOP!";
		// Close data files and open new ones
		Sim->CloseFiles();
		Sim->OpenFiles("Data/breakup_" + std::to_string(i),j);

		// Re-initialise simulation with new conditions...
		threevector StartPos = EndPositions[j]-Sample->get_position();
//		threevector StartVel = EndVelocities[j]+dvMinus;
		Sample->update_motion(StartPos,dvMinus);
		Sample->set_mass(EndMasses[j]);

		std::cout << "\nSimulating Branch " << j+1 << "\nStart Pos = " << StartPos 
			<< "\nVelocity = " << Sample->get_velocity() << "\nMass = " << Sample->get_mass(); std::cin.get();
	}
	
}
