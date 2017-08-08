#ifndef __BREAKUP_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __BREAKUP_H_INCLUDED__

#include "DTOKSU.h"

class Breakup{

	private:
		// Private member data
		std::vector<GrainData> GDvector; 
		DTOKSU *Sim;
		Matter *Sample;
	public:
//		Breakup();
		Breakup( DTOKSU *& dtoksu, Matter *& sample );

		~Breakup(){
		};
		
		int Run();
};

#endif
