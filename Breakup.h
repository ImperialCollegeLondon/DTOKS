#ifndef __BREAKUP_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __BREAKUP_H_INCLUDED__

#include "DTOKSU.h"
#include <random>	// for std::normal_distribution<> etc.
#include <chrono>	// for chrono::high_resolution_clock::now().time_since_epoch().count();

class Breakup{

	private:
		// Private member data
		std::vector<GrainData> GDvector;
		std::vector<double> HMTime;
		std::vector<double> FMTime;
		std::vector<double> CMTime;
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
