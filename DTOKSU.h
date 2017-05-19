#ifndef __DTOKSU_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __DTOKSU_H_INCLUDED__

#include "HeatingModel.h"
#include "ForceModel.h"
#include "ChargingModel.h"

class DTOKSU{

	private:
		double TimeStep;			// Seconds, the length of a particular time step
		double TotalTime;			// Seconds, total time taken to perform simulation

		Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite
		HeatingModel HM;
		ForceModel FM;
		ChargingModel CM;
		PlasmaData Pdata;
		
		std::array<bool,9> HeatingSwitch; 	// Heating Models turned on of possibly 9
		std::ofstream MyFile;			// Output data file

		void CheckTimeStep();			// Verify time step
		void Print();				// Write to output data file
		void CreateFile(std::string filename);

	public:
		DTOKSU();

		~DTOKSU(){
			delete Sample;
		};
		
		int Run();
		// Functions which generate and save data from heating the Sample.

};

#endif
