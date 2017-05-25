#ifndef __DTOKSU_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __DTOKSU_H_INCLUDED__

#include "HeatingModel.h"
#include "ForceModel.h"
#include "ChargingModel.h"

class DTOKSU{

	private:
		double TimeStep;			// Seconds, the length of a particular time step
		double TotalTime;			// Seconds, total time taken to perform simulation

		std::shared_ptr <Matter> Sample;	// Tungsten, Beryllium, Iron or Graphite
		HeatingModel HM;
		ForceModel FM;
		ChargingModel CM;
		PlasmaData Pdata;
		
		std::array<bool,9> HeatingSwitch; 	// Heating Models turned on of possibly 9
		std::array<bool,9> ForceSwitch; 	// Force Models turned on of possibly N
		std::array<bool,9> ChargingSwitch; 	// Charging Models turned on of possibly N
		std::ofstream MyFile;			// Output data file

		void CheckTimeStep();			// Verify time step
		void Print();				// Write to output data file
		void CreateFile(std::string filename);

	public:
		DTOKSU();
		DTOKSU( double timestep, std::shared_ptr<Matter> const& sample, PlasmaData const &pdata,
				std::array<bool,9> &heatmodels, std::array<bool,3> &forcemodels, std::array<bool,1> &chargemodels);

		~DTOKSU(){
		};
		
		int Run();
		// Functions which generate and save data from heating the Sample.

		const double DeltaTherm(double DustTemperature)const;	
		const double DeltaSec()const;
		const double DeltaTot(double DustTemperature)const;


};

#endif
