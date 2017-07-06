#ifndef __DTOKSU_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __DTOKSU_H_INCLUDED__

#include "PlasmaGrid.h"
#include "HeatingModel.h"
#include "ForceModel.h"
#include "ChargingModel.h"

#include <algorithm>

class DTOKSU{

	private:
		// Private member data
		double MinTimeStep;			// Seconds, the length of a particular time step
		double TotalTime;			// Seconds, total time taken to perform simulation

		std::vector<Matter*> Sample;		// Matter sample can be either Tungsten, Beryllium, Graphite or Iron

		std::vector<HeatingModel> HM;		// Heating Model 
		std::vector<ForceModel> FM;		// Force Model 
		std::vector<ChargingModel> CM;		// Charge Model

		std::ofstream MyFile;			// Output data file
	
		// Private Functions
		void CheckTimeStep();			// Verify time step
		void Print();				// Write to output data file
		void CreateFile(std::string filename);

	public:
//		DTOKSU();
		DTOKSU( double timestep, std::array<double,3> alvls, Matter *& sample, PlasmaData *&pdata,
				std::array<bool,9> &heatmodels, std::array<bool,5> &forcemodels, std::array<bool,1> &chargemodels);
		DTOKSU( double timestep, std::array<double,3> alvls, Matter *& sample, PlasmaGrid &pgrid,
				std::array<bool,9> &heatmodels, std::array<bool,5> &forcemodels, std::array<bool,1> &chargemodels);

		~DTOKSU(){
		};
		
		int Run();
		// Functions which generate and save data from heating the Sample.

		const double DeltaTherm(double DustTemperature)const;	
		const double DeltaSec()const;
		const double DeltaTot(double DustTemperature)const;
};

#endif
