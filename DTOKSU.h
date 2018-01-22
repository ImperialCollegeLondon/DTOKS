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
		double TotalTime;			// Seconds, total time taken to perform simulation

		Matter *Sample;				// Matter sample can be either Tungsten, Beryllium, Graphite or Iron

		HeatingModel HM;			// Heating Model 
		ForceModel FM;				// Force Model 
		ChargingModel CM;			// Charge Model

		std::ofstream MyFile;			// Output data file
	
		// Private Functions
		void Print();				// Write to output data file
		void CreateFile(std::string filename);

	public:
//		DTOKSU();
		DTOKSU( std::array<float,3> alvls, Matter *& sample, PlasmaData *&pdata,
				std::array<bool,9> &heatmodels, std::array<bool,6> &forcemodels, std::array<bool,3> &chargemodels);
		DTOKSU( std::array<float,3> alvls, Matter *& sample, PlasmaGrid &pgrid,
				std::array<bool,9> &heatmodels, std::array<bool,6> &forcemodels, std::array<bool,3> &chargemodels);

		~DTOKSU(){
		};
		
		int Run();
		void OpenFiles(std::string filename, unsigned int i);
		void CloseFiles();			// Close all model files
		void ResetModelTime(double HMTime, double FMTime, double CMTime);
	
		double get_HMTime()const{ return HM.get_totaltime(); }
		double get_FMTime()const{ return FM.get_totaltime(); }
		double get_CMTime()const{ return CM.get_totaltime(); }
	
		threevector get_bfielddir(){
			return (FM.get_plasmadata()->MagneticField.getunit());
		}
};

#endif
