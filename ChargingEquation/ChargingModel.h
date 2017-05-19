#ifndef __CHARGINGMODEL_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __CHARGINGMODEL_H_INCLUDED__

#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

class ChargingModel{

	private:
		// Parameters defining the Heating equation
		Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite
		PlasmaData Pdata;
		
		std::array<bool,9> UseModel; 		// Charging Models turned on of possibly 9
		std::ofstream ChargingFile;		// Output data file

		void Print();				// Write to output data file

	public:
		// Constructors
		ChargingModel();

		// Destructor
		~ChargingModel(){
			delete Sample;
		};
		
		// Functions which generate and save data from heating the Sample.
		void CreateFile(std::string filename);
};

#endif
