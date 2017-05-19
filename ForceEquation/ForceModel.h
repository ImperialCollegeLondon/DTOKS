#ifndef __FORCEMODEL_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __FORCEMODEL_H_INCLUDED__


#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

class ForceModel{

	private:
		// Parameters defining the Heating equation
		Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite
		PlasmaData Pdata;
		
		std::array<bool,9> UseModel; 		// Force Models turned on of possibly 9
		std::ofstream ForceFile;		// Output data file

		void Print();				// Write to output data file

	public:
		// Constructors
		ForceModel();

		// Destructor
		~ForceModel(){
			delete Sample;
		};
		
		// Functions which generate and save data from heating the Sample.
		void CreateFile(std::string filename);
};

#endif
