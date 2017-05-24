#ifndef __FORCEMODEL_H_INCLUDED__   // if ForceModel.h hasn't been included yet...
#define __FORCEMODEL_H_INCLUDED__

//#include "threevector.h"
#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

class ForceModel{

	private:
		// Parameters defining the Force equation
		double TimeStep;

		Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite
		PlasmaData Pdata;
		
		std::array<bool,3> UseModel; 		// Force Models turned on of possibly 3
		std::ofstream ForceFile;		// Output data file

		void Print();			// Write to output data file
		
		// Different Force terms
		threevector DTOKSIonDrag()const;
		threevector LorentzForce()const;
		threevector Centrifugal()const;

	public:
		// Constructors
		ForceModel();

		// Destructor
		~ForceModel(){
			delete Sample;
		};
		
		// Functions which generate and save data from heating the Sample.
		void CreateFile(std::string filename);

		void Force();
};

#endif
