#ifndef __FORCEMODEL_H_INCLUDED__   // if ForceModel.h hasn't been included yet...
#define __FORCEMODEL_H_INCLUDED__

#include "Model.h"

class ForceModel : public Model {

	private:
		// Parameters defining the Force equation
		double TimeStep;
		double TotalTime;		
		std::array<bool,4> UseModel; 		// Force Models turned on of possibly 4
		std::ofstream ForceFile;		// Output data file

		void Print();			// Write to output data file
		
		// Different Force terms
		threevector DTOKSIonDrag()const;
		threevector LorentzForce()const;
		threevector Centrifugal()const;

	public:
		// Constructors
		ForceModel();
		ForceModel(std::string filename, double accuracy, std::array<bool,4> models, 
				std::shared_ptr <Matter> const& sample, PlasmaData const& pdata);

		// Destructor
		~ForceModel(){
		};
		
		double CheckTimeStep();	// Verify time step

		// Functions which generate and save data from heating the Sample.
		void CreateFile(std::string filename);

		void Force();
};

#endif
