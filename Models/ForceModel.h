#ifndef __FORCEMODEL_H_INCLUDED__   // if ForceModel.h hasn't been included yet...
#define __FORCEMODEL_H_INCLUDED__

#include "Model.h"

class ForceModel : public Model {

	private:
		// Models defining the Force equation
		enum { NumModels = 6 };
		std::array<bool,NumModels> UseModel; 		// Force Models turned on of possibly 4
		std::string FileName;				// Variable to hold data file name

		void Print();			// Write to output data file

		// Different Force terms
		threevector CalculateAcceleration()const;	// Sum of Force terms
		threevector Gravity()const;
		threevector DTOKSIonDrag()const;
		threevector HybridIonDrag()const;
		threevector NeutralDrag()const;
		threevector LorentzForce()const;
		threevector Centrifugal()const;


	public:
		// Constructors
		ForceModel();
		ForceModel(std::string filename, double accuracy, std::array<bool,NumModels> models, 
				Matter *& sample, PlasmaData *& pdata);
		ForceModel(std::string filename, double accuracy, std::array<bool,NumModels> models, 
				Matter *& sample, PlasmaGrid & pgrid);

		// Destructor
		~ForceModel(){};
		
		void CreateFile(std::string filename);

		double ProbeTimeStep()const;	// Verify time step
		double UpdateTimeStep();	// Verify time step

		// Functions which generate and save data from heating the Sample.
		void Force();
		void Force(double timestep);
};

#endif
