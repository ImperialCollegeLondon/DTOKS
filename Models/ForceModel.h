#ifndef __FORCEMODEL_H_INCLUDED__   // if ForceModel.h hasn't been included yet...
#define __FORCEMODEL_H_INCLUDED__

#include "Model.h"

class ForceModel : public Model {

	private:
		// Parameters defining the Force equation
		std::array<bool,5> UseModel; 		// Force Models turned on of possibly 4

		void Print();			// Write to output data file

		// Different Force terms
		threevector CalculateAcceleration()const;	// Sum of Force terms
		threevector DTOKSIonDrag()const;
		threevector NeutralDrag()const;
		threevector LorentzForce()const;
		threevector Centrifugal()const;


	public:
		// Constructors
		ForceModel();
		ForceModel(std::string filename, double accuracy, std::array<bool,5> models, 
				Matter *& sample, PlasmaData *& pdata);
		ForceModel(std::string filename, double accuracy, std::array<bool,5> models, 
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
