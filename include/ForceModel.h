#ifndef __FORCEMODEL_H_INCLUDED__   // if ForceModel.h hasn't been included yet...
#define __FORCEMODEL_H_INCLUDED__

#include "Model.h"

const unsigned int FMN = 6;	// CHARGE MODEL NUMBER, the number of charge models

class ForceModel : public Model {

	private:
		// Models defining the Force equation
		std::array<bool,FMN> UseModel; 		// Force Models turned on of possibly 4
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
		ForceModel(std::string filename, float accuracy, std::array<bool,FMN> models, 
				Matter *& sample, PlasmaData & pdata);
		ForceModel(std::string filename, float accuracy, std::array<bool,FMN> models, 
				Matter *& sample, PlasmaData * pdata);
		ForceModel(std::string filename, float accuracy, std::array<bool,FMN> models, 
				Matter *& sample, PlasmaGrid_Data & pgrid);

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
