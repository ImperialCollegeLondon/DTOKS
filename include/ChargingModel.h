#ifndef __CHARGINGMODEL_H_INCLUDED__   // if ChargingModel.h hasn't been included yet...
#define __CHARGINGMODEL_H_INCLUDED__

#include "Model.h"

const unsigned int CMN = 10;	// CHARGE MODEL NUMBER, the number of charge models

class ChargingModel : public Model{

	private:
		// Models defining the Force equation
		std::array<bool,CMN> UseModel; 		// Charging Models turned on of possible 1
		std::string FileName;						// Variable to hold data file name
		
		void Print();			// Write to output data file

		
		double solveOML(double a, double guess);
		double solvePosOML(double a, double guess);
		double solveNegSchottkyOML(double guess);
		double solvePosSchottkyOML();
		double solvePHL(double a, double guess);
		double solveMOML();
		double solveSOML(double a, double guess);
		double solveSMOML(double a, double guess);
		double solveMOMLWEM(double );
		//double solveMOMLEM();
		double DeltaTherm()const;
		double DeltaSec()const;

	public:
		// Constructors
		ChargingModel();
		ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models, 
				Matter *& sample, PlasmaData & pdata);
		ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models, 
				Matter *& sample, PlasmaData * pdata);
		ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models, 
				Matter *& sample, PlasmaGrid_Data & pgrid);
		ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models, 
				Matter *& sample, PlasmaGrid_Data & pgrid, PlasmaData &pdata);

		void CreateFile(std::string filename);

		// Destructor
		~ChargingModel(){
		};
		
		double ProbeTimeStep()const;
		double UpdateTimeStep();
		// Functions which generate and save data from heating the Sample.

		
		void Charge();
		void Charge(double timestep);
};

#endif
