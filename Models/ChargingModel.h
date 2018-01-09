#ifndef __CHARGINGMODEL_H_INCLUDED__   // if ChargingModel.h hasn't been included yet...
#define __CHARGINGMODEL_H_INCLUDED__

#include "Model.h"

class ChargingModel : public Model{

	private:
		
		std::array<bool,3> UseModel; 		// Charging Models turned on of possible 1

		void Print();			// Write to output data file

		double solveOML(double a, double guess);
		double solveNegSchottkyOML(double guess);
		double solvePosSchottkyOML();
		double DeltaTherm()const;
		double DeltaSec()const;

		double ChargeOfGrain;		// Coulombs, Charge on dust grain
	public:
		// Constructors
		ChargingModel();
		ChargingModel(std::string filename, double accuracy, std::array<bool,3> models, 
				Matter *& sample, PlasmaData *& pdata);
		ChargingModel(std::string filename, double accuracy, std::array<bool,3> models, 
				Matter *& sample, PlasmaGrid & pgrid);

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
