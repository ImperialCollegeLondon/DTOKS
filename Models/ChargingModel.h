#ifndef __CHARGINGMODEL_H_INCLUDED__   // if ChargingModel.h hasn't been included yet...
#define __CHARGINGMODEL_H_INCLUDED__

#include "Model.h"

class ChargingModel : public Model{

	private:
		// Parameters defining the Heating equation
		double TotalTime;				
		double TimeStep;

		std::array<bool,1> UseModel; 		// Charging Models turned on of possibly 9

		void Print();			// Write to output data file
		double solveOML(double a, double guess);

	public:
		// Constructors
		ChargingModel();
		ChargingModel(std::string filename, double accuracy, std::array<bool,1> models, 
				Matter *& sample, PlasmaData & pdata);
		ChargingModel(std::string filename, double accuracy, std::array<bool,1> models, 
				Matter *& sample, PlasmaGrid & pgrid);

		// Destructor
		~ChargingModel(){
		};
		
		double CheckTimeStep();
		// Functions which generate and save data from heating the Sample.
		void CreateFile(std::string filename);
		
		void Charge();


};

#endif
