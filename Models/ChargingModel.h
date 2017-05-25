#ifndef __CHARGINGMODEL_H_INCLUDED__   // if ChargingModel.h hasn't been included yet...
#define __CHARGINGMODEL_H_INCLUDED__

#include "Model.h"

class ChargingModel : public Model{

	private:
		// Parameters defining the Heating equation
				
		std::array<bool,1> UseModel; 		// Charging Models turned on of possibly 9

		void Print();			// Write to output data file
		double solveOML(double a, double guess);

	public:
		// Constructors
		ChargingModel();
		ChargingModel(std::string filename, std::array<bool,1> models, std::shared_ptr <Matter> const& sample, 
				PlasmaData const& pdata);

		// Destructor
		~ChargingModel(){
		};
		
		// Functions which generate and save data from heating the Sample.
		void CreateFile(std::string filename);

		void Charge();


};

#endif
