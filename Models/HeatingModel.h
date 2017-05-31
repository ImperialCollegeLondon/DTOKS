#ifndef __HEATINGMODEL_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __HEATINGMODEL_H_INCLUDED__

#include "Model.h"

class HeatingModel : public Model{

	private:
		// Parameters defining the Heating equation
		const std::string Type;			// Type of Heating Model : 'constant'
		double PowerIncident;			// Kilo-Watts, 
		double TotalPower;			// Kilo-Watts,
		double TimeStep;			// Seconds, the length of a particular time step
		double TotalTime;			// Seconds, total time taken to perform simulation
		bool ForceNegative;			// If we want to force the dust grain to be negative

		std::array<bool,9> UseModel; 		// Heating Models turned on of possibly 9

		void Print();				// Write to output data file

	public:
		// Constructors
		HeatingModel();
		HeatingModel( std::string filename, double timestep, double accuracy, std::array<bool,9> &models, 
				std::shared_ptr<Matter> const& sample, PlasmaData const& pdata);

		// Destructor
		~HeatingModel(){
		};
		
		void Defaults(); // Sets default settings

		// Functions which generate and save data from heating the Sample.
		const int Vapourise();
		void Heat();

		void CreateFile(std::string filename, bool PrintPhaseData);
		double CalculatePower(double DustTemperature)const;
		double RungeKutta4();


		double CheckTimeStep();	// Verify time step

		// Heating Models
		const double EmissivityModel		(double DustTemperature)const;
		const double EvaporationModel		(double DustTemperature)const;
		const double NewtonCooling		(double DustTemperature)const;
		const double SEE			(double DustTemperature)const;
		const double TEE			(double DustTemperature)const;
		const double NeutralRecombination	(double DustTemperature)const;
		const double IonHeatFlux		(double DustTemperature)const;
		const double ElectronHeatFlux		(double DustTemperature)const;
		const double NeutralHeatFlux		()			const;

		// Fluxes of particles and coefficients
		const double EvaporationFlux		(double DustTemperature)const;
		const double IonFlux			(double DustTemperature)const;
		const double ElectronFlux		(double DustTemperature)const;
		const double NeutralFlux		()			const;

		double get_totaltime			()const{ return TotalTime;	};
};

#endif
