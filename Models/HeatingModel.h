#ifndef __HEATINGMODEL_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __HEATINGMODEL_H_INCLUDED__

#include "Model.h"

class HeatingModel : public Model{

	private:
		// Parameters defining the Heating equation
		double PowerIncident;			// Kilo-Watts, 
		double OldTemp;				// Kelvin, temperature last step, used to determine TE
		bool ForceNegative;			// If we want to force the dust grain to be negative
		bool ThermalEquilibrium;		// If the Dust grain is in Thermal Equilibrium (Constant Plasma only)

		std::array<bool,9> UseModel; 		// Heating Models turned on of possibly 9

		void Print();				// Write to output data file

	public:
		// Constructors
		HeatingModel();
		HeatingModel( std::string filename, double accuracy, std::array<bool,9> &models, 
				Matter *& sample, PlasmaData *& pdata);
		HeatingModel( std::string filename, double accuracy, std::array<bool,9> &models, 
				Matter *& sample, PlasmaGrid & pgrid);

		// Destructor
		~HeatingModel(){
		};
		
		void Defaults(); // Sets default settings

		// Functions which generate and save data from heating the Sample.
		const int Vapourise();
		void Heat();
		void Heat(double timestep);

		void CreateFile(std::string filename, bool PrintPhaseData);
		double CalculatePower(double DustTemperature)const;
		double RungeKutta4();

		double ProbeTimeStep()const;		// Check time step
		double UpdateTimeStep();		// Update time step

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

		bool get_thermalequilibrium		(){ return ThermalEquilibrium; }
		void set_PowerIncident			(double powerincident){ PowerIncident = powerincident; 		}
		void set_ThermalEquilibrium		(double thermequilib ){ ThermalEquilibrium = thermequilib; 	}
};

#endif
