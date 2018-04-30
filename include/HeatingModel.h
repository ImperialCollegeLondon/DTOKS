#ifndef __HEATINGMODEL_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __HEATINGMODEL_H_INCLUDED__

#include "Model.h"

class HeatingModel : public Model{

	private:
		// Parameters defining the Heating equation
		double PowerIncident;			// Kilo-Watts, 
		double OldTemp;				// Kelvin, temperature last step, used to determine TE
		double RE;				// Fraction of backscattered energy
		double RN;				// Fraction of backscattered particles
		bool ForceNegative;			// If we want to force the dust grain to be negative
		bool ThermalEquilibrium;		// If the Dust grain is in Thermal Equilibrium (Constant Plasma only)
		std::string FileName;			// Variable to hold data file name

		// Models defining the heating equation
		enum { NumModels = 9 };
		std::array<bool,NumModels> UseModel; 		// Heating Models turned on of possibly 9

		void Print();				// Write to output data file
		void Defaults(); // Sets default settings
		double RungeKutta4(double timestep);

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

	public:
		// Constructors
		HeatingModel();
		HeatingModel( std::string filename, float accuracy, std::array<bool,NumModels> &models, 
				Matter *& sample, PlasmaData *& pdata);
		HeatingModel( std::string filename, float accuracy, std::array<bool,NumModels> &models, 
				Matter *& sample, PlasmaGrid & pgrid);

		// Destructor
		~HeatingModel(){
		};
		
		void CreateFile(std::string filename, bool PrintPhaseData);

		// Functions which generate and save data from heating the Sample.
		const int Vapourise();
		void Heat();
		void Heat(double timestep);

		double CalculatePower(double DustTemperature)const;

		double ProbeTimeStep()const;		// Check time step
		double UpdateTimeStep();		// Update time step

		bool get_thermalequilibrium		(){ return ThermalEquilibrium; }
		void set_PowerIncident			(double powerincident){ PowerIncident = powerincident; 		}
		void set_ThermalEquilibrium		(double thermequilib ){ ThermalEquilibrium = thermequilib; 	}
};

#endif
