#ifndef __HEATINGMODEL_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __HEATINGMODEL_H_INCLUDED__

#include "Model.h"

const unsigned int HMN = 18;	// CHARGE MODEL NUMBER, the number of charge models

class HeatingModel : public Model{

	private:
		// Parameters defining the Heating equation
		double PowerIncident;			// Kilo-Watts, 
		double OldTemp;				// Kelvin, temperature last step, used to determine TE
		double RE;				// Fraction of backscattered energy
		double RN;				// Fraction of backscattered particles
		bool ThermalEquilibrium;		// If the Dust grain is in Thermal Equilibrium (Constant Plasma only)
		std::string FileName;			// Variable to hold data file name

		// Models defining the heating equation
		std::array<bool,HMN> UseModel; 		// Heating Models turned on of possibly 9

		void Print();				// Write to output data file
		void Defaults(); // Sets default settings
		double RungeKutta4(double timestep);

		// Fluxes of particles and coefficients
		const double EmissivityModel			(double DustTemperature)const;
		const double EvaporationFlux			(double DustTemperature)const;

		// Heating Models
		const double EvaporationModel			(double DustTemperature)const;
		const double NewtonCooling				(double DustTemperature)const;
		const double NeutralHeatFlux			(double DustTemperature)const;

		const double SOMLIonHeatFlux			(double DustTemperature)const;
		const double SOMLNeutralRecombination	(double DustTemperature)const;

		const double SMOMLIonHeatFlux			(double DustTemperature)const;
		const double SMOMLNeutralRecombination	(double DustTemperature)const;

		const double SEE						(double DustTemperature)const;
		const double TEE						(double DustTemperature)const;
		const double PHLElectronHeatFlux		(double DustTemperature)const;

		const double OMLElectronHeatFlux		(double DustTemperature)const;

		const double DTOKSSEE					(double DustTemperature)const;
		const double DTOKSTEE					(double DustTemperature)const;
		const double DTOKSIonHeatFlux			(double DustTemperature)const;
		const double DTOKSNeutralRecombination	(double DustTemperature)const;
		const double DTOKSElectronHeatFlux		(double DustTemperature)const;

		const double DUSTTIonHeatFlux			(double DustTemperature)const;

	public:
		// Constructors
		HeatingModel();
		HeatingModel( std::string filename, float accuracy, std::array<bool,HMN> &models, 
				Matter *& sample, PlasmaData & pdata);
		HeatingModel( std::string filename, float accuracy, std::array<bool,HMN> &models, 
				Matter *& sample, PlasmaData * pdata);
		HeatingModel( std::string filename, float accuracy, std::array<bool,HMN> &models, 
				Matter *& sample, PlasmaGrid_Data & pgrid);
		HeatingModel( std::string filename, float accuracy, std::array<bool,HMN> &models, 
				Matter *& sample, PlasmaGrid_Data & pgrid, PlasmaData &pdata);
		// Destructor
		~HeatingModel(){
		};
		
		void CreateFile(std::string filename, bool PrintPhaseData);

		// Functions which generate and save data from heating the Sample.
		const int Vapourise();
		void Heat();
		void Heat(double timestep);

		void UpdateRERN();
		double CalculatePower(double DustTemperature)const;

		double ProbeTimeStep()const;		// Check time step
		double UpdateTimeStep();		// Update time step

		bool get_thermalequilibrium		(){ return ThermalEquilibrium; }
		void set_PowerIncident			(double powerincident){ PowerIncident = powerincident; 		}
		void set_ThermalEquilibrium		(double thermequilib ){ ThermalEquilibrium = thermequilib; 	}
};

#endif
