#ifndef __HEATINGMODEL_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __HEATINGMODEL_H_INCLUDED__

#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

class HeatingModel{

	private:
		// Parameters defining the Heating equation
		const std::string Type;			// Type of Heating Model : 'constant'
		double PowerIncident;			// Kilo-Watts, 
		double TimeStep;			// Seconds, the length of a particular time step
		double TotalTime;			// Seconds, total time taken to perform simulation
		double OldTemp;				// K, Temperature last step
		bool ThermalEquilibrium;		// If the sample is in thermal equilibrium, this variable is true
		bool ForceNegative;			// If we want to force the dust grain to be negative
		Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite

		PlasmaData Pdata;
		
		std::array<bool,9> UseModel; 		// Heating Models turned on of possibly 9

		std::ofstream HeatingFile;		// Output data file

		void CheckTimeStep(double TotalEnergy,char TimeStepType);	// Verify time step
		void Print(double TotPower, std::array<char,4> &constmodels);	// Write to output data file

	public:
		// Constructors
		HeatingModel();
		HeatingModel(std::string name, double power);
		HeatingModel(std::string name, char element, double radius, double temp, double power,
				std::array<bool,9> &models, double timestep );
		HeatingModel(std::string name, char element, double radius, double temp, double power, 
				std::array<bool,9> &models, std::array<char,4> &constmodels, double timestep );
		HeatingModel(std::string name, char element, double power, std::array<bool,9> &models, double timestep, 
				Matter *&sample, PlasmaData &pdata);
		// Destructor
		~HeatingModel(){
			delete Sample;
		};
		
		void Defaults(); // Sets default settings

		// Functions which generate and save data from heating the Sample.
		const int Vapourise(std::string filename, std::array<char,4> &ConstModels, char TimeStepType);
		void Heat(std::array<char,4> &ConstModels, char TimeStepType);
		void Reset(std::string filename, double radius, double temp, double timestep, PlasmaData &pdata, Matter *&sample);
		void CreateFile(std::string filename, bool PrintPhaseData,std::array<char,4> &constmodels);
		double CalculatePower(double DustTemperature)const;
		double RungeKutta4();

		// Heating Models
		const double EmissivityModel		(double DustTemperature)	const;
		const double EvaporationModel		(double DustTemperature)	const;
		const double NewtonCooling		(double DustTemperature)	const;
		const double SEE			(double DustTemperature)	const;
		const double TEE			(double DustTemperature)	const;
		const double NeutralRecombination	(double DustTemperature)	const;
		const double IonHeatFlux		(double DustTemperature)	const;
		const double ElectronHeatFlux		(double DustTemperature)	const;
		const double NeutralHeatFlux		()				const;

		// Fluxes of particles and coefficients
		const double EvaporationFlux		(double DustTemperature)	const;
		const double DeltaTherm			(double DustTemperature)	const;
		const double DeltaSec			()				const;
		const double DeltaTot			(double DustTemperature)	const;
		const double IonFlux			(double DustTemperature)	const;
		const double ElectronFlux		(double DustTemperature)	const;
		const double NeutralFlux		()				const;

		double get_totaltime			()const{ return TotalTime;	};
};

#endif
