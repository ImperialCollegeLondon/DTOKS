#ifndef __MODEL_H_INCLUDED__   // if Model.h hasn't been included yet...
#define __MODEL_H_INCLUDED__

#include <memory>
#include <iomanip>     // std::ofstream::setprecision(), This allows for precision to be set for output

#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

static struct PlasmaData PlasmaDataDefaults = {
	1e19,				// m^-3, Neutral Density
	1e19,				// m^-3, Electron Density
	1e19,				// m^-3, Ion Density
	116045.25,			// K, Ion Temperature
	116045.25,			// K, Electron Temperature
	0.025*116045.25,	// K, Neutral Temperature
	300,				// K, Ambient Temperature
	1.66054e-27,		// kg, Mass of ions
	1.0,			// (1/e), Mean Ionization state of plasma
	threevector(),		// m s^-1, Plasma Velocity (Should eventually be normalised to sound speed cs)
	threevector(),		// m s^-2, Acceleration due to gravity
	threevector(),		// V m^-1, Electric field at dust location (Normalised later) 
	threevector(),		// T, Magnetic field at dust location (Normalised later)
};

static struct PlasmaGrid_Data PlasmaGrid_DataDefaults = {
	// Plasma Parameters
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<int> >	(),

	251,
	401,
	0,
	1.5,
	4.0,
	-2.0,
	2.0,

	1.66054e-27,
	'j',
};

class Model{
//	Could also be
 	private:
		// Feasably this could actually be const
		std::shared_ptr<PlasmaGrid_Data> PG_data;		// Holds information about what the current background plasma is.
		int i;								// x Position
		int k;								// y Position

		// Functions used accross different models (Heating, Charging, Forces etc.)
		const double DTOKSIonFlux(double DustTemperature)const;
		const double OMLIonFlux(double DustTemperature)const;
		const double DTOKSElectronFlux(double DustTemperature)const;
		const double OMLElectronFlux(double DustTemperature)const;

		const bool locate(int &i, int &k, threevector xd)const; // Locate dust particle in the plasma grid		
		const bool checkingrid(int i, int k)const;

	protected:
		// It would be nice to have these as const pointers but this is not possible as this requires all function calls
		// Be const ones that do not modify member data. The other option is to capture by reference using Lambda functions.
		// Hah, good luck with that.
		Matter *Sample;						// Tungsten, Beryllium, Iron or Graphite
		std::shared_ptr<PlasmaData> Pdata;	// Const Plasma Data structure
		const float Accuracy;				// Accuracy of model, normalised to one.
		const bool ContinuousPlasma;		// Is the plasma background is constant in space.
		double TimeStep;					// Current time step (s)
		double TotalTime;					// Total Time taken by model (s)

		std::ofstream ModelDataFile;		// Output data file for model
		std::ofstream PlasmaDataFile;		// Plasma data file

		virtual void Print()=0;			// Write to output data file
		virtual double UpdateTimeStep()=0;	// Update Time Scale of development, returns the time step
		virtual double ProbeTimeStep()const=0;	// Check Time Scale of development, returns the time step

		const double IonFlux(double DustTemperature)const;
		const double ElectronFlux(double DustTemperature)const;
		const double NeutralFlux()const;

	public:
		// Constructors
		Model();
		Model(Matter *& sample, PlasmaData &pdata, float accuracy);
		Model(Matter *& sample, PlasmaGrid_Data &pgrid, float accuracy);
		Model(Matter *& sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata, float accuracy );

		// Destructor
		virtual ~Model(){};

		const Matter *get_sample			()const{ return Sample;		}
		const threevector  get_bfield		()const{ return Pdata->MagneticField.getunit(); 		}
		double get_totaltime				()const{ return TotalTime; 	}
		double get_timestep					()const{ return TimeStep; 	}	

		bool new_cell						()const;

		void RecordPlasmadata(std::string filename);				// Record the plasma Data	

		void update_plasmadata(PlasmaData &pdata);
		const bool update_plasmadata();
		void update_fields(int i, int k);

		void close_file();

		void AddTime(double T){	TotalTime = TotalTime + T;	}

		const double get_dlx		()		const{return PG_data->dlx;}

		// Functions for printing
		void vtkcircle(double r, std::ofstream &fout); // Print the inside and the outside of the tokamak 
		void vtktoroid();
		void impurityprint(double totalmass);
		void datadump();

};

#endif
