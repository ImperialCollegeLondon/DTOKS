#ifndef __MODEL_H_INCLUDED__   // if Model.h hasn't been included yet...
#define __MODEL_H_INCLUDED__

#include <memory>

#include "PlasmaGrid.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

class Model{
//	Could also be
 	private:
		// Feasably this could actually be const
		PlasmaGrid *Pgrid;			// Holds information about what the current background plasma is.

// 	This next section could also be private but this requires changing access methods in derived classes.
//	This should be done at some later stage to ensure security in protection of Matter* Sample and Pdata.
	protected:
		// Parameters defining the Heating equation
		// It would be nice to have these as const pointers but this is not possible as this requires all function calls
		// Be const ones that do not modify member data. The other option is to capture by reference using Lambda functions.
		// Hah, good luck with that.
		Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite
		PlasmaData *Pdata;			// Const Plasma Data structure
		std::ofstream ModelDataFile;		// Output data file
		const double Accuracy;			// Accuracy of model, normalised to one.
		const bool ContinuousPlasma;		// Is the plasma background is constant in space.
		double TimeStep;			// Current time step (s)
		double TotalTime;			// Total Time taken by model (s)

	protected:
		virtual void Print()=0;			// Write to output data file
		virtual double UpdateTimeStep()=0;	// Update Time Scale of development, returns the time step
		virtual double ProbeTimeStep()const=0;	// Check Time Scale of development, returns the time step

	public:
		// Constructors
		Model();
		Model(Matter *& sample, PlasmaData *&pdata, double accuracy);
		Model(Matter *& sample, PlasmaGrid &pgrid, double accuracy);

		// Destructor
		virtual ~Model(){};

		const Matter *get_sample		()const{ return Sample;		}
		const PlasmaData *get_plasmadata	()const{ return Pdata; 		}
		double get_totaltime			()const{ return TotalTime; 	}
		double get_timestep			()const{ return TimeStep; 	}

		void update_plasmadata(PlasmaData *&pdata);
		bool update_plasmadata(threevector pos);
		void update_fields(int i, int k);

		void AddTime(double T){	TotalTime = TotalTime + T;	}
};

#endif
