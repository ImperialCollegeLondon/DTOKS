#ifndef __MODEL_H_INCLUDED__   // if Model.h hasn't been included yet...
#define __MODEL_H_INCLUDED__

#include <memory>

#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

const struct PlasmaData PlasmaDefaults = {
	1e20,		// m^-3, Neutral Density
	1e20,		// m^-3, Electron Density
	1e20,		// m^-3, Electron Density
	116045.25,	// K, Ion Temperature
	116045.25,	// K, Electron Temperature
	116045.25,	// K, Neutral Temperature
	300,		// K, Ambient Temperature
	threevector(),	// m s^-1, Plasma Velocity (Should eventually be normalised to sound speed cs)
	threevector(),	// V m^-1, Electric field at dust location (Normalised later) 
	threevector(),	// T, Magnetic field at dust location (Normalised later)
};

class Model{
//	Could also be
// 	Private:
// 	This requires changing access methods in derived classes
	protected:
		// Parameters defining the Heating equation
		Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite
//		std::shared_ptr <Matter> Sample;	// Tungsten, Beryllium, Iron or Graphite
		PlasmaData const &Pdata;		// Reference to Plasma Data structure
		
		std::ofstream ModelDataFile;		// Output data file

	protected:
		
		double Accuracy;			// Accuracy of model
		virtual void Print()=0;			// Write to output data file
		virtual double CheckTimeStep()=0;	// Check Time Scale of development, returns the time step
//		void Reset_Data( std::shared_ptr <Matter> const& sample, PlasmaData const& pdata);

	public:
		// Constructors
		Model();
		Model(Matter *& sample, PlasmaData const& pdata, double accuracy);

		// Destructor
		virtual ~Model(){};

//		std::shared_ptr<Matter> get_sample() 	{ return Sample;	}
		Matter * get_sample() 			{ return Sample;	}

};

#endif
