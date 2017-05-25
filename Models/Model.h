#ifndef __MODEL_H_INCLUDED__   // if Model.h hasn't been included yet...
#define __MODEL_H_INCLUDED__

#include <memory>

#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"

class Model{
//	Could also be
// 	Private:
// 	This requires changing access methods in derived classes
	protected:
		// Parameters defining the Heating equation
		//Matter *Sample;				// Tungsten, Beryllium, Iron or Graphite
		std::shared_ptr <Matter> Sample;		// Tungsten, Beryllium, Iron or Graphite
		PlasmaData const &Pdata;			// Reference to Plasma Data structure
		
		std::ofstream ModelDataFile;		// Output data file

	protected:

		virtual void Print()=0;			// Write to output data file
//		void Reset_Data( std::shared_ptr <Matter> const& sample, PlasmaData const& pdata);

	public:
		// Constructors
		Model();
		Model(std::shared_ptr <Matter> const& sample, PlasmaData const& pdata);

		// Destructor
		virtual ~Model(){};

		std::shared_ptr<Matter> get_sample() 	{ return Sample;	}
};

#endif
