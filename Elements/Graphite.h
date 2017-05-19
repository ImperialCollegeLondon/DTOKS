//#include <math.h>
#ifndef __GRAPHITE_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __GRAPHITE_H_INCLUDED__

#include <iostream>
#include "Matter.h"

class Graphite: public Matter{

	private:

		static const struct ElementConsts GraphiteConsts;

		// Functions called by Matter::update()
		void update_radius		();
		void update_heatcapacity 	();
		void update_vapourpressure	();

	public:
		// Constructors
		Graphite();
		Graphite(double radius);
		Graphite(double radius, double tempin);
		Graphite(double radius, double tempin, std::array<char,4> &constmodels);

		// Destructor
		~Graphite(){};
		
		// Change Properties; Mass and Temperature
		void set_defaults		();
};

#endif
