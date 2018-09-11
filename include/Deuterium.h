#ifndef __DEUTERIUM_H_INCLUDED__   // if Iron.h hasn't been included yet...
#define __DEUTERIUM_H_INCLUDED__

#include <iostream>
#include "Matter.h"

class Deuterium: public Matter{

	private:

		// Functions called by Deuterium::update()
		void update_radius		();
		void update_heatcapacity 	();
		void update_vapourpressure	();

	public:
		// Constructors
		Deuterium();
		Deuterium(double radius);
		Deuterium(double radius, double tempin);
		Deuterium(double radius, double tempin, std::array<char,CM> &constmodels);
		Deuterium(double radius, double tempin, std::array<char,CM> &constmodels,
			const threevector& position, const threevector& velocity);

		// Destructor
		~Deuterium(){};

		
		// Change Properties; Mass and Temperature
		void set_defaults		();
};

#endif
