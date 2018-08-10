#ifndef __BERYLLIUM_H_INCLUDED__   // if Iron.h hasn't been included yet...
#define __BERYLLIUM_H_INCLUDED__

#include <iostream>
#include "Matter.h"

class Beryllium: public Matter{

	private:

		// Functions called by Beryllium::update()
		void update_radius		();
		void update_heatcapacity 	();
		void update_vapourpressure	();

	public:
		// Constructors
		Beryllium();
		Beryllium(double radius);
		Beryllium(double radius, double tempin);
		Beryllium(double radius, double tempin, std::array<char,4> &constmodels);
		Beryllium(double radius, double tempin, std::array<char,4> &constmodels,
			const threevector& position, const threevector& velocity);

		// Destructor
		~Beryllium(){};

		
		// Change Properties; Mass and Temperature
		void set_defaults		();
};

#endif
