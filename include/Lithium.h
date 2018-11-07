#ifndef __LITHIUM_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __LITHIUM_H_INCLUDED__

//#include <iostream>
#include "Matter.h"

class Lithium: public Matter{

	private:

		// Functions called by Lithium::update()
		void update_radius		();
		void update_heatcapacity 	();
		void update_vapourpressure	();
		
	public:
		// Constructors
		Lithium();
		Lithium(double radius);
		Lithium(double radius, double tempin);
		Lithium(double radius, double tempin, std::array<char,CM> &constmodels);
		Lithium(double radius, double tempin, std::array<char,CM> &constmodels,
			const threevector &position, const threevector &velocity);

		// Destructor
		~Lithium(){};
		
		// Change Properties; Mass and Temperature
		void lithium_defaults		();
		double probe_vapourpressure(double Temperature)const;
};

#endif
