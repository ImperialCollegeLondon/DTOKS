#ifndef __MOLYBDENUM_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __MOLYBDENUM_H_INCLUDED__

//#include <iostream>
#include "Matter.h"

class Molybdenum: public Matter{

	private:

		// Functions called by Molybdenum::update()
		void update_radius		();
		void update_heatcapacity 	();
		void update_vapourpressure	();
		
	public:
		// Constructors
		Molybdenum();
		Molybdenum(double radius);
		Molybdenum(double radius, double tempin);
		Molybdenum(double radius, double tempin, std::array<char,CM> &constmodels);
		Molybdenum(double radius, double tempin, std::array<char,CM> &constmodels,
			const threevector &position, const threevector &velocity);

		// Destructor
		~Molybdenum(){};
		
		// Change Properties; Mass and Temperature
		void molybdenum_defaults		();
		double probe_vapourpressure(double Temperature)const;
};

#endif
