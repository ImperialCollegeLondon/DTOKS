//#include <math.h>
#ifndef __IRON_H_INCLUDED__   // if Iron.h hasn't been included yet...
#define __IRON_H_INCLUDED__

#include <iostream>
#include "Matter.h"

class Iron : public Matter{

	private:
		
		// Functions called by Iron::update()
		void update_radius	 	();
		void update_heatcapacity 	();
		void update_vapourpressure	();
		

	public:
		// Constructors
		Iron();
		Iron(double radius);
		Iron(double radius, double tempin);
		Iron(double radius, double tempin, std::array<char,CM> &constmodels);
		Iron(double radius, double tempin, std::array<char,CM> &constmodels,
				const threevector &position, const threevector &velocity);

		// Destructor
		~Iron(){};
		
		// Change Properties; Mass and Temperature
		void set_defaults		();
		double probe_vapourpressure	(double Temperature)const;
};

#endif
