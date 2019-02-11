/** @class Tungsten.h
 *  @brief Class defining the elemental properties of tungsten
 *  
 *  A class which defines the elemental properties of tungsten relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */


#ifndef __TUNGSTEN_H_INCLUDED__
#define __TUNGSTEN_H_INCLUDED__

#include "Matter.h"

/** @class Tungsten
 *  @brief Class defining the elemental properties of tungsten
 *  
 *  A class which defines the elemental properties of tungsten relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */
class Tungsten: public Matter{

    private:
        void update_radius        ()override;
        void update_heatcapacity  ()override;
        void update_vapourpressure()override;
        void set_defaults         ()override;
        
    public:
        Tungsten();
        Tungsten(double radius);
        Tungsten(double radius, double tempin);
        Tungsten(double radius, double tempin, 
            std::array<char,CM> &constmodels);
        Tungsten(double radius, double tempin, std::array<char,CM> &constmodels,
            const threevector &position, const threevector &velocity);

        ~Tungsten(){};
        
        double probe_vapourpressure(double Temperature)const;
};

#endif /* __TUNGSTEN_H_INCLUDED__ */
