/** @class Molybdenum.h
 *  @brief Class defining the elemental properties of molybdenum
 *  
 *  A class which defines the elemental properties of molybdenum relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */

#ifndef __MOLYBDENUM_H_INCLUDED__
#define __MOLYBDENUM_H_INCLUDED__

#include "Matter.h"

/** @class Molybdenum
 *  @brief Class defining the elemental properties of molybdenum
 *  
 *  A class which defines the elemental properties of molybdenum relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */
class Molybdenum: public Matter{

    private:
        void update_radius        ()override;
        void update_heatcapacity  ()override;
        void update_vapourpressure()override;
        void set_defaults         ()override;
        
    public:
        Molybdenum();
        Molybdenum(double radius);
        Molybdenum(double radius, double tempin);
        Molybdenum(double radius, double tempin, 
            std::array<char,CM> &constmodels);
        Molybdenum(double radius, double tempin, 
            std::array<char,CM> &constmodels, const threevector &position, 
            const threevector &velocity);

        ~Molybdenum(){};
        
        double probe_vapourpressure(double Temperature)const;
};

#endif /* __MOLYBDENUM_H_INCLUDED__ */
