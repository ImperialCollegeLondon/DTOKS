/** @class Lithium.h
 *  @brief Class defining the elemental properties of lithium
 *  
 *  A class which defines the elemental properties of lithium relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */

#ifndef __LITHIUM_H_INCLUDED__
#define __LITHIUM_H_INCLUDED__

#include "Matter.h"

/** @class Iron
 *  @brief Class defining the elemental properties of lithium
 *  
 *  A class which defines the elemental properties of lithium relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */
class Lithium: public Matter{

    private:
        void update_radius        ()override;
        void update_heatcapacity  ()override;
        void update_vapourpressure()override;
        void set_defaults         ()override;
        
    public:
        Lithium();
        Lithium(double radius);
        Lithium(double radius, double tempin);
        Lithium(double radius, double tempin, std::array<char,CM> &constmodels);
        Lithium(double radius, double tempin, std::array<char,CM> &constmodels,
            const threevector &position, const threevector &velocity);

        ~Lithium(){};
        
        double probe_vapourpressure(double Temperature)const;
};

#endif /* __LITHIUM_H_INCLUDED__ */
