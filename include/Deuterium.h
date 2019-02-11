/** @class Deuterium.h
 *  @brief Class defining the elemental properties of deuterium
 *  
 *  A class which defines the elemental properties of deuterium relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */

#ifndef __DEUTERIUM_H_INCLUDED__
#define __DEUTERIUM_H_INCLUDED__

#include "Matter.h"

/** @class Deuterium
 *  @brief Class defining the elemental properties of deuterium
 *  
 *  A class which defines the elemental properties of deuterium relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */
class Deuterium: public Matter{

    private:
        void update_radius        ()override;
        void update_heatcapacity  ()override;
        void update_vapourpressure()override;
        void set_defaults         ()override;
        
    public:
        Deuterium();
        Deuterium(double radius);
        Deuterium(double radius, double tempin);
        Deuterium(double radius, double tempin, 
            std::array<char,CM> &constmodels);
        Deuterium(double radius, double tempin, 
            std::array<char,CM> &constmodels, const threevector& position, 
            const threevector& velocity);

        ~Deuterium(){};

        double probe_vapourpressure(double Temperature)const;
};

#endif /* __DEUTERIUM_H_INCLUDED__ */