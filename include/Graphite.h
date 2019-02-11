/** @class Graphite.h
 *  @brief Class defining the elemental properties of graphite
 *  
 *  A class which defines the elemental properties of graphite relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */

#ifndef __GRAPHITE_H_INCLUDED__
#define __GRAPHITE_H_INCLUDED__

#include "Matter.h"

/** @class Graphite
 *  @brief Class defining the elemental properties of graphite
 *  
 *  A class which defines the elemental properties of graphite relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */
class Graphite: public Matter{

    private:
        void update_radius        ()override;
        void update_heatcapacity  ()override;
        void update_vapourpressure()override;
        void set_defaults         ()override;
        
    public:
        Graphite();
        Graphite(double radius);
        Graphite(double radius, double tempin);
        Graphite(double radius, double tempin, 
            std::array<char,CM> &constmodels);
        Graphite(double radius, double tempin, std::array<char,CM> &constmodels,
            const threevector &position, const threevector &velocity);

        ~Graphite(){};
        
        double probe_vapourpressure(double Temperature)const;
};

#endif /* __GRAPHITE_H_INCLUDED__ */
