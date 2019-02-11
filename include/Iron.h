/** @class Iron.h
 *  @brief Class defining the elemental properties of iron
 *  
 *  A class which defines the elemental properties of iron relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */

#ifndef __IRON_H_INCLUDED__
#define __IRON_H_INCLUDED__

#include "Matter.h"

/** @class Iron
 *  @brief Class defining the elemental properties of iron
 *  
 *  A class which defines the elemental properties of iron relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */
class Iron : public Matter{

    private:
        void update_radius        ()override;
        void update_heatcapacity  ()override;
        void update_vapourpressure()override;
        void set_defaults         ()override;
        

    public:
        Iron();
        Iron(double radius);
        Iron(double radius, double tempin);
        Iron(double radius, double tempin, std::array<char,CM> &constmodels);
        Iron(double radius, double tempin, std::array<char,CM> &constmodels,
                const threevector &position, const threevector &velocity);

        ~Iron(){};
        
        double probe_vapourpressure(double Temperature)const;
};

#endif /* __IRON_H_INCLUDED__ */
