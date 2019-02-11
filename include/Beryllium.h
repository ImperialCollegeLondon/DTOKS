/** @file Beryllium.h
 *  @brief Class defining the elemental properties of beryllium
 *  
 *  Contains a class which defines the elemental properties of beryllium 
 *  relevant to the simulation of dust in plasmas. An instance of this class is
 *  used to represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#ifndef __BERYLLIUM_H_INCLUDED__
#define __BERYLLIUM_H_INCLUDED__

#include "Matter.h"

/** @class Beryllium
 *  @brief Class defining the elemental properties of beryllium
 *  
 *  A class which defines the elemental properties of beryllium relevant to the 
 *  simulation of dust in plasmas. An instance of this class is used to 
 *  represent a spherical mass of this element, with data structures and
 *  functionality derived from the Matter class
 *  @see Matter
 */
class Beryllium: public Matter{

    private:
        void update_radius        ()override;
        void update_heatcapacity  ()override;
        void update_vapourpressure()override;
        void set_defaults         ()override;

    public:
        Beryllium();
        Beryllium(double radius);
        Beryllium(double radius, double tempin);
        Beryllium(double radius, double tempin,
            std::array<char,CM> &constmodels);
        Beryllium(double radius, double tempin, 
            std::array<char,CM> &constmodels, const threevector& position, 
            const threevector& velocity);

        ~Beryllium(){};

        double probe_vapourpressure(double Temperature)const;
};

#endif /* __BERYLLIUM_H_INCLUDED__ */
