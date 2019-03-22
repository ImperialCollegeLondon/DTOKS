/** @file Matter.h
 *  @brief Class defining the matter object for use in physics models
 *  
 *  An abstract base class dealing with the generic physical behaviour of matter.
 *  Two structures defined by \p St and \p Ec define all the physical behaviour.
 *  The vector \p ConstModels defines the const-ness of the heat capacity, 
 *  emissivity and thermal expansion as implemented by child classes.
 *  This class is inherited by specific elemental classes which specify the 
 *  elemental constants of \p Ec and specify the dependency of the constants on
 *  other parameters.
 *  
 *
 *  Assumptions:
 *  The mass is assumed to be perfectly conducting, with neglidgeable fluid 
 *  effects in the liquid phase such that it's internal dynamics can be treated 
 *  as a rigid body. Rotational dynamics are neglected except in a single 
 *  dimension for a niche application. The matter is perfectly spherical.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 *  @bug MinMass should exist within it's own namespace
 *  @bug PreBoilMas, can it be removed?
 *  @bug ConstModels, change it to an enum?
 *  @bug ConstModels, better documentation
 */

#ifndef __MATTER_H_INCLUDED__
#define __MATTER_H_INCLUDED__

//!< The lower bound on mass, below which simulations terminate
#define MinMass 10e-25 

#include <iostream>       //!< I/O operations and debugging
#include <fstream>        //!< Open emissivity filestream
#include <sstream>        //!< Convert temp to file number
#include <assert.h>       //!< Assertion errors
#include <math.h>         //!< Round
#include <array>          //!< std::array
#include <limits>         //!< std::numeric_limits<double>::max()
#include <string.h>       //!< strchr("",std::string)

#include "GrainStructs.h" //!< Contains the structures for material properties
#include "Constants.h"    //!< Contains general physical constants
#include "Functions.h"    //!< sec(Te,'f') function used by HeatingModel.cpp

//!< Constant model number, the number of constant models
const unsigned int CM = 5;

/** @class Matter
 *  @brief Class defining the matter object for use in physics models
 *  
 *  Pure virtual base class defining the functionality of matter relevant to 
 *  simulating dust grains in plasmas.
 */
class Matter{

    private:
        
        /** @name Private Functions
         *  @brief functions to update material properties
         *
         * Update geometric properties of matter:
         * Volume, Surface area, Density, emissiviy and boiling temperature
         */
        ///@{
        /** @brief Update information relevant to the physical extent.
         *
         *  Update the volume, surface area, mass, density and radius taking 
         *  account for the constness
         */
        void update_dim();
        /** @brief Update emissivity of material from data tables.
         *
         *  Emissivity is altered as a function of temperature and radius if the
         *  data tables are available for the element
         */
        void update_emissivity();
        /** @brief Update emissivity of material from data tables.
         *
         *  Four different models for boiling termed (y)es, (n)o, (s)uper and 
         *  (t)homson which produce different behaviour for the boiling 
         *  temperature
         */
        void update_boilingtemp();
        ///@}

        //<! Record the mass before boiling
        double PreBoilMass;              
        //<! Constant Models variation with Temperature turned on of possibly CM
        std::array<char,CM> ConstModels;
        
    protected: //<! Functions used by the elements inheriting Matter.

        /** @name Constructors
         *  @brief functions to construct Matter class
         */
         ///@{
        /** @brief Single Parameter constructor.
         *
         *  Assume all models constant, set element constants and use 
         *  MatterDefaults.
         *  @param elementconsts data structure defining element properties
         */
        Matter(const ElementConsts &elementconsts);

        /** @brief Two Parameter constructor.
         *
         *  Assume all models constant, set element constants and use 
         *  MatterDefaults. Set radius to \p rad
         *  @param elementconsts data structure defining element properties
         *  @param rad defines the radius of the sphere
         */
        Matter(double rad, const ElementConsts &elementconsts);

        /** @brief Three Parameter constructor.
         *
         *  Assume all models constant, set element constants and use 
         *  MatterDefaults. Set radius to \p rad and temperature to \p temp.
         *  @param elementconsts data structure defining element properties
         *  @param rad defines the radius of the sphere
         *  @param temp defines the initial temperature of the sphere
         */
        Matter(double rad, double temp, const ElementConsts &elementconsts);

        /** @brief Five Parameter constructor.
         *
         *  Set element constants and use MatterDefaults. Set radius to \p rad,
         *  temperature to \p temp and the variability of models to 
         *  \p constmodels
         *  @param elementconsts data structure defining element properties
         *  @param rad defines the radius of the sphere
         *  @param temp defines the initial temperature of the sphere
         *  @param constmodels defines the variability of models
         */
        Matter(double rad, double temp, const ElementConsts &elementconsts, 
            std::array <char,CM> &constmodels);

        /** @brief Six Parameter constructor.
         *
         *  Set element constants and use MatterDefaults. Set radius to \p rad,
         *  temperature to \p temp, the variability of models to 
         *  \p constmodels and the position and velocity to \p Position and 
         *  \p Velocity.
         *  @param elementconsts data structure defining element properties
         *  @param rad defines the radius of the sphere
         *  @param temp defines the initial temperature of the sphere
         *  @param constmodels defines the variability of models
         *  @param Position three dimensional vector giving initial position
         *  @param Velocity three dimensional vector giving initial velocity
         */
        Matter(double rad, double temp, const ElementConsts &elementconsts, 
            std::array <char,CM> &constmodels, const threevector &Position, 
            const threevector &Velocity);
        ///@}

        /** @name Member Data
         *  @brief Protected data defining the element
         *  
         *  The GrainData structure contains non-constant information about the
         *  matter which varies during the simulation.
         *  The ElementConsts structure contains constant information about the
         *  matter which is instantiated on construction of Matter by the 
         *  derived classes.
         *  ConstModels is a character array of length \p CM which defines the 
         *  variablility of the Emissivity, thermal expansion, heat capacity,
         *  boiling and breakup. See configuration file for details.
         */
        ///@{
        struct GrainData           St;
        const struct ElementConsts Ec;
        ///@}

        /** @name Pure virtual function
         *  @brief Functions defining element specific variability of constants
         *  
         *  The GrainData structure contains non-constant information about the
         *  matter which varies during the simulation.
         *  The ElementConsts structure contains constant information about the
         *  matter which is instantiated on construction of Matter by the 
         *  derived classes.
         *  ConstModels is a character array of length \p CM which defines the 
         *  variablility of the Emissivity, thermal expansion, heat capacity,
         *  boiling and breakup. See configuration file for details.
         *  the functions with the prefix update_ called by Matter::update(). 
         *  These are element dependant and must be defined in the child class.
         */
        ///@{
        virtual void set_defaults          ()=0;

        /** @brief Change radius according to thermal expansion
         *  
         *  Mutates the radius in \p St following thermal expansion of material
         */
        virtual void update_radius         ()=0;

        /** @brief Change heat capacity as a function of matter temperature
         *  
         *  Mutates the heat capacity in \p St in units of J mol^-1 K^-2
         */
        virtual void update_heatcapacity   ()=0; 

        /** @brief Change vapour pressure as a function of matter temperature
         *  
         *  Mutates the vapour pressure in \p St in units of pascals.
         *  @see probe_vapourpressure()
         */
        virtual void update_vapourpressure ()=0;
        ///@}

        /** @name Protected setter methods for state and models
         *  @brief Set the correct material state and \p ConstModels
         *  
         *  These funcitons are implemented in the constructors of derived 
         *  classes to ensure that the material is in the correct state and 
         *  operating with the correct values for \p ConstModels 
         */
        ///@{
        /** @brief Update the state of the element
         *  
         *  Change the temperature of the matter and account for changes in
         *  state when reaching either the boiling or melting temperature.
         *  @param EnergyIn amount of energy added to material in J
         */
        void update_state(double EnergyIn);
        /** @brief Update ConstModels
         *  
         *  Update the ConstModels array with new values.
         *  @param emissivmodel new value of \p ConstModels[0]
         *  @param linexpanmodel new value of \p ConstModels[1]
         *  @param heatcapacity new value of \p ConstModels[2]
         *  @param boilingmodel new value of \p ConstModels[3]
         *  @param breakupmodel new value of \p ConstModels[4]
         */
        void update_models(char emissivmodel, char linexpanmodel, 
            char heatcapacity, char boilingmodel, char breakupmodel);
        /** @brief Update ConstModels
         *  
         *  Update the ConstModels array with new values.
         *  @param constmodels new values for \p ConstModels
         */
        void update_models(std::array<char,CM> &constmodels);
        ///@}

    public:

        virtual ~Matter         (){M_Debug("M_Debug was on.\n");};

        /** @name Update functions
         *  @brief Implemented by derived model classes to mutate member data
         *  
         *  The following functions are used by the ForceModel, HeatModel and
         *  ChargeModel classes to muatate the \p St data relevant to the dust
         *  grain. As the physics of temperature change for a given input energy
         *  is fixed, that is implemented here.
         */
        ///@{
        /** @brief update geometric and variable properties of matter
         *  
         *  Call the functions update_dim(), update_heatcapacity(), 
         *  update_emissivity(), update_vapourpressure() and 
         *  update_boilingtemp() in the order. Also check if rotational 
         *  breakup has occured
         *  @see update_dim()
         *  @see update_heatcapacity()
         *  @see update_emissivity()
         *  @see update_vapourpressure()
         *  @see update_boilingtemp()
         */
        void update();
        /** @brief Remove the amount of mass lost from some process
         *  
         *  Implemented by HeatingModel after calculating total evaporation
         *  Change the amount of mass in \p St whilst making appropriate checks
         *  @param LostMass amount of matter lost in kg
         */
        void update_mass(double LostMass);
        /** @brief Change temperature and state of matter
         *  
         *  Implemented by HeatingModel after calculating total heat energy
         *  
         *  @param EnergyIn amount of energy added to material in J
         */
        void update_temperature(double EnergyIn);
        /** @brief Change the position, velocity and rotational velocity
         *  
         *  Implemented by ForceModel after calculating total acceleration
         *  Perform appropriate checks on the input before calling update_state
         *  @see update_state()
         *  @param changeinposition threevector distance moved in m
         *  @param changeinvelocity threevector change in velocity in m/s
         *  @param Rotation increase in rotational velocity
         */
        void update_motion(const threevector &changeinposition, 
            const threevector &changeinvelocity, double Rotation);
        /** @brief Change the normalised potential, charge and deltas
         *  
         *  Implemented by ChargingModel after calculating the charge
         *  @param charge the number of electronic charges on the sphere
         *  @param potential the normalised potential of the sphere
         *  @param deltas Proportion of secondary electron emission current
         *  @param deltats Proportion of thermionic electron emission current
         */
        void update_charge(double charge, double potential, double deltas, 
            double deltat);
        ///@}

        /** @name Public setter methods
         *  @brief functions required to set member data
         *  
         *  Minimal setter methods providing functionality required for various
         *  mechanisms implemented by DTOKSU, principally rotational breakup.
         */
        ///@{
        /** @brief Set the GrainData
         *  
         *  Implemented by DTOKSU_Manager to save the dust grain state before
         *  breakup
         *  @param NewData the new state of the matter being considered
         */
        void set_graindata(GrainData &NewData){ St = NewData; };
        /** @brief Set the potential
         *  
         *  Fix potential to the value \p potential
         *  @param potential force the normalised potential to be \p potential
         */
        void set_potential      (double potential){ St.Potential = potential; };
        /** @brief Set the mass
         *  
         *  Fix mass to the value \p mass
         *  @param mass force the St.Mass to be \p mass
         */
        void set_mass           (double mass)
        { 
            assert(mass > MinMass*10); 
            St.Mass = mass;   
        };
        /** @brief Set the restitution coefficients
         * 
         *  Set the values of \p RE and \p RN in \p St
         *  @param RE the fraction of backscattered energy
         *  @param RN the fraction of backscattered particles
         */
        void set_rern(double re, double rn){ St.RE = re; St.RN = rn; };
        /** @brief Set Breakup to false
         *  
         *  Implemented by DTOKSU_Manager to set breakup flag to false
         */
        void reset_breakup      (){ St.Breakup = false;         };
        ///@}

        /** @name Public getter methods
         *  @brief functions required to get member data
         *  
         *  Minimal getter methods providing functionality required for access 
         *  to various member data.
         */
        ///@{
        char get_elem               ()const{ return Ec.Elem;                };
        double get_meltingtemp      ()const{ return Ec.MeltingTemp;         };
        double get_boilingtemp      ()const{ return Ec.BoilingTemp;         };
        double get_bondenergy       ()const{ return Ec.BondEnergy;          };
        double get_workfunction     ()const{ return Ec.WorkFunction;        };
        double get_heattransair     ()const{ return Ec.HeatTransAir;        };
        double get_atomicmass       ()const{ return Ec.AtomicMass;          };
        double get_surfacetension   ()const{ return Ec.SurfaceTension;      };
        double get_superboilingtemp ()const{ return St.SuperBoilingTemp;    };
        double get_mass             ()const{ return St.Mass;                };
        double get_density          ()const{ return St.Density;             };
        double get_radius           ()const{ return St.Radius;              };
        double get_heatcapacity     ()const{ return St.HeatCapacity;        };
        double get_temperature      ()const{ return St.Temperature;         };
        double get_fusionenergy     ()const{ return St.FusionEnergy;        };
        double get_vapourenergy     ()const{ return St.VapourEnergy;        };
        double get_latentvapour     ()const{ return Ec.LatentVapour;        };
        double get_latentfusion     ()const{ return Ec.LatentFusion;        };
        double get_emissivity       ()const{ return St.Emissivity;          };
        double get_surfacearea      ()const{ return St.SurfaceArea;         };
        double get_vapourpressure   ()const{ return St.VapourPressure;      };
        double get_linearexpansion  ()const{ return St.LinearExpansion;     };
        double get_deltasec         ()const{ return St.DeltaSec;            };
        double get_deltatherm       ()const{ return St.DeltaTherm;          };
        double get_re               ()const{ return St.RE;                  };
        double get_rn               ()const{ return St.RN;                  };
        double get_deltatot         ()const
        { 
            return (St.DeltaTherm + St.DeltaSec);  
        };
        double get_potential        ()const{ return St.Potential;           };
        double get_rotationalfreq   ()const{ return St.RotationalFrequency; };
        threevector get_velocity    ()const{ return St.DustVelocity;        };
        threevector get_position    ()const{ return St.DustPosition;        };
        GrainData get_graindata     ()const{ return St;                     };
        bool is_gas                 ()const{ return St.Gas;                 };
        bool is_liquid              ()const{ return St.Liquid;              };
        bool is_split               ()const{ return St.Breakup;             };
        bool is_positive            ()const{ return St.Positive;            };
        char get_c                  (int i)const
        { 
            assert(i < CM); 
            return ConstModels[i];   
        };
        /** @brief Get the vapour pressure without mutating member data
         *  
         *  This is used by the force model to measure vapour pressure for 
         *  determining the rocket force. 
         */
        virtual double probe_vapourpressure (double Temperature)const=0;
        ///@}
};

#endif /* __MATTER_H_INCLUDED__ */

