/** @file HeatingModel.h
 *  @brief Contains a class defining the physics model for heating dust
 *  
 *  Class defining the heating models which affect the temperature of a 
 *  dust grain in a magnetically confined fusion tokamak relevant plasma
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug Heating Model number should be in a namespace
 *  @bug FileName needs to be moved into Model.
 *  @bug investigate if UseModel can be moved to Model class
 */

#ifndef __HEATINGMODEL_H_INCLUDED__
#define __HEATINGMODEL_H_INCLUDED__

#include "Model.h"

// HEATING MODEL NUMBER, the number of charge models
const unsigned int HMN = 18;

/** @class Model
 *  @brief Defines the heating models which affect the temperature of dust
 *  
 *  
 */
class HeatingModel : public Model{

    private:
        /** @name Private member data
         *  @brief Information required by multiple functions
         *
         *  \p OldTemp and \p ThermalEquilibrium are used to determine the 
         *  thermal equilibrium condition while RE and RN are intensive 
         *  quantities to calculate used by multiple functions and so are only 
         *  calculated once per time step and stored here. PowerIncident defines
         *  the background power present.
         */
        ///@{
        double OldTemp;
        double RE;
        double RN;
        double PowerIncident;
        bool ThermalEquilibrium;
        ///@}

        /** @brief The filename to which the results are written to
         */
        std::string FileName;

        /** @brief Array defining the heating models being used
         */
        std::array<bool,HMN> UseModel;

        /** @brief Print model data to ModelDataFile
         */
        void Print();

        /** @brief Set default settings for private member data
         */
        void Defaults();

        
        /** @brief Implement RK4 method to calculate the change in temperature
         *  @param timestep is the time over which heating models are active
         *  @return Total energy transfer to dust in \p timestep kJ
         */
        double RungeKutta4(double timestep);

        /** @brief flux of particles away from liquid surface due to evaporation
         */
        const double EvaporationFlux            (double DustTemperature)const;
        /** @name Heating Terms
         *  @brief functions for different terms in heating equation
         *
         *  These functions all return the power due to a particular heat term
         *  in the equation of heating for a dust grain in a tokamak plasma.
         *  The number of functions is equal to HMN+1
         *  @return double the power to the dust grain in kW
         *  @param DustTemperature the temperature of the dust, needed for RK4
         */
        ///@{
        /** @brief Power to surface due to black body radiation
         */
        const double EmissivityModel            (double DustTemperature)const;
        /** @brief Power to surface due evaporative loss of particles
         */
        const double EvaporationModel           (double DustTemperature)const;
        /** @brief Power to surface due to contact with air
         */
        const double NewtonCooling              (double DustTemperature)const;
        /** @brief Power to surface due to neutral bombardment
         */
        const double NeutralHeatFlux            (double DustTemperature)const;
        /** @brief Power to surface due to SOML ion bombardment 
         */
        const double SOMLIonHeatFlux            (double DustTemperature)const;
        /** @brief Power to surface due to SOML neutral recombination
         */
        const double SOMLNeutralRecombination   (double DustTemperature)const;
        /** @brief Power to surface due to SMOML ion bombardment 
         */
        const double SMOMLIonHeatFlux           (double DustTemperature)const;
        /** @brief Power to surface due to SMOML neutral recombination
         */
        const double SMOMLNeutralRecombination  (double DustTemperature)const;

        /** @brief Power to surface due to secondary electron emission
         */
        const double SEE                        (double DustTemperature)const;
        /** @brief Power to surface due to thermionic electron emission
         */
        const double TEE                        (double DustTemperature)const;
        /** @brief Power to surface due to PHL electron bombardment 
         */
        const double PHLElectronHeatFlux        (double DustTemperature)const;
        /** @brief Power to surface due to OML electron bombardment 
         */
        const double OMLElectronHeatFlux        (double DustTemperature)const;
        /** @brief Power to surface due to DTOKS secondary electron emission
         */
        const double DTOKSSEE                   (double DustTemperature)const;
        /** @brief Power to surface due to DTOKS secondary electron emission
         */
        const double DTOKSTEE                   (double DustTemperature)const;
        /** @brief Power to surface due to DTOKS ion bombardment 
         */
        const double DTOKSIonHeatFlux           (double DustTemperature)const;
        /** @brief Power to surface due to DTOKS neutral recombination
         */
        const double DTOKSNeutralRecombination  (double DustTemperature)const;
        /** @brief Power to surface due to DTOKS electron bombardment 
         */
        const double DTOKSElectronHeatFlux      (double DustTemperature)const;
        /** @brief Power to surface due to DUSTT ion bombardment 
         */
        const double DUSTTIonHeatFlux           (double DustTemperature)const;
        /** @brief Calculate the sum of all the heating models
         */
        double CalculatePower(double DustTemperature)const;
        ///@}

    public:

        HeatingModel();
        HeatingModel( std::string filename, float accuracy, 
            std::array<bool,HMN> &models, Matter *& sample, PlasmaData & pdata);
        HeatingModel( std::string filename, float accuracy, 
            std::array<bool,HMN> &models, Matter *& sample, PlasmaData * pdata);
        HeatingModel( std::string filename, float accuracy, 
            std::array<bool,HMN> &models, Matter *& sample, 
            PlasmaGrid_Data & pgrid);
        HeatingModel( std::string filename, float accuracy, 
            std::array<bool,HMN> &models, Matter *& sample, 
            PlasmaGrid_Data & pgrid, PlasmaData &pdata);

        ~HeatingModel(){
        };


        /** @brief Continuously heat the sample in steps of size \p TimeStep
         *  @return 1 for boiled, 2 for evaporated and 3 for thermal equilibrium
         */
        const int Vapourise();

        /** @brief Heat the sample for a time period of \p TimeStep
         */
        void Heat();

        /** @brief Heat the sample for a time period of \p timestep
         *  @param timestep the time period for which the sample is heated
         */
        void Heat(double timestep);

        /** @brief Update the values of RE and RN using external functions
         */
        void UpdateRERN();

        
        /** @brief Create file and by default don't print phase data
         */
        void CreateFile(std::string filename, bool PrintPhaseData);
        void CreateFile(std::string filename);
        double ProbeTimeStep()const;
        double UpdateTimeStep();

        /** @name Public getter methods
         *  @brief functions required to get member data
         *  
         *  Minimal getter methods providing functionality required for various
         *  tests and mechanisms implemented by DTOKSU.
         */
        bool get_thermalequilibrium     (){ return ThermalEquilibrium; }

        /** @name Public setter methods
         *  @brief functions required to set member data
         *  
         *  Minimal setter methods providing functionality required for various
         *  tests and mechanisms implemented by DTOKSU.
         */
        ///@{
        void set_PowerIncident          (double powerincident)
        { 
            PowerIncident = powerincident;      
        }
        void set_ThermalEquilibrium     (double thermequilib )
        { 
            ThermalEquilibrium = thermequilib;  
        }
        ///@}
};

#endif /* __HEATINGMODEL_H_INCLUDED__ */

