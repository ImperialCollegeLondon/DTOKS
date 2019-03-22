/** @file HeatingModel.h
 *  @brief Contains a class defining the physics model for heating dust
 *  
 *  Class defining the heating models which affect the temperature of a 
 *  dust grain in a magnetically confined fusion tokamak relevant plasma
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __HEATINGMODEL_H_INCLUDED__
#define __HEATINGMODEL_H_INCLUDED__

#include "Model.h"
#include "HeatTerms.h"

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
         *  thermal equilibrium condition. \p PowerIncident defines the
         *  background power present.
         */
        ///@{
        double OldTemp;
        double PowerIncident;
        bool ThermalEquilibrium;
        ///@}

        /** @brief vector of heating terms defining the heating terms used
         */
        std::vector<HeatTerm*> HeatTerms;

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

        /** @brief Calculate the sum of all the heating models
         */
        double CalculatePower(double DustTemperature)const;

    public:

        HeatingModel();
        HeatingModel( std::string filename, float accuracy, 
            std::vector<HeatTerm*> heatterms, Matter *& sample, 
            PlasmaData & pdata);
        HeatingModel( std::string filename, float accuracy, 
            std::vector<HeatTerm*> heatterms, Matter *& sample, 
            PlasmaData * pdata);
        HeatingModel( std::string filename, float accuracy, 
            std::vector<HeatTerm*> heatterms, Matter *& sample, 
            PlasmaGrid_Data & pgrid);
        HeatingModel( std::string filename, float accuracy, 
            std::vector<HeatTerm*> heatterms, Matter *& sample, 
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

