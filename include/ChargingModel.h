/** @file ChargingModel.h
 *  @brief Contains a class defining the physics model for charging dust
 *  
 *  Class defining the charging models which affect the potential of a 
 *  dust grain in a magnetically confined fusion tokamak relevant plasma
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __CHARGINGMODEL_H_INCLUDED__
#define __CHARGINGMODEL_H_INCLUDED__

#include "Model.h"
#include "solveMOMLEM.h"
#include "ChargingTerms.h"

/** @class Model
 *  @brief Defines the charging models which affect the potential of dust
 *  
 *  
 */
class ChargingModel : public Model{

    private:
        /** @brief Array defining the charging models being used
         */
        std::vector<CurrentTerm*> CurrentTerms; 

        /** @brief Print model data to ModelDataFile
         */
        void Print();


    public:
        ChargingModel();
        ChargingModel(std::string filename, float accuracy, 
            std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, PlasmaData & pdata);
        ChargingModel(std::string filename, float accuracy, 
            std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, PlasmaData * pdata);
        ChargingModel(std::string filename, float accuracy, 
            std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, 
            PlasmaGrid_Data & pgrid);
        ChargingModel(std::string filename, float accuracy, 
            std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, 
            PlasmaGrid_Data & pgrid, PlasmaData &pdata);


        ~ChargingModel(){
        };

        void CreateFile(std::string filename);
        double ProbeTimeStep()const;
        double UpdateTimeStep();
        
        /** @brief Charge the sample for a time period of \p TimeStep
         */
        void Charge();
        /** @brief Charge the sample for a time period of \p timestep
         *  @param timestep the time period for which the sample is charged
         */
        void Charge(double timestep);
};

#endif /* __CHARGINGMODEL_H_INCLUDED__ */
