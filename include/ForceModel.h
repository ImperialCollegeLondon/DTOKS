/** @file ForceModel.h
 *  @brief Contains a class defining the physics model for equation of motion
 *  
 *  Class defining the force models which affect the equation of motion of a 
 *  dust grain in a magnetically confined fusion tokamak relevant plasma
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __FORCEMODEL_H_INCLUDED__
#define __FORCEMODEL_H_INCLUDED__

#include "Model.h"
#include "ForceTerms.h"

/** @class ForceModel
 *  @brief Defines the force models which affect the equation of motion of dust
 *  
 *  
 */
class ForceModel : public Model {

    private:
        /** @brief vector of force terms defining the force terms used
         */
        std::vector<ForceTerm*> ForceTerms; 

        /** @brief The temperature of the dust grain in the previous step
         */
        double OldTemp;

        /** @brief Print model data to ModelDataFile
         */
        void Print();

        /** @brief The sum of all the FMN+1 force terms
         */
        threevector CalculateAcceleration()const;

    public:
        ForceModel();
        ForceModel(std::string filename, float accuracy, 
            std::vector<ForceTerm*> forceterms, Matter *& sample, 
            PlasmaData & pdata);
        ForceModel(std::string filename, float accuracy, 
            std::vector<ForceTerm*> forceterms, Matter *& sample, 
            PlasmaData * pdata);
        ForceModel(std::string filename, float accuracy, 
            std::vector<ForceTerm*> forceterms, 
                Matter *& sample, PlasmaGrid_Data & pgrid);
        ForceModel(std::string filename, float accuracy, 
            std::vector<ForceTerm*> forceterms, Matter *& sample, 
            PlasmaGrid_Data & pgrid, PlasmaData &pdata);

        ~ForceModel(){};
        
        void CreateFile(std::string filename);
        double ProbeTimeStep()const;
        double UpdateTimeStep();

        /** @brief Implement Euler method to calculate motion
         *   
         *  @see Force(double timestep)
         */
        void Force();

        /** @brief Implement Euler method to calculate motion
         *   
         *  @see CalculateAcceleration()
         */
        void Force(double timestep);
};

#endif /* __FORCEMODEL_H_INCLUDED__ */
