/** @file ForceModel.h
 *  @brief Contains a class defining the physics model for equation of motion
 *  
 *  Class defining the force models which affect the equation of motion of a 
 *  dust grain in a magnetically confined fusion tokamak relevant plasma
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug Force Model number should be in a namespace
 *  @bug FileName needs to be moved into Model.
 *  @bug investigate if UseModel can be moved to Model class
 */

#ifndef __FORCEMODEL_H_INCLUDED__
#define __FORCEMODEL_H_INCLUDED__

#include "Model.h"
#include "MathHeader.h"

//!< FORCE MODEL NUMBER, the number of charge models
const unsigned int FMN = 10;

/** @class ForceModel
 *  @brief Defines the force models which affect the equation of motion of dust
 *  
 *  
 */
class ForceModel : public Model {

    private:
        /** @brief Array defining the force models being used
         */
        std::array<bool,FMN> UseModel; 
        /** @brief The filename to which the results are written to
         */
        std::string FileName;
        /** @brief The temperature of the dust grain in the previous step
         */
        double OldTemp;

        /** @brief Print model data to ModelDataFile
         */
        void Print();

        /** @name Force Terms
         *  @brief functions for different terms in force equation
         *
         *  These functions all return the acceleration due to this force term
         *  in the equation of motion for a dust grain in a tokamak plasma.
         *  The number of functions is equal to FMN+1
         *  @return threevector acceleration of dust due to force term
         */
        ///@{
        /** @brief The acceleration due to gravity
         */
        threevector Gravity()const;
        /** @brief SOML ion drag model due to collisions of dust with ions
         */
        threevector SOMLIonDrag()const;
        /** @brief SMOML ion drag model due to collisions of dust with ions
         */
        threevector SMOMLIonDrag()const;
        /** @brief DTOKS ion drag model due to collisions of dust with ions
         */
        threevector DTOKSIonDrag()const;
        /** @brief DUSTT ion drag model due to collisions of dust with ions
         */
        threevector DUSTTIonDrag()const;
        /** @brief Hybrid ion drag model due to collisions of dust with ions
         */
        threevector HybridIonDrag()const;
        /** @brief Neutral draf due to collisions of dust with neutrals
         */
        threevector NeutralDrag()const;
        /** @brief Lorentz force due to electric and magnetic fields
         */
        threevector LorentzForce()const;
        /** @brief Centrifugal force due to cylindrical coordinate system
         */
        threevector Centrifugal()const;
        /** @brief Rocket force due to assymetric ablation
         */
        threevector RocketForce()const;
        /** @brief The sum of all the FMN+1 force terms
         */
        threevector CalculateAcceleration()const;
        ///@}

    public:
        ForceModel();
        ForceModel(std::string filename, float accuracy, 
            std::array<bool,FMN> models, Matter *& sample, PlasmaData & pdata);
        ForceModel(std::string filename, float accuracy, 
            std::array<bool,FMN> models, Matter *& sample, PlasmaData * pdata);
        ForceModel(std::string filename, float accuracy, 
            std::array<bool,FMN> models, 
                Matter *& sample, PlasmaGrid_Data & pgrid);
        ForceModel(std::string filename, float accuracy, 
            std::array<bool,FMN> models, Matter *& sample, 
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
