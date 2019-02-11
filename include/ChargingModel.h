/** @file ChargingModel.h
 *  @brief Contains a class defining the physics model for charging dust
 *  
 *  Class defining the charging models which affect the potential of a 
 *  dust grain in a magnetically confined fusion tokamak relevant plasma
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug Charging Model number should be in a namespace
 *  @bug FileName needs to be moved into Model.
 *  @bug investigate if UseModel can be moved to Model class
 */

#ifndef __CHARGINGMODEL_H_INCLUDED__
#define __CHARGINGMODEL_H_INCLUDED__

#include "Model.h"
#include "solveMOMLEM.h"

// CHARGE MODEL NUMBER, the number of charge models
const unsigned int CMN = 11;

/** @class Model
 *  @brief Defines the charging models which affect the potential of dust
 *  
 *  
 */
class ChargingModel : public Model{

    private:
        /** @brief Array defining the charging models being used
         */
        std::array<bool,CMN> UseModel;
        /** @brief The filename to which the results are written to
         */
        std::string FileName;
        /** @brief Print model data to ModelDataFile
         */
        void Print();
        /** @brief Function used to compare the results of the charging models
         */
        void Test();
        
        /** @name Charging Models
         *  @brief different models for determining the dust grain charge
         *
         *  These functions all return the normalised potential of the dust 
         *  grain following a current balance to the grain surface or semi-
         *  empirical model. The number of functions is equal to CMN+1
         *  @return double the normalised surface potential of dust grain
         *  @param a the total electron emission yield
         *  @param guess an initial value guess for root finding algorithm
         */
        ///@{
        /** @brief calculate normalised potential following OML charging model
         */
        double solveOML(double a, double guess);
        /** @brief follow OML charging model for positive dust grains
         */
        double solvePosOML(double a, double guess);
        /** @brief calculate normalised potential following PHL charging model
         */
        double solvePHL(double guess);
        /** @brief calculate normalised potential following CW charging model
         */
        double solveCW(double guess);
        /** @brief calculate normalised potential following THS charging model
         */
        double solveTHS();
        /** @brief calculate normalised potential following MOML charging model
         */
        double solveMOML();
        /** @brief calculate normalised potential following SOML charging model
         */
        double solveSOML(double guess);
        /** @brief calculate normalised potential following SMOML charging model
         */
        double solveSMOML(double guess);
        /** @brief calculate normalised potential following MOMLWEM model
         */
        double solveMOMLWEM(double guess);
        /** @brief calculate normalised potential following MOMLEM model
         */
        double solveMOMLEM();
        ///@}
        /** @brief Thermionic electron emission yield
         *  @return the thermionic electron emission yield coefficient
         */
        double DeltaTherm()const;
        /** @brief Thermionic electron emission yield with schottky correction
         *  @param Potential the potential of the dust grain surface
         *  @return the thermionic electron emission yield coefficient
         */
        double ThermFluxSchottky(double Potential)const;
        /** @brief Secondary electron emission yield
         *  @return the secondary electron emission yield coefficient
         */
        double DeltaSec()const;


    public:
        ChargingModel();
        ChargingModel(std::string filename, float accuracy, 
            std::array<bool,CMN> models, Matter *& sample, PlasmaData & pdata);
        ChargingModel(std::string filename, float accuracy, 
            std::array<bool,CMN> models, Matter *& sample, PlasmaData * pdata);
        ChargingModel(std::string filename, float accuracy, 
            std::array<bool,CMN> models, Matter *& sample, 
            PlasmaGrid_Data & pgrid);
        ChargingModel(std::string filename, float accuracy, 
            std::array<bool,CMN> models, Matter *& sample, 
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
