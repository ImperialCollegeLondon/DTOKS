/** @file DTOKSU.h
 *  @brief Contains a class which controls three physics models with matter
 *  
 *  This file contains the DTOKSU class which defines the DTOKSU program which
 *  is used to simulate solid and liquid dust grains in tokamak plasmas.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __DTOKSU_H_INCLUDED__
#define __DTOKSU_H_INCLUDED__

//#define DTOKSU_DEBUG
//#define DTOKSU_DEEP_DEBUG

#include <algorithm>

#include "HeatingModel.h"
#include "ForceModel.h"
#include "ChargingModel.h"

/** @brief default boundary data is an empty vector of pairs, i.e no data
 */
static struct Boundary_Data BoundaryDefaults = {
    std::vector<std::pair<double,double>> ()
};

/** @class DTOKSU
 *  @brief Class bringing together dust grain data and physical models
 *  
 *  This class simulates the motion, charging and heating of a dust grain in a 
 *  tokamak plasma. The charging model is assumed to act on the fastest time 
 *  scale and is updated every global timestep. The heating and force models
 *  compete and are solved to a relevant accuracy. Optional data specifying the
 *  tokamak wall and core boundaries can also be provided and methods for 
 *  performing reflection on these boundaries are provided.
 */
class DTOKSU{

    private:
        /** @name Private Member data
         *  @brief Objects for physical models, plasma and matter data
         *
         *  The physics models are encapsulated in the \p HM, \p FM and \p CM
         *  which represent the heating, force and charging models. These act
         *  upon the \p Sample which contains data structures for the dust grain
         *  sample. DTOKSU maintains a pointer to this here for easier access to
         *  this information. The \p WallBound and \p CoreBound are two vectors
         *  of pairs which are a series of points that define boundaries.
         *  \p TotalTime is used to record the total time taken to perform a 
         *  simulation and \p MyFile is a output file
         */
        ///@{
        double TotalTime;
        Matter *Sample;
        HeatingModel HM;
        ForceModel FM;
        ChargingModel CM;
        Boundary_Data WallBound, CoreBound;
        std::ofstream MyFile;
        ///@}

        /** @name Printing functions
         *  @brief Functions used to print data to a master simulation file
         *  @param filename is the name of the file to write to
         */
        ///@{
        void print();
        void create_file(std::string filename);
        ///@}

        /** @name Boundary functions
         *  @brief Functions used to check particle interaction with boundary
         */
        ///@{
        /** @brief perform specular reflection on particle incident on boundary
         */
        void SpecularReflection();
        /** @brief determine if the particle is inside or outside the boundary
         *  @param InOrOut specifies if we are checking whether it is in or out
         *  @return true if inside \p WallBound and false if outside.
         *  viceversa for \p CoreBound
         */
        bool Boundary_Check(bool InOrOut);
        ///@}

    public:
        //!< MODEL NUMBER, the number of physical models in DTOKS
        static const unsigned int MN = 3;   

        /** @name Constructors
         *  @brief functions to construct DTOKSU class
         *
         *  In all cases, we specify \p MN number of accuracies to solve the 
         *  physics models to, as well as the \p sample, \p heatmodels, 
         *  \p forcemodels and \p chargemodels to specify the simulation.
         */
        ///@{
        /** @brief pdata constructor.
         *
         *  @param alvls the accuracy levels for each of MN models
         *  @param sample pointer to reference of Matter object data
         *  @param pdata the plasma data used in the simulation
         *  @param HeatTerms pointers to Heating Terms used by HeatModel
         *  @param ForceTerms pointers to Force Terms used by ForceModel
         *  @param CurrentTerms pointers to Current Terms used by ChargingModel
         */
        DTOKSU( std::array<float,MN> alvls, Matter *& sample, PlasmaData &pdata,
            std::vector<HeatTerm*> HeatTerms, 
            std::vector<ForceTerm*> ForceTerms, 
            std::vector<CurrentTerm*> CurrentTerms);

        /** @brief pgrid constructor.
         *
         *  @param alvls the accuracy levels for each of MN models
         *  @param sample pointer to reference of Matter object data
         *  @param pgrid the plasma grid containing all plasma data
         *  @param HeatTerms pointers to Heating Terms used by HeatModel
         *  @param ForceTerms pointers to Force Terms used by ForceModel
         *  @param CurrentTerms pointers to Current Terms used by ChargingModel
         */
        DTOKSU( std::array<float,MN> alvls, Matter *& sample, 
            PlasmaGrid_Data &pgrid, std::vector<HeatTerm*> HeatTerms, 
            std::vector<ForceTerm*> ForceTerms, 
            std::vector<CurrentTerm*> CurrentTerms);

        /** @brief pdata and pgrid constructor.
         *
         *  @param alvls the accuracy levels for each of MN models
         *  @param sample pointer to reference of Matter object data
         *  @param pgrid the plasma grid containing all plasma data
         *  @param pdata the plasma data used in the simulation
         *  @param HeatTerms pointers to Heating Terms used by HeatModel
         *  @param ForceTerms pointers to Force Terms used by ForceModel
         *  @param CurrentTerms pointers to Current Terms used by ChargingModel
         */
        DTOKSU( std::array<float,MN> alvls, Matter *& sample, 
            PlasmaGrid_Data &pgrid, PlasmaData &pdata, 
            std::vector<HeatTerm*> HeatTerms, 
            std::vector<ForceTerm*> ForceTerms, 
            std::vector<CurrentTerm*> CurrentTerms);

        /** @brief boundary constructor.
         *
         *  @param alvls the accuracy levels for each of MN models
         *  @param sample pointer to reference of Matter object data
         *  @param pgrid the plasma grid containing all plasma data
         *  @param pdata the plasma data used in the simulation
         *  @param cbound the list of points defining the core boundary
         *  @param wbound the list of points defining the wall boundary
         *  @param HeatTerms pointers to Heating Terms used by HeatModel
         *  @param ForceTerms pointers to Force Terms used by ForceModel
         *  @param CurrentTerms pointers to Current Terms used by ChargingModel
         */
        DTOKSU( std::array<float,MN> alvls, Matter *& sample, 
            PlasmaGrid_Data &pgrid, PlasmaData &pdata, Boundary_Data &wbound, 
            Boundary_Data &cbound, std::vector<HeatTerm*> HeatTerms, 
            std::vector<ForceTerm*> ForceTerms, 
            std::vector<CurrentTerm*> CurrentTerms);
        ///@}

        ~DTOKSU(){
        };

        /** @brief Used to execute the program and solve the physics models
         *
         *  This function will continue to execute whilst the particle has not
         *  satisified any of the following conditions:
         *  1) Exited the simulation domain (return 1)
         *  2) Reached thermal equilibrium in a continuous plasma (return 2)
         *  3) Undergone a breakup process (return 3)
         *  4) Evaporated, boiled or fallen below lower mass limit (return 4)
         *  5) An error occures (return 5)
         *  @return An integer code refering to a particular exit condition
         */
        int Run();

        /** @brief Used to open all the model data files
         *
         *  @param filename files opened have the prefix filename
         *  @param i files opened have the suffix i
         */
        void OpenFiles(std::string filename, unsigned int i);
        /** @brief Used to close all the model data files
         */
        void CloseFiles();
        /** @brief Reset the time as recorded by each model if necessary
         */
        void ResetModelTime(double HMTime, double FMTime, double CMTime);
    
        /** @brief Prints to a file the data accumulated about impurity
         *  deposition in plasma grid
         */
        void ImpurityPrint();
        
        /** @name Public getter methods
         *  @brief functions required to get member data
         */
        ///@{
        double      get_HMTime()const   {   return HM.get_totaltime(); }
        double      get_FMTime()const   {   return FM.get_totaltime(); }
        double      get_CMTime()const   {   return CM.get_totaltime(); }
        threevector get_bfielddir()const{   return (FM.get_bfield());  }
        ///@}
};

#endif
