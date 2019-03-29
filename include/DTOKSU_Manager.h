/** @file DTOKSU_Manager.h
 *  @brief Class Wrapper for configuring and running DTOKSU 
 *  
 *  This contains the class DTOKSU_Manager used for configuring and running 
 *  instances of DTOKSU correctly and safely. This class makes significant use
 *  of the config4cpp configuration library and optionally also netcdfcpp.h 
 *  libraries for reading .netcdfcpp files.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#ifndef __DTOKSU_Manager_H_INCLUDED__
#define __DTOKSU_Manager_H_INCLUDED__


#include <random>                     //!< for std::normal_distribution<> etc.
#include <chrono>                     //!< for chrono::high_resolution_clock;
#include <config4cpp/Configuration.h> //!< For Configuration
#ifdef NETCDF_SWITCH
#include <netcdfcpp.h>                //!< for reading netcdf files
#endif
#include <exception>                  //!< for Exception handling
#include <stdlib.h>

#include "DTOKSU.h"

struct PlasmaFileReadFailure : public std::exception {
   const char * what () const throw () {
      return "Exception raised opening/reading/closing PlasmaData file\n";
   }
};

/** @class DTOKSU_Manager
 *  @brief Class wrapping DTOKSU class for configuring and running simulations
 *  
 *  This class is used to accumulate all the user inputs together to construct
 *  DTOKSU simulations. This requires input of plasma data and dust grain data
 *  and specifying the models which are to be incorporated.
 */
class DTOKSU_Manager{

    private:
        /** @name Private Member data
         *  @brief DTOKSU object is constructed with material and plasma data
         *
         *  A pointer to an instance of DTOKSU is constructed by passing the
         *  \p Sample, \p Pgrid and/or \p Pdata information to it. Using this,
         *  it is then possible to run multiple simulations with varying initial
         *  conditions or continue a simulation using some of the terminating
         *  conditions of the previous simulation such is required in breakup.
         */
        ///@{
        DTOKSU *Sim;
        Matter *Sample;
        PlasmaGrid_Data Pgrid;
        PlasmaData Pdata;
        //!< FORCE MODEL NUMBER, the number of charge models
        const static unsigned int FMN = 9;
        // HEATING MODEL NUMBER, the number of charge models
        const static unsigned int HMN = 18;
        // CHARGE MODEL NUMBER, the number of charge models
        const static unsigned int CMN = 12;

        /** @brief Defines the status of configuration
         *
         *  The value of the configuration status is used for error handling and 
         *  controlling program flow.
         *  @see config_message()
         */
        int Config_Status;
        ///@}

        
        /** @brief Used to handle the input given by user
         */
        template<typename T> int input_function(int &argc, char* argv[], int &i, 
            std::stringstream &ss0, T &Temp)const;
        /** @brief Remind user of the command line options
         */
        void show_usage(std::string name)const;
        /** @brief Displays string relevant to \p Config_Status value
         */
        void config_message()const;
        /** @brief check that the plasma data is within expected range
         *  @see Overflows
         *  @see Underflows
         *  @return false, plasma is badly configured else true
         */
        bool check_pdata_range()const;

        /** @name Configuration functions
         *  @brief functions for configuring plasma grid and boundaries
         *  
         *  These functions take directories and filenames as inputs to look for
         *  files containing information about the background plasma and 
         *  simulation boundaries.
         *  @param plasma_dirname directory containing plasma data file
         *  @param dirname directory containing generic boundary data file
         *  @param filename name of file containing generic boundary data
         *  @param BD the variable in which boundary data is stored
         *  @param wall_dirname directory containing wall boundary data file
         *  @return 1 if configured correctly, 0 if not
         */
        ///@{
        int configure_plasmagrid(std::string plasma_dirname);
        int configure_boundary(std::string dirname, std::string filename, 
            Boundary_Data& BD);
        int configure_coregrid(std::string wall_dirname);
        ///@}
        /** @brief function to read plasma data in from a formatted text file
         * 
         *  @param plasma_dirname directory containing plasma data file
         *  @return a value corresponding to the status of file reading
         */
        int read_data(std::string plasma_dirname);
        #ifdef NETCDF_SWITCH
        /** @brief function to read plasma data for MPSI in .netcdf file format
         * 
         *  @param plasma_dirname directory containing plasma data file
         *  @return a value corresponding to the status of file reading
         */
        int read_MPSIdata(std::string plasma_dirname);
        #endif

        /** @brief called by DTOKSU_Manager::Run(), operates DTOKSU with breakup
         */
        void Breakup();

    public:

        /** @name Constructors
         *  @brief functions to construct DTOKSU_Manager class
         */
        ///@{
        /** @brief Default constructor
         *
         *  Sets \p Config_Status = -1, unconfigured
         */
        DTOKSU_Manager();
        /** @brief Command line constructor
         *
         *  Sets \p Config_Status = -1, unconfigured, taking command line inputs
         */
        DTOKSU_Manager( int argc, char* argv[] );
        /** @brief Command line and config file constructor
         *
         *  Sets \p Config_Status = -1, unconfigured, taking command line inputs
         *  which are over-ridden by a configuration file.
         */
        DTOKSU_Manager( int argc, char* argv[], std::string = "Config_Files/DTOKSU_Config.cfg" );
        ///@}

        ~DTOKSU_Manager(){
        };
        
        int Configure(std::string = "Config_Files/DTOKSU_Config.cfg");
        int Configure(int argc, char* argv[], std::string = "Config_Files/DTOKSU_Config.cfg");

        /** @name Model Tests
         *  @brief Small tests which instantiate the models and do one step
         *
         *  1) First, create Model object of relevant type, 
         *  2) Check sample is inside simulation domain
         *  3) Retrieve timescale of the physical model 
         *  4) Perform a single timestep of same size as timescale.
         */
        ///@{
        int ChargeTest(double accuracy,std::vector<CurrentTerm*> CurrentTerms);
        int ForceTest(double accuracy,std::vector<ForceTerm*> ForceTerms);
        int HeatTest(double accuracy,std::vector<HeatTerm*> HeatTerms);
        ///@}

        /** @brief If correctly configured, run DTOKSU 
         *  @return the result of DTOKSU::Run() or 1 if not configured.
         */
        int Run();
};

#endif /* __DTOKSU_Manager_H_INCLUDED__ */
