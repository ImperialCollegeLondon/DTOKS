/** @file Model.h
 *  @brief Contains a class defining the behaviour of a physics model for dust
 *  
 *  Abstract base class which defines the behaviour of a physical model acting
 *  on the Matter class. Contains two structures which define default plasma 
 *  conditions.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 *  @bug Default structures should exist within a namespace
 *  @bug Consider encapsulating plasma data information in a separate class
 *  @bug InterpolatePlasmadata() is broken and only works for square grid
 */

#ifndef __MODEL_H_INCLUDED__
#define __MODEL_H_INCLUDED__

#include <memory>
#include <iomanip> //!< std::ofstream::setprecision()

#include "PlasmaData.h"
#include "Iron.h"
#include "Tungsten.h"
#include "Graphite.h"
#include "Beryllium.h"
#include "Deuterium.h"
#include "Lithium.h"
#include "Molybdenum.h"

/** @brief PlasmaData structure defining default values for plasma
 *
 *  This structure is used to construct the model class in the case where no 
 *  plasma data is provided
 */
static struct PlasmaData PlasmaDataDefaults = {
    1e19,            //!< m^-3, Neutral Density
    1e19,            //!< m^-3, Electron Density
    1e19,            //!< m^-3, Ion Density
    116045.25,       //!< K, Ion Temperature
    116045.25,       //!< K, Electron Temperature
    0.025*116045.25, //!< K, Neutral Temperature
    300,             //!< K, Ambient Temperature
    1.66054e-27,     //!< kg, Mass of ions
    1.0,             //!< (1/e), Mean Ionization state of plasma
    1.0,             //!< (1/mu), Mean Atomic Mass of plasma
    threevector(),   //!< m s^-1, Plasma Velocity
    threevector(),   //!< m s^-2, Acceleration due to gravity
    threevector(),   //!< V m^-1, Electric field at dust location
    threevector(),   //!< T, Magnetic field at dust location
};

/** @brief PlasmaGrid structure defining default values for plasma
 *
 *  This structure is used to construct the model class in the case where no 
 *  plasma grid is provided
 */
static struct PlasmaGrid_Data PlasmaGrid_DataDefaults = {
    // Plasma Parameters
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<double> >(),
    std::vector< std::vector<int> > (),

    251,
    401,
    0,
    1.5,
    4.0,
    -2.0,
    2.0,

    1.66054e-27,
    'j',
};

/** @class Model
 *  @brief Class defining the behaviour of a physics model for dust
 *  
 *  Abstract base class which defines the behaviour of a physical model acting 
 *  on the Matter class.
 */
class Model{
    private:
        /** @brief Shared_ptr to structure containing all plasma grid data
         *
         *  This structure contains the plasma parameters for a regular 
         *  rectangular grid within which the dust grain exists
         */
        std::shared_ptr<PlasmaGrid_Data> PG_data;
        
        /** @name Grid coordinates of dust
         *  @brief information regarding dust coordinates
         *
         *  Dust coordinates within two cylindrical dimensions of plasma grid
         */
        ///@{
        int i; // r Position of dust
        int k; // z Position of dust
        ///@}

        /** @name Private Member Functions
         *  @brief Functions for dust grain position to plasma grid 
         *
         *  These functions are used to identify dust position within the plasma
         *  grid.
         */
        ///@{
        /** @brief Determine the dust position in the grid
         *
         *  Mutate the values of \p i and \p k to be the x and z grid 
         *  coordinates of the dust grain
         *  return result of checkingrid()
         *  @see checkingrid()
         *  @return true if dust is within grid and false if it's not
         *  @param i grid coordinate of dust in x direction
         *  @param k grid coordinate of dust in z direction
         *  @param xd Threevector position of the dust 
         */
        const bool locate(int &i, int &k, threevector xd)const;

        /** @brief Determine whether dust position is within grid
         *
         *  return true if dust is within the grid and false if it's not
         *  @return true if dust is within grid and false if it's not
         *  @param i grid coordinate of dust in x direction
         *  @param k grid coordinate of dust in z direction
         *  @param xd Threevector position of the dust 
         */
        const bool checkingrid(int i, int k)const;

        /** @brief update fields in PlasmaData from PlasmaGrid at dust grain
         *  @param i grid coordinate of dust in x direction
         *  @param k grid coordinate of dust in z direction
         */
        void update_fields(int i, int k);
        ///@}

    protected:
        /** @name Protected Member Data
         *  @brief Functions for dust grain position to plasma grid 
         *
         *  These functions are used to identify dust position within the plasma
         *  grid and interpolate the plasma data
         */
        ///@{
        /** @brief Pointer to the Matter class containing material data
         *
         *  Contains all information relevant to dust grain sample
         */
        Matter *Sample;                    
        /** @brief Data structure containing local plasma parameters
         */
        std::shared_ptr<PlasmaData> Pdata;
        /** @brief Determines the accuracy to which the model is calculated
         */
        const float Accuracy;
        /** @brief Determines if plasma varies spatially (true) or not (false)
         */
        const bool ContinuousPlasma;
        /** @brief The timescale of the model
         */
        double TimeStep;
        /** @brief The total physics time for which simulation has been running
         */
        double TotalTime;
        /** @brief Data file where data relevant to physics model is printed
         */
        std::ofstream ModelDataFile;
        /** @brief Data file where plasma data is printed
         */
        std::ofstream PlasmaDataFile;
        ///@}

        /** @name Pure virtual functions
         *  @brief Functions for dust grain position to plasma grid 
         *
         *  These functions are used to identify dust position within the plasma
         *  grid and interpolate the plasma data
         */
        ///@{
        /** @brief Print the results of the model to \p ModelDataFile
         */
        virtual void Print()=0;
        /** @brief Update \p TimeStep of model to an appropriate timescale
         *  @return the time scale of the process
         */
        virtual double UpdateTimeStep()=0;
        /** @brief Getter method for \p TimeStep of model
         *  @return the time scale of the process
         */
        virtual double ProbeTimeStep()const=0;
        ///@}

        /** @name Plasma Fluxes
         *  @brief Functions which calculate the plasma flux to the sphere
         *
         *  These functions calculate the approximate flux to a charged sphere
         *  under different regimes relevant to different charging models. These
         *  are used to provide consistency in calculations made by different 
         *  models
         */
        ///@{
        /** @brief Calculate the Ion Flux as specified originally in DTOKS
         *  @param Potential the normalised potential on the dust grain
         *  @return the ion flux as originally given by DTOKS
         */
        const double DTOKSIonFlux(double Potential)const;
        /** @brief Calculate the Ion Flux as specified by OML
         *  @param Potential the normalised potential on the dust grain
         *  @return the ion flux as originally given by DTOKS
         */
        const double OMLIonFlux(double Potential)const;
        /** @brief Calculate the Ion Flux as specified by SOML
         *  @param Potential the normalised potential on the dust grain
         *  @return the ion flux as originally given by SOML
         */
        const double SOMLIonFlux(double Potential)const;
        /** @brief Calculate the Ion Flux as specified by SMOML
         *  @param Potential the normalised potential on the dust grain
         *  @return the ion flux as originally given by SMOML
         */
        const double SMOMLIonFlux(double Potential)const;
        /** @brief Calculate the Ion Flux as specified by PHL
         *  @param Potential the normalised potential on the dust grain
         *  @return the ion flux as originally given by PHL
         */
        const double PHLElectronFlux(double Potential)const;
        /** @brief Calculate the electron flux as specified originally in DTOKS
         *  @param Potential the normalised potential on the dust grain
         *  @return the electron flux as originally given by DTOKS
         */
        const double DTOKSElectronFlux(double Potential)const;
        /** @brief Calculate the electron flux as specified by OML
         *  @param Potential the normalised potential on the dust grain
         *  @return the ion flux as originally given by OML
         */
        const double OMLElectronFlux(double Potential)const;
        /** @brief Calculate the thermal neutral flux 
         *  @return the  thermal neutral flux
         */
        const double NeutralFlux()const;
        ///@}

    public:
         /** @name Constructors
         *  @brief functions to construct Model class
         */
        ///@{
        /** @brief Default constructor.
         *
         *  Create a tungsten matter sample with default plasma parameters,
         *  Accuracy = 1.0 and continuous plasma (=true).
         */
        Model();

        /** @brief PlasmaData constructor.
         *
         *  Set matter pointer to refer to \p sample, with plasma background 
         *  given by fixed spatially continuous values of \p pdata and with
         *  \p accuracy.
         *  @param sample the matter class which this model is acting on
         *  @param pdata the spatially continuous plasma paramater data
         *  @param accuracy the accuracy to which the model is calculated
         */
        Model(Matter *& sample, PlasmaData &pdata, float accuracy);

        /** @brief PlasmaGrid constructor.
         *
         *  Set matter pointer to refer to \p sample, with plasma background 
         *  given by the plasma grid specified by the \p pgrid structure with
         *  \p accuracy.
         *  @param sample the matter class which this model is acting on
         *  @param pgrid the spatial grid over which plasma parameters are given
         *  @param accuracy the accuracy to which the model is calculated
         */
        Model(Matter *& sample, PlasmaGrid_Data &pgrid, float accuracy);

        /** @brief PlasmaData and PlasmaGrid constructor.
         *
         *  Set matter pointer to refer to \p sample, with plasma background 
         *  given by the plasma grid specified by the \p pgrid structure with
         *  \p accuracy and \p pdata used as the initial conditions
         *  @param sample the matter class which this model is acting on
         *  @param pgrid the spatial grid over which plasma parameters are given
         *  @param pdata the spatially continuous plasma paramater data
         *  @param accuracy the accuracy to which the model is calculated
         */
        Model(Matter *& sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata, 
            float accuracy );
        ///@}

        virtual ~Model(){};

        /** @name Public getter methods
         *  @brief functions required to get member data
         */
        ///@{
        const threevector  get_bfield ()const
        { 
            return Pdata->MagneticField.getunit();         
        }
        double get_totaltime          ()const{ return TotalTime;    }
        double get_timestep           ()const{ return TimeStep;     }   
        const double get_dlx          ()const{ return PG_data->dlx; }
        ///@}

        /** @brief Determine whether the particle has entered a new cell
         *  @return True if the particle has entered a new cell, else false
         */
        bool new_cell                 ()const;

        /** @brief Write local plasma parameters to \p PlasmaDataFile
         *  @param filename location to write to
         */
        void RecordPlasmadata(std::string filename); // Record the plasma Data   

        /** @brief set the plasma data to the data structure passed
         *  @param pdata set the \p PlasmaData to \p pdata
         */
        void set_plasmadata(PlasmaData &pdata);

        /** @brief set the plasma grid to the data structure passed
         *  @param pgrid set the \p PG_data to \p pgrid
         */
        void set_plasmagrid(PlasmaGrid_Data &pgrid);
        
        /** @brief if not a continuous plasma, update the plasma from PlasmaGrid
         *  @see locate()
         *  @see update_fields()
         *  @return return the result of locate()
         */
        const bool update_plasmadata();

        /** @brief Create new file at \p filename to write model data to
         */
        virtual void CreateFile(std::string filename)=0;
        /** @brief close the \p ModelDataFile after writing
         */
        void close_file();

        /** @brief add to \p TotalTime, used to get correct timing across models
         *  @param T the time to be added to \p TotalTime
         */
        void AddTime(double T){ TotalTime = TotalTime + T;  }

        // Functions for printing, these haven't been validated yet
        // Print the inside and the outside of the tokamak 
        //void vtkcircle(double r, std::ofstream &fout); 
        //void vtktoroid();
        //impurityprint(double totalmass);
        //void datadump();

};

#endif /* __MODEL_H_INCLUDED__ */

        /** @brief Calculate interpolated plasma parameters at dust position
         *
         *  Attempt to mutate the mean value plasma parameters at the dust 
         *  position
         *  @bug This doesn't work for a non-square grid, or generally
         *  @param i grid coordinate of dust in x direction
         *  @param k grid coordinate of dust in z direction
         */
        //const void interpolatepdata(const int i,const int k)const;
