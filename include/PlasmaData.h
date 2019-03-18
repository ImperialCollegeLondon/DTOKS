/** @file PlasmaData.h
 *  @brief Data structures defining Plasma paramaters and variables
 *  
 *  Structures containing all the information relevant to dust particle
 *  simulation. Three structures, one defining current plasma data, one for
 *  defining the plasma parameters across an irregular rectangular grid and one
 *  for defining an additional boundary for terminating the simulation
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#ifndef __PLASMADATA_H_INCLUDED__
#define __PLASMADATA_H_INCLUDED__

#include <vector>

#include "threevector.h"

/** @brief Structure containing all the information defining the plasma 
 * conditions in the immediate surroundings of the dust grain 
 */
struct PlasmaData{

    // Parameters used to calculate particle fluxes
    double NeutralDensity;     //!< m^-3
    double ElectronDensity;    //!< m^-3
    double IonDensity;         //!< m^-3
    double IonTemp;            //!< K
    double ElectronTemp;       //!< K
    double NeutralTemp;        //!< K
    double AmbientTemp;        //!< K
    double mi;                 //!< kg, mi is the mass of the atoms
    double Z;                  //!< (1/e), Z is the mean ionisation of the gas
    double A;                  //!< (1/mu), A is the Atomic Number of the gas
    threevector PlasmaVel;     //!< m s^-1, Plasma Velocity
    threevector Gravity;       //!< m s^-2, Acceleration due to gravity
    threevector ElectricField; //!< V m^-1, Electric field at dust location
    threevector MagneticField; //!< T, Magnetic field at dust location
};

/** @brief Defines the data to completely parameterise a rectangular grid of 
 * plasma parameters. This is made up of the values for the plasma parameters,
 * and variables defining the dimensions of the grid over which they're 
 * defined 
 */
struct PlasmaGrid_Data{

    /* Plasma Parameters */
    std::vector<std::vector<double>> Te;    //!< K, Electron temperature
    std::vector<std::vector<double>> Ti;    //!< K, Ion temperature
    std::vector<std::vector<double>> Tn;    //!< K, Neutral temperature
    std::vector<std::vector<double>> Ta;    //!< K, Ambient temperature
    std::vector<std::vector<double>> na0;   //!< m^-3, Ion density
    std::vector<std::vector<double>> na1;   //!< m^-3, Electron density
    std::vector<std::vector<double>> na2;   //!< m^-3, Neutral density
    std::vector<std::vector<double>> po;    //!< V, Potential
    std::vector<std::vector<double>> ua0;   //!< m s^-1, Ion drift vel
    std::vector<std::vector<double>> ua1;   //!< m s^-1, Electron drift vel
    std::vector<std::vector<double>> bx;    //!< T, Mang field, x direction
    std::vector<std::vector<double>> by;    //!< T, Mang field, y direction
    std::vector<std::vector<double>> bz;    //!< T, Mang field, z direction
    std::vector<std::vector<double>> x;     //!< m, Position of x cells
    std::vector<std::vector<double>> z;     //!< m, Position of z cells
    std::vector<std::vector<double>> dm;    //!< kg, Lost mass in every cell
    std::vector<std::vector<int>> gridflag; //!< Determine if cell is empty

    /* Plasma Simulation Domain */
    int gridx;       //!< the number of grid cells in x direction
    int gridz;       //!< the number of grid cells in z direction
    int gridtheta;   //!< the number of grid cells in theta direction
    double gridxmin; //!< the minimum grid cell number in x direction
    double gridxmax; //!< the maximum grid cell number in x direction
    double gridzmin; //!< the minimum grid cell number in z direction
    double gridzmax; //!< the maximum grid cell number in z direction
    double dlx;      //!< the grid spacing in the x direction
    double dlz;      //!< the grid spacing in the z direction

    /* Basic Parameters defining plasma type. */
    double mi;   //!< mi is the mass of the gas (kg)
    char device; //!< Specify the machine type ('m', 'j', 'i', 'p' or 'd')
};

/** @brief Two dimensional positional information defining a boundary which
 * exists in cylindrical r, z space. Used by DTOKSU to terminate simulations. 
 */
struct Boundary_Data{
    std::vector< std::pair<double,double> > Grid_Pos; //!< Defines Boundary
};

#endif /* __PLASMADATA_H_INCLUDED__ */

