#ifndef __PLASMADATA_H_INCLUDED__   // if PlasmaData.h hasn't been included yet...
#define __PLASMADATA_H_INCLUDED__

#include "threevector.h"
#include <vector>

// Structure containing all the information defining the plasma conditions
struct PlasmaData{

	// Parameters used to calculate particle fluxes
	double NeutralDensity;		// m^-3, Neutral Density
	double ElectronDensity;		// m^-3, Electron Density
	double IonDensity;			// m^-3, Electron Density
	double IonTemp;				// K, Ion Temperature
	double ElectronTemp;		// K, Electron Temperature
	double NeutralTemp;			// K, Neutral Temperature
	double AmbientTemp;			// K, Ambient Temperature
	double mi;					// mi is the mass of the gas
	threevector PlasmaVel;		// m s^-1, Plasma Velocity (Should eventually be normalised to sound speed cs)
	threevector Gravity;		// m s^-2, Acceleration due to gravity 
	threevector ElectricField;	// V m^-1, Electric field at dust location (Normalised later) 
	threevector MagneticField;	// T, Magnetic field at dust location (Normalised later)
};

struct PlasmaGrid_Data{

	// Plasma Parameters
	std::vector<std::vector<double>> Te;	// Te is the electron temperature stored in units of (J)
	std::vector<std::vector<double>> Ti;	// Ti is the ion temperature stored in units of (J)
	std::vector<std::vector<double>> na0;	// na0 is the ion density stored in units of (m^-3)
	std::vector<std::vector<double>> na1;	// na1 is the electron density stored in units of (m^-3)
	std::vector<std::vector<double>> po;	// po is the potential
	std::vector<std::vector<double>> ua0;	// ua0 is the ion drift velocity stored in units of (m s^-1)
	std::vector<std::vector<double>> ua1;	// ua1 is the electron drift velocity stored in units of (m s^-1)
	std::vector<std::vector<double>> bx;	// bx is the magnitude of the mangetic field in the x direction (T)
	std::vector<std::vector<double>> by;	// by is the magnitude of the mangetic field in the y direction (T)
	std::vector<std::vector<double>> bz;	// bz is the magnitude of the mangetic field in the z direction (T)
	std::vector<std::vector<double>> x;		// x is the position of the cells in x direction
	std::vector<std::vector<double>> z;		// z is the position of the cells in z direction
	std::vector<std::vector<int>> gridflag;	// gridflag is used by some of the plasma data to determine if cell is empty

	// Plasma Simulation Domain
	int gridx;			// gridx is the number of grid cells in x direction
	int gridz;			// gridz is the number of grid cells in z direction
	int gridtheta;		// gridtheta is the number of grid cells in theta direction
	double gridxmin;	// gridxmin is the minimum grid cell number in x direction
	double gridxmax;	// gridxmax is the maximum grid cell number in x direction
	double gridzmin;	// gridzmin is the minimum grid cell number in z direction
	double gridzmax;	// gridzmax is the maximum grid cell number in z direction
	double dlx;			// dlx is the grid spacing in the x direction
	double dlz;			// dlz is the grid spacing in the z direction

	// Basic Parameters defining plasma type.
	double mi;			// mi is the mass of the gas
	char gas;			// gas is a character specifying the plasma type (i.e Hydrogen)
	char device;		// device is a character specifying the machine type ('m', 'j', 'i', 'p' or 'd')

};
// PlasmaData():NeutralDensity(0.0), ElectronDensity(0.0), IonDensity(0.0),
// 		IonTemp(0.0), ElectronTemp(0.0), NeutralTemp(0.0),
// 		AmbientTemp(0.0), PlasmaVel(threevector()), Gravity(threevector()),
// 		ElectricField(threevector()), MagneticField(threevector()){}
// PlasmaGrid_Data():Te(std::vector< std::vector<double> >()),Ti(std::vector< std::vector<double> >()),
// 				na0(std::vector< std::vector<double> >()), na1(std::vector< std::vector<double> >()),
// 				po(std::vector< std::vector<double> >()), ua1(std::vector< std::vector<double> >()),
// 				ua0(std::vector< std::vector<double> >()), bx(std::vector< std::vector<double> >()),
// 				by(std::vector< std::vector<double> >()), bz(std::vector< std::vector<double> >()),
// 				x(std::vector< std::vector<double> >()), z(std::vector< std::vector<double> >()),
// 				gridflag(std::vector< std::vector<int> >()), gridx(0), gridz(0), gridtheta(0),
// 				gridxmin(0.0), gridxmax(0.0), gridzmin(0.0), gridzmax(0.0), dlx(0.0), dlz(0.0),
// 				mi(0.0), gas('n'), device('n'){}

#endif
