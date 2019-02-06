/** @file GrainStructs.h
 *  @brief Data structures defining material constants and variables
 *  
 *  Structure containing all the intrinsic and extrinsic information about a 
 *  material comprised of a particular element.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug No known bugs.
 */

#ifndef __GRAINSTRUCTS_H_INCLUDED__
#define __GRAINSTRUCTS_H_INCLUDED__

#include "threevector.h"

/** @brief Data that varies with time and is generally a result of the 
 * simulation */
struct GrainData{

    bool Liquid, Gas; /* Define if the material is a liquid or a gas */
    bool Breakup;     /* True if the dust has undergone liquid breakup */

    /* ********************* Do vary with Temperature  ********************* */
    /* Geometric Data */
    double UnheatedRadius;   /* m */
    double Mass;             /* kg */
    double Radius;           /* m */
    double SurfaceArea;      /* m^2 */
    double Volume;           /* m^3 */
    double Density;          /* kg / m^3 */

    // Thermal Data
    double SuperBoilingTemp; /* K, Super heated boiling temperature */
    double Temperature;      /* K */
    double VapourPressure;   /* N/m^2 or Pa, Pascals */
    double Emissivity;       /* Arb, deviation from black body spectrum */
    double LinearExpansion;  /* m / K */
    double HeatCapacity;     /* kJ/(kg·K)    (Constant Pressure!) */

    // Charging Data
    double DeltaSec;         /* Secondary Electron Emission Yield (Arb) */
    double DeltaTherm;       /* Thermionic Electron Emission Yield (Arb) */
    double Potential;        /* Normalised Potential (-(e*phi) / (kB * Te)) */
    bool Positive;           /* Defines the sign of the grain charge. */

    // Force/Motion Data
    threevector DustPosition;   /* m, Cylindrical Coordinates */
    threevector DustVelocity;   /* m/s */
    double RotationalFrequency; /* s^-1, in B Field direction */

    // Latent Heat
    double FusionEnergy;        /* kJ, Energy of formed solid bonds */
    double VapourEnergy;        /* kJ, Energy of formed liquid bonds */
};

/** @brief Input parameters which are specific to the type of matter being 
 * tested, these are fixed once this structure is initialised. */
struct ElementConsts{
    /* Variable identifying the element material
     *(W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron. */
    char Elem;

    double MeltingTemp;    /* K */
    double BoilingTemp;    /* K */
    double BondEnergy;     /* kJ/mol, Energy to break an atomic bond */
    double WorkFunction;   /* eV, Energy to remove an electron */
    double HeatTransAir;   /* W/m^2 K */
    double AtomicMass;     /* kg/mol */
    double LatentFusion;   /* kJ/kg, Energy to melt 1kg of solid */
    double LatentVapour;   /* kJ/kg, Energy to vapourize 1kg of liquid */
    double SurfaceTension; /* N/m */
    double RTDensity;      /* Kg /m^3, Denisty at room temperature */
    double ThermConduct;   /* kW /m K */ 
};

#endif /* __GRAINSTRUCTS_H_INCLUDED__ */
