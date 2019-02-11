/** @file Constants.cpp
 *  @brief Implementation physical constants used in DTOKSU
 *  
 *  File defining physics constants which are required frequently enough to 
 *  justify their global definition
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug Constants should be within a namespace
 */

double Kb         = 1.38064852e-23;     //!< (kg m^2 s^-2 K^1) || (J K^-1)
double R          = 8.3144598;          //!< https://en.wikipedia.org/
                                        //!< wiki/Gas_constant 
double echarge    = 1.60217662e-19;     //!< C 
double Me         = 9.10938356e-31;     //!< kg, mass of electron
double Mp         = 1.66054e-27;        //!< kg, mass of ion (And Neutrons)
double AvNo       = 6.0221409e+23;      //!< mol^-1
double PI         = 3.14159265359;
double AMU        = 1.660539040e-27;    //!< kg, Atomic Mass unit
double Richardson = 1.20173e6;          //!< A/(m K^2), Richardson constant
double c          = 299792458;          //!< m/s, Speed of light
double h          = 6.62607004e-34;     //!< m^2 kg / s, Planck's Constant
double epsilon0   = 8.854187817620e-12; //!< F/m, vacuum permittivity
double Sigma      = 5.670373e-8;        //!< W/(m^2 K^4), Boltsmann constant
