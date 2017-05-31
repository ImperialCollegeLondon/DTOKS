#ifdef MODEL_DEBUG
#define Mo_Debug(x) std::cout << x
#else
#define Mo_Debug(x) 
#endif 

#ifdef CHARGING_DEBUG
#define C_Debug(x) std::cout << x
#else
#define C_Debug(x) 
#endif 

#ifdef HEATING_DEBUG
#define H_Debug(x) std::cout << x
#else
#define H_Debug(x) 
#endif 

#ifdef HEATING_DEEP_DEBUG
#define H1_Debug(x) std::cout << x
#else
#define H1_Debug(x) 
#endif 

#ifdef FORCE_DEBUG
#define F_Debug(x) std::cout << x
#else
#define F_Debug(x) 
#endif 

#ifdef MATTER_DEBUG
#define M_Debug(x) std::cout << x
#else
#define M_Debug(x) 
#endif 

#ifdef MATTER_DEEP_DEBUG
#define M2_Debug(x) std::cout << x
#else
#define M2_Debug(x) 
#endif 

#ifdef ELEMENT_DEEP_DEBUG
#define E1_Debug(x) std::cout << x
#else
#define E1_Debug(x) 
#endif 

#ifdef ELEMENT_DEBUG
#define E_Debug(x) std::cout << x
#else
#define E_Debug(x) 
#endif 

#ifdef DTOKSU_DEBUG
#define D_Debug(x) std::cout << x
#else
#define D_Debug(x) 
#endif 

#ifdef PAUSE
#define Pause(); std::cin.get();
#else
#define Pause();
#endif

#define Kb 1.38064852e-23  		// (kg m^2 s^-2 K^1) || (J K^-1)
#define R 8.3144598 			// https://en.wikipedia.org/wiki/Gas_constant 
#define echarge 1.60217662e-19 		// C 
#define Me 9.11e-31			// kg, mass of electron
#define Mi 1.66054e-27			// kg, mass of ion (And Neutrons)
#define AvNo 6.0221409e+23 		// mol^-1
#define PI 3.14159265359
#define AMU 1.660539040e-27	 	// kg, Atomic Mass unit
#define Richardson 1.20173e6		// A/(metres K^2), Richardson constant
#define c 299792458 			// m/s, Speed of light
#define h 6.62607004e-34 		// m^2 kg / s, Planck's Constant
#define epsilon0 8.854187817620e-12 	// F/m, vacuum permittivity
#define Sigma 5.670373e-8 		// Boltsmann constant: W/(m^2 K^4) 
