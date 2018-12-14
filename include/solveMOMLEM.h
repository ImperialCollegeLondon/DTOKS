#ifndef __SOLVEMOMLEM_H_INCLUDED__   // if ChargingModel.h hasn't been included yet...
#define __SOLVEMOMLEM_H_INCLUDED__

// Equation for PI_1 as defined in equation (47)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double Pi_1(double x);

// Equation for PI_2 as defined in equation (50)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double Pi_2(double x);

// Equation for PI_2 as defined in equation (22)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
// NWSP_f -> No Well Sheath Potential Function, solve for sheath potential
// x 		: Integration variable, sheath potential
// Tau 		: Themperature Ratio Ti/Te
// Chi 		: Density Ratio of Emitted to sheath Nem/Nse
// Delta 	: Temperature ratio of emitted to electron Tem/Te
// Phic		: Potential of sheath wall
double NWSP_f(double x, double Tau, double Chi, double Delta, double Phic);

// Equation for Density balance in case with a potential well as defined in equation (53), (29), (31) [region B] and (33)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double WDB_f(double phis, double phic, double phiw, double Tau, double Chi, double Delta);

// Equation for flux balance in case with a potential well as defined in equation (54), (35), (37) [region B], (39) and (41)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double WFB_f(double phic, double phiw, double Tau, double Chi, double Delta, double Mu);


// Equation for density integral balance in case with a potential well as defined in equation (55), (46), (49) and (52)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double WIB_f(double phis, double phiw, double Tau, double Chi, double Delta);

// Function defined in equation (55) using (57), (59) & (61)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
// CWP_f -> Critical Well Potential Function,
// x 		: Integration variable
// Tau 		: Themperature Ratio Ti/Te
// Chi 		: Density Ratio of Emitted to sheath Nem/Nse
// Delta 	: Temperature ratio of emitted to electron Tem/Te
// Mu		: Mass Ratio Mi/Me
// CrtiValue	: Character defining which critical value is being found
double CWP_f(double x, double Tau, double Chi, double Delta, double Mu, char CritValue);

double FindCriticalVal(const double FirstFixedValue, const double SecondFixedValue, double MassRatio, char CritValue );

double solveWellCase( double &phis, double & phic, double &phiw, double tau, double chi, double delta, double mu);

// Solve the Modified orbital motion limited potential for large emitting dust grains.
// See the paper by Minas and Nikoleta, equation (1) and (2)
// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
double solveDeltaMOMLEM(double Tau, double MassRatio, double Chi, double Delta);

#endif
