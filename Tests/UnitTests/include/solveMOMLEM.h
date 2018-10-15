#include "Constants.h"
#include "Functions.h"
#include <math.h>
#include <assert.h>	// for assert()

// Equation for PI_1 as defined in equation (47)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double Pi_1(double x){
	return exp(x*x)-exp(x*x)*erf(x)+(2.0*x)/sqrt(PI);
}

// Equation for PI_2 as defined in equation (50)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double Pi_2(double x){
	return exp(x*x)+exp(x*x)*erf(x)-(2.0*x)/sqrt(PI);
}

// Equation for PI_2 as defined in equation (22)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double NWSP_f(double x, double Tau, double Chi, double Delta, double Phic){
	double DensityRatio = 1.0;
	double T1(0.0), T2(0.0), T3(0.0);
	if( x >= 0.0 ){
		T1 = DensityRatio*exp(-x/Tau);
		T2 = 1.0*exp(x)*(1.0+erf(sqrt(x-Phic)));
		T3 = 1.0*Chi*exp((x-Phic)/Delta)*(1.0-erf(sqrt((x-Phic)/Delta)));
	}else{
		T1 = DensityRatio*exp(-x/Tau)*(1.0+erf(-sqrt(-x)/sqrt(Tau)));
		T2 = 1.0*exp(x)*(1.0+erf(sqrt(x-Phic)));
		T3 = 1.0*Chi*exp((x-Phic)/Delta)*(1.0-erf(sqrt((x-Phic)/Delta)));
	}
//	std::cout << "\nerf(-sqrt(-x)/sqrt(Tau)) = \t" << erf(-sqrt(-x)/sqrt(Tau));
//	std::cout << "\nerf(sqrt(x-Phic))) = \t" << erf(sqrt(x-Phic));
//	std::cout << "\nerf(sqrt((x-Phic)/Delta)) = \t" << erf(sqrt((x-Phic)/Delta));
	
//	std::cout << "\n\n* NWSP *\n\tPhic:\t" << Phic;
//	std::cout << "\tT:\t" << Tau << "\tC:\t" << Chi << "\tD:\t" << Delta << "\tx:\t" << x;
//	std::cout << "\nPhic = " << Phic;
//	std::cout << "\nx-Phic = " << x-Phic;
//	std::cout << "\n(x-Phic)/Delta = " << (x-Phic)/Delta;
//	
//	std::cout << "\n\nT1 = \t" << T1;
//	std::cout << "\nT2 = \t" << T2;
//	std::cout << "\nT3 = \t" << T3;
//	std::cout << "\n\n* Result = " << DensityRatio*T1-T2-T3 << " *\n";
	
	return DensityRatio*T1-T2-T3;
}

// Equation for Density balance in case with a potential well as defined in equation (53), (29), (31) [region B] and (33)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double WDB_f(double phis, double phic, double phiw, double Tau, double Chi, double Delta){
	double DensityRatio = 1.0;
	assert( phis>phiw );
	assert( phis*phiw >= 0.0 );
	double T1 = exp(phis)*(1.0-erf(sqrt(phis-phiw)));
	double T2 = Chi*exp((phis-phic)/Delta)*(1.0+erf(sqrt((phis-phiw)/Delta)));
	double T3(0.0);
	
	if( phis >= 0.0 ){
		T3 = exp(-phis/Tau)*(1.0+erf(sqrt(phis/Tau)));
	}else{
		T3 = exp(-phis/Tau)*(1.0+erf(-sqrt(-phis/Tau)));
	}

	return T1+T2-DensityRatio*T3;
}

// Equation for flux balance in case with a potential well as defined in equation (54), (35), (37) [region B], (39) and (41)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double WFB_f(double phic, double phiw, double Tau, double Chi, double Delta, double Mu){
	double DensityRatio = 1.0;
	double T1 = exp(phiw);
	double T2 = Chi*sqrt(Delta)*exp((phiw-phic)/Delta);
	double T3(0.0);
	
	if( phic >= 0.0 ){
		T3 = sqrt(Tau*Mu)*exp(-phic/Tau);
	}else{
		T3 = sqrt(Tau*Mu);
	}

	return T1-T2-DensityRatio*T3;
}

// Equation for density integral balance in case with a potential well as defined in equation (55), (46), (49) and (52)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double WIB_f(double phis, double phic, double phiw, double Tau, double Chi, double Delta){
	double DensityRatio = 1.0;
	assert( phis < 0.0 );
	assert( phiw < 0.0 );
	assert( phis>phiw );
	assert( phis*phiw >= 0.0 );
	double T1 = Tau*(Pi_1(sqrt(-phis/Tau))-Pi_1(sqrt(-phiw/Tau)));
	double T2 = exp(phiw)*(1.0-Pi_2(sqrt(phis-phiw)));
	double T3 = Chi*Delta*exp(phiw)*(1.0-Pi_1(sqrt((phis-phiw)/Delta)));
	
	return DensityRatio*T1-T2-T3;
}

// Function defined in equation (55) using (57), (59) & (61)
// See the paper by Minas and Nikoleta,
// N. Rizopoulou, A. P. L. Robinson, M. Coppins, and M. Bacharis, Phys. Plasmas 21, (2014).
double CWP_f(double x, double Tau, double Chi, double Delta, double Mu, char CritValue){


	double DensityRatio = 1.0;
	double Phic = log(sqrt(Tau*Mu)+Chi*sqrt(Delta));

	double Phis(-0.0001), NWSP(0.0);
	assert(Phic <= Phis);
	do{
//		std::cout << "\n" << NWSP << "\t" << Phis << "\t" << Phic;
		NWSP = NWSP_f(Phis,Tau,Chi,Delta,Phic);
		Phis = Phis-0.0001;
	}while( fabs(NWSP) > 0.001 && Phis > Phic );
	if( fabs(NWSP) > 0.001 ){
//		std::cerr << "\nWARNING! ROOT NOT FOUND!";
		return 0.0;
	}
//	std::cout << "\n\n* CWP *\n\tPhic:\t" << Phic << "\tPhis:\t" << Phis << "\tx:\t" << x;

	if( CritValue == 'c' ){
		std::cout << "\tT:\t" << Tau << "\tC:\t" << x << "\tD:\t" << Delta;
		return DensityRatio*Tau*exp(-Phic)*(Pi_1(sqrt(-Phis/Tau))-Pi_1(sqrt(-Phic/Tau)))
			+(1.0-Pi_2(sqrt(Phis-Phic)))+x*Delta*(1.0-Pi_1(sqrt((Phis-Phic)/Delta)));
	}else if( CritValue == 'd' ){
//		std::cout << "\tT:\t" << Tau << "\tC:\t" << Chi << "\tD:\t" << x;
//		std::cout << "\nPhis = " << Phis;
//		std::cout << "\tPhic = " << Phic;
//		std::cout << "\n\nT1 = " << DensityRatio*Tau*exp(-Phic)*(Pi_1(sqrt(-Phis/Tau))-Pi_1(sqrt(-Phic/Tau)));
//		std::cout << "\nT2 = " << (1.0-Pi_2(sqrt(Phis-Phic)));
//		std::cout << "\nT3 = " << x*Chi*(1.0-Pi_1(sqrt((Phis-Phic)/x)));

		return DensityRatio*Tau*exp(-Phic)*(Pi_1(sqrt(-Phis/Tau))-Pi_1(sqrt(-Phic/Tau)))
                        +(1.0-Pi_2(sqrt(Phis-Phic)))+x*Chi*(1.0-Pi_1(sqrt((Phis-Phic)/x)));
	}else if( CritValue == 't' ){
		std::cout << "\tT:\t" << x << "\tC:\t" << Chi << "\tD:\t" << Delta;
		return DensityRatio*x*exp(-Phic)*(Pi_1(sqrt(-Phis/x))-Pi_1(sqrt(-Phic/x)))
                        +(1.0-Pi_2(sqrt(Phis-Phic)))+Chi*Delta*(1.0-Pi_1(sqrt((Phis-Phic)/Delta)));
	}else{
		std::cerr << "\nUnrecognised CritValue. Please enter 'c' for Chi, 'd' for delta or 't' for tau";
	}

}

double FindCriticalVal(const double FirstFixedValue, const double SecondFixedValue, double MassRatio, char CritValue ){
//	std::cout << "\n\n* IN double FindCriticalVal(const double FirstFixedValue, const double SecondFixedValue, double MassRatio, char CritValue ) *\n";

	double DensityRatio = 1.0;
	double Tau(0.01), Chi(0.01), Delta(0.01), x(0.01);	
	if( CritValue == 'c' ){
		Tau = FirstFixedValue;
		Delta = SecondFixedValue;
	}else if( CritValue == 'd' ){
		Tau = FirstFixedValue;
		Chi = SecondFixedValue;

		double CWP_fResult(0.0);	
		do{
			CWP_fResult = fabs(CWP_f(x,Tau,Chi,x,1.0/MassRatio,CritValue));
			if( CWP_fResult == 0.0 ){
//				std::cerr << "\nINVALID REGIME!";
			}
//			std::cout << "\n\n" << x << "\t" << CWP_fResult;
			x = x+0.01;
	
//		}while( CWP_fResult > 0.001 && x < 1.0 );
		}while( x < 0.2 );
//		std::cout << "\nCWP_f = " << CWP_f(x,Tau,Chi,x,1.0/MassRatio,CritValue);
//		std::cout << "\nDelta = " << x;
	}else if( CritValue == 't' ){
		Chi = FirstFixedValue;
		Delta = SecondFixedValue;
	}else{
		std::cerr << "\nUnrecognised CritValue. Please enter 'c' for Chi, 'd' for delta or 't' for tau";
	}
//	std::cout << "\nT:\t" << Tau << "\tC:\t" << Chi << "\tD:\t" << Delta << "\tCV = " << CritValue;

	return x;

}

double solveWellCase( double &phis, double & phic, double &phiw, double tau, double chi, double delta, double mu){
//	std::cout << "\n\n* IN double FindCriticalVal(const double FirstFixedValue, const double SecondFixedValue, double MassRatio, char CritValue ) *\n";

	phis = -1.85;
	phic = -2.1;
	double steps = 0.0001;
	phiw = phic-steps;
	double dx = steps;
	double dy = steps;
	double dz = steps;
	double s = fabs(WDB_f(phis,phic,phiw,tau,chi,delta) + WIB_f(phis,phic,phiw,tau,chi,delta)
                + WFB_f(phic,phiw,tau,chi,delta,mu));
	double Diffs[6];
	int StepNumber(0.0);

	// CURRENTLY:
	// phis can only get more positive starting from -1.85
	// phiw can only get more negative starting from phic

	while( s > 0.0001 ){

		phis = phis+dx;
		double s1 = fabs(WDB_f(phis,phic,phiw,tau,chi,delta)) + fabs(WIB_f(phis,phic,phiw,tau,chi,delta))
				+ fabs(WFB_f(phic,phiw,tau,chi,delta,mu));
		if( fabs(s1) < fabs(s) ){
			Diffs[0] = fabs(s-s1);
		}else{
			Diffs[0] = 0.0;
		}
		phis = phis-dx;
		phic = phic+dy;
		double s2 = fabs(WDB_f(phis,phic,phiw,tau,chi,delta) + WIB_f(phis,phic,phiw,tau,chi,delta)
				+ WFB_f(phic,phiw,tau,chi,delta,mu));
		if( fabs(s2) < fabs(s) ){
			Diffs[1] = fabs(s-s2);
		}else{
			Diffs[1] = 0.0;
		}
		phic = phic-dy;

		phiw = phiw+dz;
		double s3 = fabs(WDB_f(phis,phic,phiw,tau,chi,delta) + WIB_f(phis,phic,phiw,tau,chi,delta)
				+ WFB_f(phic,phiw,tau,chi,delta,mu));
                if( fabs(s3) < fabs(s) ){
                        Diffs[2] = fabs(s-s3);
                }else{
			Diffs[2] = 0.0;
		}	
		phiw = phiw-dz;

/*
		phis = phis-dx;
		double s4 = fabs(WDB_f(phis,phic,phiw,tau,chi,delta) + WIB_f(phis,phic,phiw,tau,chi,delta)
				+ WFB_f(phic,phiw,tau,chi,delta,mu));
                if( fabs(s4) < fabs(s) ){
                        Diffs[3] = fabs(s-s4);
                }else{
			Diffs[3] = 0.0;
		}
		phis = phis+dx;
*/
		double s4 = 0.0;
		Diffs[4] = 0.0;
		phic = phic-dy;
		double s5 = fabs(WDB_f(phis,phic,phiw,tau,chi,delta) + WIB_f(phis,phic,phiw,tau,chi,delta)
				+ WFB_f(phic,phiw,tau,chi,delta,mu));
                if( fabs(s5) < fabs(s) ){
                        Diffs[4] = fabs(s-s5);
                }else{
			Diffs[4] = 0.0;
		}

		phic = phic+dy;
                phiw = phiw-dz;
                double s6 = fabs(WDB_f(phis,phic,phiw,tau,chi,delta) + WIB_f(phis,phic,phiw,tau,chi,delta)
                                + WFB_f(phic,phiw,tau,chi,delta,mu));
                if( fabs(s6) < fabs(s) ){
                        Diffs[5] = fabs(s-s6);
                }else{
			Diffs[5] = 0.0;
		}
		phiw = phiw+dz;
		unsigned int index(10);
		for( int i = 0; i < 6; i ++){
			double temp(0.0);
			if( Diffs[i] > temp ){
				temp = Diffs[i];
				index = i;
			}
		}
		if( index == 0 ){
			phis = phis+dx;
			s = s1;
		}else if( index == 1 ){
			phic = phic+dy;
			s = s2;
		}else if( index == 2 ){
			phiw = phiw+dz;
			s = s3;
		}else if( index == 3 ){
			phis = phis-dx;
			s = s4;
		}else if( index == 4 ){
			phic = phic-dy;
			s = s5;
		}else if( index == 5 ){
			phiw = phiw-dz;
			s = s6;
		}else{
			std::cerr << "\nERROR! MINIMISATION FAILED!!";
		}
	/*	if( phiw < -2.2 ){
			if( phic <= phis ){
				phiw = phic-dz;
			}else{
				phiw = phis-dz;
			}
		}*/
		std::cout << "\n" << StepNumber << " " << index << " " << phis << " " << phic << " " << phiw;
		std::cout << " " << s << " " << s1 << " " << s2 << " " << s3 << " " << s4;
		std::cout << " " << s5 << " " << s6;
		StepNumber ++;
	}

	return s;
}

void exploreWellCase( double tau, double chi, double delta, double mu){
//	std::cout << "\n\n* IN double FindCriticalVal(const double FirstFixedValue, const double SecondFixedValue, double MassRatio, char CritValue ) *\n";

	double steps = 0.01;
	double dx_max = -1.7;
	double dy_max = -1.6;
	double dz_max = -1.8;
	int StepNumber(0.0);

	// CURRENTLY:
	// phis can only get more positive starting from -1.85
	// phiw can only get more negative starting from phic

	for( double phis(-2.0); phis < dx_max; phis += steps){
		for( double phic(-2.1); phic < dy_max; phic += steps){
			for( double phiw(-2.2); phiw < dz_max; phiw += steps){
				double s=0.0;
				double s1=0.0;
				double s2=0.0;
				double s3=0.0;
				if( phis > phiw ){
					s = fabs(WDB_f(phis,phic,phiw,tau,chi,delta) + WIB_f(phis,phic,phiw,tau,chi,delta)
				                + WFB_f(phic,phiw,tau,chi,delta,mu));
					s1 = WDB_f(phis,phic,phiw,tau,chi,delta);
					s3 = WIB_f(phis,phic,phiw,tau,chi,delta);
				}

				s2 = WFB_f(phic,phiw,tau,chi,delta,mu);

				std::cout << "\n" << StepNumber << " "
					<< phis << " " << phic << " " << phiw << " " 
					<< s << " " << s1 << " " << s2 << " " << s3;
				StepNumber ++;
			}
		}
	}
}



// Solve the Modified orbital motion limited potential for large emitting dust grains.
// See the paper by Minas and Nikoleta, equation (1) and (2)
// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
double solveMOMLEM(double Tau, double MassRatio, double Chi, double Delta){
	double Mu = 1.0/MassRatio;
	double DensityRatio = 1.0;
	double Phic = log(sqrt(Tau*Mu)+Chi*sqrt(Delta));

	double Phis(-0.0001), NWSP(0.0);
	assert(Phic <= Phis);
	do{
		NWSP = NWSP_f(Phis,Tau,Chi,Delta,Phic);
		Phis = Phis-0.0001;
	}while( fabs(NWSP) > 0.001 && Phis > Phic );
	if( fabs(NWSP) > 0.001 ){
		Phis = Phic;
		double Phiw = Phic;
		NWSP = solveWellCase(Phis,Phic,Phiw,Tau,Chi,Delta,Mu);
//		std::cout << "\n" << Delta << " " << Phis << " " << Phic << " " << Phiw;
//		std::cerr << "\nWARNING! ROOT NOT FOUND!";
		return 0.0;
	}


	return Phis;
}

