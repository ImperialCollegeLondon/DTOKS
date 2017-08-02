#include <boost/math/tools/roots.hpp>

double solvePosSchottkyOML(double Ti, double Te, double Td, double ni, double ne, double Dsec, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");

	double WorkFunction = 3.4*echarge;	
	double Ve = sqrt(Kb*Te/Me);
	double Vi = sqrt(Kb*Ti/Mp);

	double Esee = 3*echarge;
	// For a negative dust grain...
	double K = Kb*Te/Esee;
	double A = ne*Ve*Dsec;
	double B = (Te*Dsec)/Ti;
	double C = -ne*Ve;
	double D = -Te/Ti;
	double H = ni*Vi;
	double I = 4*Richardson*pow(Td,2)/echarge;
	double f = WorkFunction/(Kb*Td);


	double fx 	= C+A*exp(K*guess)+B*guess*exp(K*guess)+D*guess+H*exp(-guess)+I*exp(guess-f);
	double fxprime 	= K*A*exp(K*guess)+B*exp(K*guess)+K*B*guess*exp(K*guess)+D-H*exp(-guess)+I*exp(guess-f);

	double x1 = guess - fx / fxprime;

//	std::cout << "\nx1 = " << x1 << "\tfx = " << fx << "\tfxprime = " << fxprime << "\tguess = " << guess;
//	std::cout << "\nT1 = " << C << "\tT2 = " << A*exp(K*guess) << "\tT3 = " << B*guess*exp(K*guess) << "\tT4 = " << D*guess << "\tT5 = " << H*exp(-guess) << "\tT6 = " << I*exp(guess-f); 
//	std::cout << "\nK = " << K << "\tguess = " << guess << "\tK*guess = " << K*guess; std::cin.get();

/*
	while(fabs(guess-x1)>1e-4){// >1e-2){

		guess = x1;
		fx 		= C+A*exp(K*guess)+B*guess*exp(K*guess)+D*guess+H*exp(-guess)+I*exp(guess-f);
		fxprime 	= K*A*exp(K*guess)+B*exp(K*guess)+K*B*guess*exp(K*guess)+D-H*exp(-guess)+I*exp(guess-f);

		x1 = guess - fx/fxprime;

//		std::cout << "\nx1 = " << x1 << "\tfx = " << fx << "\tfxprime = " << fxprime << "\tguess = " << guess << "\tf = " << f << "\texp(guess-f) = " << exp(guess-f);
//		std::cout << "\n\na*exp(-guess) = " << a*exp(-guess) << "\tb*guess = " << b*guess << "\td*exp(guess-f) = " << d*exp(guess-f) << "\n\n"; std::cin.get();
//		std::cout << "\nItot = " << a*exp(-guess) + b*guess + c + d*exp(guess-f); std::cin.get();

	}

*/
	return guess;
}

