// Solve the orbital motion limited potential for small dust grains accounting for electron emission.
// WARNING: This is only valid for negatively charged dust grains.

double solveNegSchottkyOML(double Ti, double Te, double Td, double ni, double ne, double Dsec, double guess){

	double WorkFunction = 3.4*echarge;
	double Vi = sqrt(Kb*Ti/(2*PI*Mp));
	double Ve = sqrt(Kb*Te/(2*PI*Me));

	// For a negative dust grain...
//	Dsec = 0.0;
	double a = ne*Ve*(Dsec-1);
	double b = ni*Vi*(Te/Ti);
	double C = ni*Vi;
	double d = Richardson*pow(Td,2)/echarge;
	double f = WorkFunction/(Kb*Td);

	double fx 	= a*exp(-guess)+b*guess+C+d*exp(guess-f);
	double fxprime 	= -a*exp(-guess)+b+d*exp(guess-f);

	double x1 = guess - fx / fxprime;

//	std::cout << "\nx1 = " << x1 << "\tfx = " << fx << "\tfxprime = " << fxprime << "\tguess = " << guess << "\tf = " << f << "\texp(guess-f) = " << exp(guess-f);
//	std::cout << "\n\na*exp(-guess) = " << a*exp(-guess) << "\tb*guess = " << b*guess << "\tC = " << C << "\td*exp(guess-f) = " << d*exp(guess-f) << "\td = " << d << "\n\n"; 
//	std::cout << "\nItot = " << a*exp(-guess) + b*guess + c + d*exp(guess-f); std::cin.get();

	while(fabs(guess-x1)>1e-4){// >1e-2){

		guess = x1;

		fx = a*exp(-guess)+b*guess+C+d*exp(guess-f);	
		fxprime = (-a*exp(-guess)+b+d*exp(guess-f));
		x1 = guess - fx/fxprime;
		
//		std::cout << "\nx1 = " << x1 << "\tfx = " << fx << "\tfxprime = " << fxprime << "\tguess = " << guess << "\tf = " << f << "\texp(guess-f) = " << exp(guess-f);
//		std::cout << "\n\na*exp(-guess) = " << a*exp(-guess) << "\tb*guess = " << b*guess << "\td*exp(guess-f) = " << d*exp(guess-f) << "\n\n"; std::cin.get();
//		std::cout << "\nItot = " << a*exp(-guess) + b*guess + c + d*exp(guess-f); std::cin.get();
	}
	return guess;
}

