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

	double x1 = guess;
	do{
		guess = x1;

		double fx = a*exp(-guess)+b*guess+C+d*exp(guess-f);	
		double fxprime = (-a*exp(-guess)+b+d*exp(guess-f));
		x1 = guess - fx/fxprime;

	}while(fabs(guess-x1)>1e-4);// >1e-2){

	return guess;
}

