/*
double solveFlippedOML(double Ti, double Te, double Td, double ni, double ne, double Dsec, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
	double Esee = 3*echarge;
	double Ve = sqrt(Kb*Te/Me);
	double Vi = sqrt(Kb*Ti/Mp);
	double WorkFunction = 3.4*echarge;	

	double xold(0.0);
	double a = (Te*ne*Ve*Dsec)/Ti;
	double b = (Kb*Te)/Esee;
	double C = (Te*ne*Ve)/Ti;
	double d = ni*Vi;
	double e = Richardson*pow(Td,2);
	double f = Te/Td;
	double g = WorkFunction/(Kb*Td);

//	double fx = (Dsec*exp(guess*Kb*Te/Esee)*(1-guess*Kb*Te/Esee)-1)*Ve*ne*(1-guess*Te/Ti)+ni*vi*exp(guess)+(Richardson*pow(Td,2)-guess*Te*Td)*exp(-(WorkFunction+guess*kB*Te)/(Kb*Td));
	double fx = -a*exp(-guess*b)*(2+b*guess)-C*guess+d*exp(-guess)+e*(1-f*guess)*exp(-g+guess*f)-Ve*ne;
	double fxprime = -a*exp(b*guess)*(pow(b*guess,2)+4*b*guess+2)-C-d*exp(guess)+e*pow(f*guess,2)*exp(-g+f*guess);
	double x1 = guess - fx/fxprime;
	
//	std::cout << -a*exp(-guess*b)*(2+b*guess)-C*guess << "\t" << d*exp(-guess) << "\t" << e*(1-f*guess)*exp(-g+guess*f) << "\t" << Ve*ne << "\t" << exp(-g) << "\t" << exp(guess*f);
//	std::cout << "\nfx = " << fx << "\nfxprime = " << fxprime;
//	std::cout << "\nguess = " << guess << "\nx1 = " << x1; std::cin.get();
	while(fabs(guess-x1)>1e-2){// >1e-2){

		if( x1 != xold ){
			xold = x1;
		}else if( xold == guess ) break;

		guess = x1;
		x1 = guess - (-a*exp(-guess*b)*(2+b*guess)-C*guess+d*exp(-guess)+e*(1-f*guess)*exp(-g+guess*f)-Ve*ne)/(-a*exp(b*guess)*(pow(b*guess,2)+4*b*guess+2)-C-d*exp(guess)+e*pow(f*guess,2)*exp(-g+f*guess));
		std::cout << "\nguess = " << guess << "\nx1 = " << x1; std::cin.get();

	}
	return guess;
}
*/

double solveFlippedOML(double Ti, double Te, double Td, double ni, double ne, double Dsec, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");

	double WorkFunction = 3.4*echarge;	
	double Phi = guess*Kb*Td/echarge;
	double Ve = sqrt(Kb*Te/Me);
	double Vi = sqrt(Kb*Ti/Mp);

	// For a negative dust grain...
	double a = ne*Ve*(Dsec-1);	
	double b = ni*Vi*(guess*Te)/Ti;
	double C = ni*Vi;
	double d = Richardson*pow(Td,2);
	double f = WorkFunction/(Kb*Td);
	double g = Te/Td;	

	double fx 	= a*exp(-guess)+b*guess+C+d*exp((echarge*Phi)/(Kb*Td)-f);
	double fxprime 	= -a*exp(-guess)+b+d*g*exp((echarge*Phi)/(Kb*Td)-f);

	double x1 = guess - fx / fxprime;

//	std::cout << "\nx1 = " << x1 << "\tguess = " << guess << "\tfx = " << fx << "\tfxprime = " << fxprime << "\texp((echarge*Phi)/(Kb*Td)-f) = " << exp((echarge*Phi)/(Kb*Td)-f); std::cin.get();
	while(fabs(guess-x1)>1e-2){// >1e-2){

		guess = x1;

		fx = a*exp(-guess)+b*guess+C+d*exp((echarge*Phi)/(Kb*Td)-f);	
		fxprime = (-a*exp(-guess)+b+d*g*exp((echarge*Phi)/(Kb*Td)-f));
		x1 = guess - fx/fxprime;
		


		Phi = guess*Kb*Td/echarge;
//		std::cout << "\nx1 = " << x1 << "\tguess = " << guess << "\tfx = " << fx << "\tfxprime = " << fxprime << "\texp((echarge*Phi)/(Kb*Td)-f) = " << exp((echarge*Phi)/(Kb*Td)-f);
//		std::cout << "\n(e*Phi) = " << echarge*Phi << "\tKb*Td = " << Kb*Td << "\t(e*Phi)/(Kb*Td) = " << (echarge*Phi)/(Kb*Td) << "\tf = " << f << "\t(echarge*Phi)/(Kb*Td)-f = " << (echarge*Phi)/(Kb*Td)-f; std::cin.get();
	}
	return guess;
}
