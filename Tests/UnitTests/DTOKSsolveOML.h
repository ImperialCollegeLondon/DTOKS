
double DTOKSsolveOML(double Ti, double Te, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
	double a(0.0);
	double b(0);
	try{
		b = Ti/Te;
	}catch(std::overflow_error e){
		std::cout << e.what();
		b = 0;
	}
	double C = Me/Mp;

	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/((a-1)*exp(-guess) - sqrt(C/b) ) );

	while(fabs(guess-x1)>1e-2){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}
