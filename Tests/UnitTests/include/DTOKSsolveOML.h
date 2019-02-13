#ifndef __DTOKSSOLVEOML_H_INCLUDED__
#define __DTOKSSOLVEOML_H_INCLUDED__

double DTOKSsolveOML(double a, double Ti, double Te, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
	double b(0);
	try{
		b = Ti/Te;
	}catch(std::overflow_error e){
		std::cout << e.what();
		b = 0;
	}
	double C = Me/Mp;

	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/
		((a-1)*exp(-guess) - sqrt(C/b) ) );

	// With 1e-2 as the limit, we see some discontinuous behaviour around 
	// deltatot = 1.0
//	while(fabs(guess-x1)>1e-2){ 
	while(fabs(guess-x1)>1e-3){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /
			( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}
#endif
