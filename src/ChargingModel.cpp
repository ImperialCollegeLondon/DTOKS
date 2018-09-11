//#define PAUSE
//#define CHARGING_DEBUG

#include "ChargingModel.h"

ChargingModel::ChargingModel():Model(){
	C_Debug("\n\nIn ChargingModel::ChargingModel():Model()\n\n");
	// Charging Models turned on of possibl3 3
	UseModel[0] = true;
	UseModel[1] = false;
	UseModel[2] = false;
	UseModel[3] = false;
	CreateFile("Default_Charging_Filename.txt");
}

ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,
				Matter *& sample, PlasmaData &pdata) : Model(sample,pdata,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,Matter *& sample, PlasmaData *&pdata) : Model(sample,pdata,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,
				Matter *& sample, PlasmaData *pdata) : Model(sample,*pdata,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,Matter *& sample, PlasmaData *&pdata) : Model(sample,pdata,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,
				Matter *& sample, PlasmaGrid_Data &pgrid) : Model(sample,pgrid,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,Matter *& sample, PlasmaGrid_Data &pgrid) : Model(sample,pgrid,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,
				Matter *& sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata) 
				: Model(sample,pgrid,pdata,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,CMN> models,Matter *& sample, PlasmaGrid_Data &pgrid) : Model(sample,pgrid,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}

void ChargingModel::CreateFile(std::string filename){
	C_Debug("\tIn ChargingModel::CreateFile(std::string filename)\n\n");
	FileName = filename;
	ModelDataFile.open(FileName);
	ModelDataFile << std::fixed << std::setprecision(16) << std::endl;
	ModelDataFile << "Time\tChargeOfGrain";
	ModelDataFile << "\tPositive\tDeltaTot\tPotential";
	ModelDataFile << "\n";
	ModelDataFile.close();
	ModelDataFile.clear();
	Print();
}

void ChargingModel::Print(){
	C_Debug("\tIn ChargingModel::Print()\n\n");
	ModelDataFile.open(FileName,std::ofstream::app);
	ModelDataFile << TotalTime << "\t" << -(4.0*PI*epsilon0*Sample->get_radius()*Sample->get_potential()*Kb*Pdata->ElectronTemp)/(echarge*echarge);
	if( Sample->is_positive() )  ModelDataFile << "\tPos";
	if( !Sample->is_positive() ) ModelDataFile << "\tNeg";
	ModelDataFile << "\t" << Sample->get_deltatot();
	ModelDataFile << "\t" << Sample->get_potential();
	ModelDataFile << "\n";
	ModelDataFile.close();
	ModelDataFile.clear();
}

double ChargingModel::ProbeTimeStep()const{
	C_Debug( "\tIn ChargingModel::ProbeTimeStep()\n\n" );

	double timestep(1.0);

	// Tests have shown that
	// Time step based on the electron plasma frequency. 
	if( Pdata->ElectronDensity != 0 )
		timestep = Accuracy*sqrt((epsilon0*Me)/(2*PI*Pdata->ElectronDensity*pow(echarge,2)));

	// Time step based on the formulation by Krasheninnikov
	// Smirnov, R. D., Pigarov, A. Y., Rosenberg, M., Krasheninnikov, S. I., & Mendis, D. a. (2007). 
	// Modelling of dynamics and transport of carbon dust particles in tokamaks. 
	// Plasma Physics and Controlled Fusion, 49(4), 347–371.
//	if( Pdata->ElectronDensity != 0 && Pdata->IonTemp != 0 ){
		// Calcualte the time scale of the behaviour from Krashinnenikovs equation
//		double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));
//		double PlasmaFreq = sqrt((Pdata->ElectronDensity*pow(echarge,2))/(epsilon0*Me));
//		timestep = sqrt(2*PI) * ((DebyeLength)/Sample->get_radius()) 
//				* (1/(PlasmaFreq*(1+Pdata->ElectronTemp/Pdata->IonTemp+fabs(Sample->get_potential()))))*Accuracy;
//		C_Debug("\n\t\tDebyeLength = " << DebyeLength << "\n\t\tPlasmaFreq = " << PlasmaFreq 
//			<< "\n\t\ttimestep = " << timestep << "\n\n");
//	}else{	timestep = 1; } // In region of no plasma


	assert(timestep == timestep);
	assert(timestep > 0);

	return timestep;
}

double ChargingModel::UpdateTimeStep(){
	C_Debug( "\tIn ChargingModel::UpdateTimeStep()\n\n" );
	TimeStep = ProbeTimeStep();
	return TimeStep;
}

void ChargingModel::Charge(double timestep){
	C_Debug("\tIn ChargingModel::Charge(double timestep)\n\n");

	// Make sure timestep input time is valid. Shouldn't exceed the timescale of the process.
	assert(timestep > 0);// && timestep <= TimeStep );
	double DSec = 0.0;
	double DTherm = 0.0;

	double Potential;
	if( UseModel[0] ){ // Original DTOKS Charging scheme
//		WARNING! THIS SCHEME FOR THE CHARGING MODEL CREATES DISCONTINUITIES WHEN FORMING A WELL!
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		if( (DSec + DTherm) >= 1.0 ){ 	// If electron emission yield exceeds unity, we have potential well...
						// Take away factor of temperature ratios for depth of well
						// This is not explained further...
			Potential = solveOML( 0.0, Sample->get_potential()) 
					- (Kb*Sample->get_temperature())/(Pdata->ElectronTemp*Kb); 
		}else{ // If the grain is negative...
			// Calculate the potential with the normal current balance including electron emission
			Potential = solveOML( DSec + DTherm,Sample->get_potential());
			if( Potential < 0.0 ){
				// But! If it's now positive, our assumptions must be wrong!
				// So now we assume it's positive and calculate the potential with a well.
				Potential = solveOML(0.0,Sample->get_potential())-(Kb*Sample->get_temperature())
						/(Pdata->ElectronTemp*Kb);
			}
		}
	}else if( UseModel[1] ){ // 
		// Assume the grain is negative and calculate potential
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		if ( (DSec + DTherm) >= 1.0 ){ // If dust grain is actually positive, (negative normalised potential) recalculate it...
			Potential = solvePosSchottkyOML(); // Quasi-neutrality assumed, we take ni as ni=ne
		}else{
			Potential = solveNegSchottkyOML(Potential);
			if( Potential < 0.0 ){
				Potential = solvePosSchottkyOML();
			}
		}
	}else if( UseModel[2] ){ // In this case, maintain a potential well for entire temperature range
		Potential = solveOML( 0.0, Sample->get_potential()) 
			- Kb*Sample->get_temperature()/(Pdata->ElectronTemp*Kb);
//		Potential = solveOML_LambertW(DSec+DTherm);//-Sample->get_temperature()/Pdata->ElectronTemp;
//		Potential = solveOML( DSec + DTherm, Sample->get_potential());
//		if( (DSec+ DTherm) >= 1.0 )
	}else if( UseModel[3] ){	// In this case, use Hutchinson, Patterchini and Lapenta
		Potential = solvePHL(Potential);
	}else if( UseModel[4] ){	// MOML-EM Charging model
		Potential = solveBIBHAS(Potential);
	}else if( UseModel[5] ){	// MOML Charging model
		Potential = solveMOML();
	}else if( UseModel[6] ){	// SOML Charging model
		Potential = solveSOML();
	}else if( UseModel[7] ){	// MOML-EM Charging model
		Potential = solveSMOML();
	}else if( UseModel[8] ){	// MOML-EM Charging model
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		Potential = solveMOMLWEM(DSec + DTherm);
	}
	/*else if( UseModel[8] ){	// MOML-EM Charging model
		Potential = solveMOMLEM(DSec + DTherm);
	}*/
	// Have to calculate charge of grain here since it doesn't know about the Electron Temp and since potential is normalised.
	// This information has to be passed to the grain.
	double charge = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*Pdata->ElectronTemp)/echarge;
	Sample->update_charge(charge,Potential,DTherm,DSec);
	TotalTime += timestep;

	C_Debug("\t"); Print();
}

void ChargingModel::Charge(){
	C_Debug("\tIn ChargingModel::Charge()\n\n");
	Charge(TimeStep);
}

double ChargingModel::solveOML(double a, double guess){
	C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
	if( a >= 1.0 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"DeltaTot >= 1.0. DeltaTot being set equal to unity.");
		a = 1.0;
	}
	double b = Pdata->IonTemp/Pdata->ElectronTemp;

	double C = Me/Pdata->mi;

	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/((a-1)*exp(-guess) - sqrt(C/b) ) );
	
	while( fabs(guess-x1) > Accuracy ){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}

// Solve the Modified orbital motion limited potential for large dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solveMOML(){
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double MassRatio = Pdata->mi/Me;
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double HeatCapacityRatio = 1.0;
    return -1.0*(TemperatureRatio
    		-LambertW(sqrt(2*PI*TemperatureRatio*(1+HeatCapacityRatio*TemperatureRatio))
    		*exp(TemperatureRatio))+log(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio)/2);
}

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solveSOML(){
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double MassRatio = Pdata->mi/Me;
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double PlasmaFlowSpeed = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	double s1 = sqrt(PI)*(1+2*pow(PlasmaFlowSpeed,2))*erf(PlasmaFlowSpeed)/(4*PlasmaFlowSpeed)+exp(-pow(PlasmaFlowSpeed,2))/2;
	double s2 = sqrt(PI)*erf(PlasmaFlowSpeed)/(2*PlasmaFlowSpeed);

	return TemperatureRatio*s1/s2-LambertW(sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio*s1/s2)/s2);
}

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solveSMOML(){
	double HeatCapacityRatio = 1.0;
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double MassRatio = Pdata->mi/Me;
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double PlasmaFlowSpeed = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;

    double s1 = sqrt(PI)*(1+2*pow(PlasmaFlowSpeed,2))*erf(PlasmaFlowSpeed)/(4*PlasmaFlowSpeed)+exp(-pow(PlasmaFlowSpeed,2))/2;
    double s2 = sqrt(PI)*erf(PlasmaFlowSpeed)/(2*PlasmaFlowSpeed);

    return TemperatureRatio*s1/s2-LambertW(sqrt(2*PI*TemperatureRatio
    	*(1+HeatCapacityRatio*TemperatureRatio))*exp(TemperatureRatio*s1/s2)/s2)
    	+log(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio)/2;
}

// Solve the Modified orbital motion limited potential for large emitting dust grains.
// See the paper by Minas and Nikoleta, equation (1) and (2)
// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
double ChargingModel::solveMOMLWEM(double DeltaTot){
	double HeatCapacityRatio = 1.0;
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double MassRatio = Pdata->mi/Me;
	double Ionization = Pdata->Z; 		// Ionization
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double PlasmaFlowSpeed = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	double Delta_Phi_em = 0.5*log((2.0*PI/MassRatio)*(1+HeatCapacityRatio*TemperatureRatio)/pow(1.0-DeltaTot,2.0));
    // Uncomment following line to compare MOMLWEM results with figure (1) of paper:
    // N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
    //std::cout << DeltaTot << "\t" << Delta_Phi_em << "\n";
    double Arg = (1.0-DeltaTot)*sqrt(MassRatio*TemperatureRatio)
                    *exp(TemperatureRatio/Ionization+Delta_Phi_em/Ionization)/Ionization;

    return TemperatureRatio/Ionization+Delta_Phi_em/Ionization-LambertW(Arg);
}

// Solve the orbital motion limited potential for small dust grains accounting for electron emission.
// WARNING: This is only valid for positively charged dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solvePosSchottkyOML(){

	double Wf = echarge*Sample->get_workfunction();
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double Ionization = Pdata->Z; 		// Ionization
	double vi = sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi));
	double ve = sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
	double Arg = TemperatureRatio*exp(TemperatureRatio)*(vi+Richardson*pow(Sample->get_temperature(),2)
			*exp(-Wf/(Kb*Sample->get_temperature()))/(Ionization*echarge*Pdata->IonDensity))/(ve*(1-DeltaSec()));

	if( std::isinf(Arg) || Arg < 0.0 ){
		std::cout << "\nLambertW(Arg == " << Arg << ")!";
		throw LambertWFailure();
		return 0;
	}
	double Coeff = TemperatureRatio;
	return +Coeff - LambertW(Arg); 	// NOTE! Should be "return -Coeff + LambertW(Arg);" but we've reversed sign as equation
					// as potential definition is reversed ( Positive potential for positive dust in this case)
}

// Solve the orbital motion limited potential for small dust grains accounting for electron emission.
// WARNING: This is only valid for negatively charged dust grains.
double ChargingModel::solveNegSchottkyOML(double guess){
	double Vi = sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi));
	double Ve = sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));

	double a = Pdata->ElectronDensity*Ve*(DeltaSec()-1);
	double b = Pdata->IonDensity*Vi*(Pdata->ElectronTemp/Pdata->IonTemp);
	double C = Pdata->IonDensity*Vi;
	double d = Richardson*pow(Sample->get_temperature(),2)/echarge;
	double f = (echarge*Sample->get_workfunction())/(Kb*Sample->get_temperature());

	double fx = a*exp(-guess)+b*guess+C+d*exp(guess-f);	
	double fxprime = (-a*exp(-guess)+b+d*exp(guess-f));

	double x1 = guess - fx/fxprime;
	while(fabs(guess-x1)>Accuracy){
		guess = x1;

		fx = a*exp(-guess)+b*guess+C+d*exp(guess-f);	
		fxprime = (-a*exp(-guess)+b+d*exp(guess-f));
		x1 = guess - fx/fxprime;

	}// >1e-2){

	return guess;
}


double ChargingModel::solveOML_LambertW(double DeltaTot){
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;

	double Ionization = Pdata->Z;
	double MassRatio = Pdata->mi/Me;
	double LW(0);
	try{
		double Arg = (1.0-DeltaTot)*sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/Ionization)/Ionization;
		if( std::isinf(Arg) ){
			std::cout << "\nLambertW(Arg == " << Arg << ")! Using solveOML";
			return solveOML( -0.5, Sample->get_potential()) 
					- Sample->get_temperature()/Pdata->ElectronTemp;
		}
		LW = LambertW((1.0-DeltaTot)*sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/Ionization)/Ionization);
		
	}catch( LambertWFailure &e ){
		std::cout << e.what();
	}
	return -(TemperatureRatio/Ionization-LW);
}

// Solve the Potential for a sphere in a collisionless magnetoplasma following
// L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).
// Assuming that Lambda_s >> r_p, i.e dust grain is small.
double ChargingModel::solvePHL(double Phi){ 
	double TemperatureRatio = Pdata->ElectronTemp/Pdata->IonTemp;
	double Beta = Sample->get_radius()
			/(sqrt(PI*Pdata->ElectronTemp*Me)/(2.0*echarge*echarge*Pdata->MagneticField*Pdata->MagneticField));
	double MassRatio = Pdata->mi/Me;

	if( Beta/MassRatio > 0.01 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"Beta/MassRatio > 0.01 in solvePHL! Model may not be valid in this range! see Fig 11. of  L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).");
	}
	double AtomicNumber = 1.0;	
	double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));

	double z = Beta/(1+Beta);
	double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+1.150*z*z*z*z*z;

	double eta = -(Phi/Beta)*(1+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
	double w(1.0);
	if( Beta == 0.0 || eta == -1.0 ){
		w = 1.0;
	}else if( std::isnan(eta) ){
		std::cout << "\nWarning! w being set to 1.0 (Assuming high B field limit) but Phi/Beta is nan while Beta != 0.";
		w = 1.0;
	}else{
		w = eta/(1+eta);
	} 
	
	double A = 0.678*w+1.543*w*w-1.212*w*w*w;

	double Arg = (A+(1.0-A)*i_star)*sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/AtomicNumber)/AtomicNumber;

	if( std::isinf(Arg) ){
		std::cout << "\nLambertW(Arg == " << Arg << ")!";
		throw LambertWFailure();
		return 0.0;
	}
	double Potential = TemperatureRatio/AtomicNumber
		-LambertW(Arg);
	
	while( fabs(Phi - Potential) >= Accuracy ){
		Phi = solvePHL(Potential);
	
		eta = -(Phi/Beta)*(1+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
		double w(1.0);
		if( Beta == 0.0 || eta == -1.0 ){
			w = 1.0;
		}else if( std::isnan(eta) ){
			std::cout << "\nWarning! w being set to 1.0 (Assuming high B field limit) but Phi/Beta is nan while Beta != 0.";
			w = 1.0;
		}else{
			w = eta/(1+eta);
		} 

	        A = 0.678*w+1.543*w*w-1.212*w*w*w;
		
		Arg = (A+(1.0-A)*i_star)*sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/AtomicNumber)/AtomicNumber;
		if( std::isinf(Arg) ){
			std::cout << "\nLambertW(Arg == " << Arg << ")!";
			throw LambertWFailure();
			return 0.0;
		}

		Potential = TemperatureRatio/AtomicNumber
			-LambertW(Arg);
	}
	return Potential;
}


// BIBHAS Charging Test:
// This test is designed to find the floating potential for arbitary sized dust grain
// The calculation follows the work by R. DE Bibhas,
// see R. DE Bibhas, Astrophys. Space Sci. 30, (1974).
double ChargingModel::solveBIBHAS(double guess){ 
        double MassRatio = Pdata->mi/Me;
        double Ionization = Pdata->Z;
        double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));
        double Radius = Sample->get_radius(); // Normalised
        double s = Radius+DebyeLength;
        
        double x1 = guess;

//      double r = Radius*(1.0-0.17/sqrt(Ionization))*sqrt(1.0/(2.0*Ionization))*(1.0+sqrt(1.0+2.0*Ionization));
//      double Pot_m = (Ionization*echarge*echarge/r - echarge*echarge*Radius/(2.0*(r*r-Radius*Radius))
//                      +ehcarge*echarge*Radius/(2*r*r))/(echarge*Pdata->ElectronTemp);
/*
        double Coeff = Ionization*sqrt(Pdata->Ion/temp/(Pdata->ElectronTemp*MassRatio))*(s*s/(Radius*Radius));
        double a = Coeff;
        double b = Coeff*((s*s-1.0)/(s*s));
        do{
                guess = x1;

                double fx = a-b*exp(guess/(s*s-1.0*1.0))-exp(guess);
                double fxprime = -(b/(s*s-1.0*1.0))*exp(guess*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-1.0*1.0)))-exp(guess);
                x1 = guess - fx/fxprime;

        }while(fabs(guess-x1)>1e-4);// >1e-2){
*/
        if( Pdata->IonTemp != Pdata->ElectronTemp ){
                std::cerr << "\nWarning! Mopel is invalid for Pdata->IonTemp != Pdata->ElectronTemp";
        }
        do{
                guess = x1;
                double u_d = sqrt(Kb*Pdata->IonTemp/(2.0*PI*MassRatio*Me))*exp(-guess)*(s*s/(Radius*Radius))
                        *(1.0-((s*s-Radius*Radius)/(s*s))*exp(guess*Radius*Radius*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-Radius*Radius))));
                double u_dprime = -1.0*u_d-sqrt(Kb*Pdata->IonTemp/(2.0*PI*MassRatio*Me))
                        *exp(-guess)*exp(guess*Radius*Radius*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-Radius*Radius)));

//              std::cout << "\nu_dprime = " << u_dprime;
//              std::cout << "\nexp(-x) = " << exp(-guess);
//              std::cout << "\nexp(x*a^2/(s^2-a^2)) = " << exp(guess*Radius*Radius*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-Radius*Radius)));
//              std::cout << "\nU = " << sqrt(Kb*Pdata->IonTemp/(2.0*PI*MassRatio*Me));

                double p = sqrt(MassRatio*Me/(2.0*Kb*Pdata->IonTemp))*u_d;
                double p_prime = sqrt(MassRatio*Me/(2.0*Kb*Pdata->IonTemp))*u_dprime;
                double q = sqrt((s*s/(s*s-Radius*Radius))*MassRatio*Me/(2.0*Kb*Pdata->IonTemp))*u_d;
                double q_prime = sqrt((s*s/(s*s-Radius*Radius))*MassRatio*Me/(2.0*Kb*Pdata->IonTemp))*u_dprime;

//              std::cout << "\nsqrt(MassRatio*Me/(2.0*Kb*Pdata->IonTemp)) = " << sqrt(MassRatio*Me/(2.0*echarge*Pdata->IonTemp));
//              std::cout << "\nu_d = " << u_d;
//              std::cout << "\nu_d_lim = " << sqrt(Kb*Pdata->IonTemp/(2.0*PI*MassRatio*Me))*exp(-guess)*(1.0-guess);
//              std::cout << "\n\nPot = " << guess;
//              std::cout << "\nPot_lim = " << 1.11;
//              std::cout << "\n\nP = " << p;
//              std::cout << "\nQ = " << q;
//              std::cout << "\n\nP' = " << p_prime;
//              std::cout << "\nQ' = " << q_prime;

                double Coeff = Ionization*sqrt(Pdata->IonTemp/(Pdata->ElectronTemp*MassRatio))*(s*s/(Radius*Radius));
                double a = Coeff*(exp(-p*p)+sqrt(PI)*p*(1.0+erf(p)));
                double b = Coeff*(((s*s-Radius*Radius)/(s*s))*exp(Radius*Radius*guess*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-Radius*Radius)))
                                *(exp(-q*q)+sqrt(PI)*q*(1.0+erf(q))));

                double a_prime = Coeff*(sqrt(PI)*p_prime*(1.0+erf(p)));
                double b_prime = Coeff*(((s*s-Radius*Radius)/(s*s))*exp(Radius*Radius*guess*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-Radius*Radius)))
                                *sqrt(PI)*q_prime*(1.0+erf(q))+b*Radius*Radius*Pdata->IonTemp/(Coeff*Pdata->ElectronTemp*(s*s-Radius*Radius)));
//              double b_prime = Coeff*(((s*s-Radius*Radius)/(s*s))*exp(Radius*Radius*guess*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-Radius*Radius)))
//                              *sqrt(PI)*q_prime*(1.0+erf(q))
//                              +(Radius*Radius*Pdata->IonTemp/(Pdata->ElectronTemp*s*s))*exp(Radius*Radius*guess*Pdata->IonTemp/(Pdata->ElectronTemp*(s*s-Radius*Radius)))
//                              *(exp(-q*q)+sqrt(PI)*q*(1.0+erf(q))));

                double fx = a-b-exp(guess);
                double fxprime = a_prime-b_prime-exp(guess);
                x1 = guess - fx/fxprime;

        }while(fabs(guess-x1)>1e-4);// >1e-2){


        return guess;
}

double ChargingModel::DeltaTherm()const{
	C_Debug("\tIn ChargingModel::DeltaTherm()\n\n");

	double dtherm = (Richardson*pow(Sample->get_temperature(),2)*exp(-(Sample->get_workfunction()*echarge)
					/(Kb*Sample->get_temperature())))/(echarge*ElectronFlux(Sample->get_temperature()));

	assert(dtherm >= 0.0 && dtherm == dtherm && dtherm != INFINITY );
	return dtherm;
}

double ChargingModel::DeltaSec()const{
	C_Debug("\tIn ChargingModel::DeltaSec()\n\n");
	double ConvertKtoev(8.6173303e-5);
	return sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem());
}
