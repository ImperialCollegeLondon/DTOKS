//#define PAUSE
//#define CHARGING_DEBUG

#include "ChargingModel.h"

ChargingModel::ChargingModel():Model(){
	C_Debug("\n\nIn ChargingModel::ChargingModel():Model()\n\n");
	// Charging Models turned on of possibl3 3
	UseModel[0] = true;
	UseModel[1] = false;
	UseModel[2] = false;
	ChargeOfGrain = 0;
	CreateFile("Default_Charging_Filename.txt");
}

ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,NumModels> models,
				Matter *& sample, PlasmaData *&pdata) : Model(sample,pdata,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename,float accuracy,std::array<bool,1> models,Matter *& sample, PlasmaData const& pdata) : Model(sample,pdata,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
	ChargeOfGrain = 0;
}

ChargingModel::ChargingModel(std::string filename, float accuracy, std::array<bool,NumModels> models,
				Matter *& sample, PlasmaGrid &pgrid) : Model(sample,pgrid,accuracy){
	C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename,float accuracy,std::array<bool,1> models,Matter *& sample, PlasmaGrid const& pgrid) : Model(sample,pgrid,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
	ChargeOfGrain = 0;
}

void ChargingModel::CreateFile(std::string filename){
	C_Debug("\tIn ChargingModel::CreateFile(std::string filename)\n\n");
	FileName = filename;
	ModelDataFile.open(FileName);
	ModelDataFile << std::fixed << std::setprecision(16) << std::endl;
	ModelDataFile << "Time\tChargeOfGrain";
	ModelDataFile << "\tPositive\tPotential";

	ModelDataFile << "\n";
	ModelDataFile.close();
	ModelDataFile.clear();
	Print();
}

void ChargingModel::Print(){
	C_Debug("\tIn ChargingModel::Print()\n\n");
	ModelDataFile.open(FileName,std::ofstream::app);
	ModelDataFile << TotalTime << "\t" << ChargeOfGrain/echarge;
	if( Sample->is_positive() )  ModelDataFile << "\tPos";
	if( !Sample->is_positive() ) ModelDataFile << "\tNeg";
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
	double DSec = DeltaSec();
	double DTherm = DeltaTherm();

	double Potential;
	if( Pdata->ElectronTemp <= 0 ){ // Avoid dividing by zero
/*		if( UseModel [0] ){	// This model is supposed to provide a steady loss of charge
			ChargeOfGrain = ChargeOfGrain + 4*PI*pow(Sample->get_radius()*Sample->get_temperature(),2)*Richardson*
			exp(-Sample->get_workfunction()*echarge/(Kb*Sample->get_temperature()))*timestep;
		}else if( UseModel [1] ){
			ChargeOfGrain += 4*PI*pow(Sample->get_radius()*Sample->get_temperature(),2)*Richardson*
			exp(Sample->get_potential()-Sample->get_workfunction()*echarge/(Kb*Sample->get_temperature()))*TimeStep;
		} */
		Potential = 0.0;
	}else if( UseModel[0] ){
//		WARNING! THIS SCHEME FOR THE CHARGING MODEL CREATES DISCONTINUITIES WHEN FORMING A WELL!
		if( (DSec + DTherm) >= 1.0 ){ 	// If electron emission yield exceeds unity, we have potential well...
						// Take away factor of temperature ratios for depth of well
						// This is not explained further...
			Potential = solveOML( 0.0, Sample->get_potential()) 
					- (Kb*Sample->get_temperature())/(Pdata->ElectronTemp*echarge); 
		}else{ // If the grain is negative...
			// Calculate the potential with the normal current balance including electron emission
			Potential = solveOML( DSec + DTherm,Sample->get_potential());
			if( Potential < 0.0 ){
				// But! If it's now positive, our assumptions must be wrong!
				// So now we assume it's positive and calculate the potential with a well.
				Potential = solveOML(0.0,Sample->get_potential())-(Kb*Sample->get_temperature())
						/(Pdata->ElectronTemp*echarge);
			}
		}
		ChargeOfGrain = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*Pdata->ElectronTemp)/echarge;
	}else if( UseModel[1] ){
		// Assume the grain is negative and calculate potential
		Potential = solveNegSchottkyOML(Potential);
		ChargeOfGrain = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*Pdata->ElectronTemp)/echarge;
		if ( Potential < 0 ){ // If dust grain is actually positive, (negative normalised potential) recalculate it...
			Potential = solvePosSchottkyOML(); // Quasi-neutrality assumed, we take ni as ni=ne
			ChargeOfGrain = (4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*Pdata->ElectronTemp)/echarge;
		}
	}else if( UseModel[2] ){	// In this case, maintain a potential well for entire temperature range
//		Potential = solveOML( -0.5, Sample->get_potential()) 
//			- Sample->get_temperature()/Pdata->ElectronTemp;
		Potential = solveOML2();//-Sample->get_temperature()/Pdata->ElectronTemp;
//		Potential = solveOML( DSec + DTherm, Sample->get_potential());
//		if( (DSec+ DTherm) >= 1.0 )
			
		ChargeOfGrain = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*Pdata->ElectronTemp)/echarge;
	}
	// Have to calculate charge of grain here since it doesn't know about the Electron Temp and since potential is normalised.
	// This information has to be passed to the grain.
	Sample->update_charge(ChargeOfGrain,Potential,DTherm,DSec);
	TotalTime += timestep;

	C_Debug("\t"); Print();
}

void ChargingModel::Charge(){
	C_Debug("\tIn ChargingModel::Charge()\n\n");
	Charge(TimeStep);
}

double ChargingModel::solveOML2(){
	double TemperatureRatio = Pdata->ElectronTemp/Pdata->IonTemp;
	double Ionization = 1.0;
	double MassRatio = Mp/Me;
	return -(TemperatureRatio/Ionization-LambertW(sqrt(MassRatio*TemperatureRatio)*exp(TemperatureRatio/Ionization)/Ionization));
}
double ChargingModel::solveOML(double a, double guess){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
        if( a >= 1.0 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"DeltaTot >= 1.0. DeltaTot being set equal to unity.");
		a = 1.0;
	}
	double b(0);
	if( Pdata->ElectronTemp != 0 ) b = Pdata->IonTemp/Pdata->ElectronTemp;

	double C = Me/Mp;

	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/((a-1)*exp(-guess) - sqrt(C/b) ) );

	while(fabs(guess-x1)>1e-2){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}

// Solve the orbital motion limited potential for small dust grains accounting for electron emission.
// WARNING: This is only valid for positively charged dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solvePosSchottkyOML(){

	double Wf = echarge*Sample->get_workfunction();
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double Z = 1.0; 		// Ionization
	double vi = sqrt(Kb*Pdata->IonTemp/(2*PI*Mp));
	double ve = sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
	double Arg = TemperatureRatio*exp(TemperatureRatio)*(vi+Richardson*pow(Sample->get_temperature(),2)
			*exp(-Wf/(Kb*Sample->get_temperature()))/(Z*echarge*Pdata->IonDensity))/(ve*(1-DeltaSec()));
	assert( Arg != INFINITY && Arg != -INFINITY );

	double Coeff = TemperatureRatio;
	return +Coeff - LambertW(Arg); 	// NOTE! Should be "return -Coeff + LambertW(Arg);" but we've reversed sign as equation
					// as potential definition is reversed ( Positive potential for positive dust in this case)
}

// Solve the orbital motion limited potential for small dust grains accounting for electron emission.
// WARNING: This is only valid for negatively charged dust grains.
double ChargingModel::solveNegSchottkyOML(double guess){
        double Vi = sqrt(Kb*Pdata->IonTemp/(2*PI*Mp));
        double Ve = sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));

        double a = Pdata->ElectronDensity*Ve*(DeltaSec()-1);
        double b = Pdata->IonDensity*Vi*(Pdata->ElectronTemp/Pdata->IonTemp);
        double C = Pdata->IonDensity*Vi;
        double d = Richardson*pow(Sample->get_temperature(),2)/echarge;
        double f = (echarge*Sample->get_workfunction())/(Kb*Sample->get_temperature());

	double x1 = guess;
	do{
		guess = x1;

		double fx = a*exp(-guess)+b*guess+C+d*exp(guess-f);	
		double fxprime = (-a*exp(-guess)+b+d*exp(guess-f));
		x1 = guess - fx/fxprime;

	}while(fabs(guess-x1)>1e-4);// >1e-2){

        return guess;
}

double ChargingModel::DeltaTherm()const{
	C_Debug("\tIn ChargingModel::DeltaTherm()\n\n");

	double dtherm(0.0);

	if( ElectronFlux(Sample->get_temperature()) > 0.0 )
		dtherm = (Richardson*pow(Sample->get_temperature(),2)*exp(-(Sample->get_workfunction()*echarge)
                                /(Kb*Sample->get_temperature())))/(echarge*ElectronFlux(Sample->get_temperature()));

	assert(dtherm >= 0.0 && dtherm == dtherm && dtherm != INFINITY );
	return dtherm;
}

double ChargingModel::DeltaSec()const{
	C_Debug("\tIn ChargingModel::DeltaSec()\n\n");
	double ConvertKtoev(8.6173303e-5);
	return sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem());
}
