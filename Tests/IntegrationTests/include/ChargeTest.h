#include "ChargingModel.h"

double solveOML(double a, double guess, double iontemp, double etemp){
        C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
        if( a >= 1.0 ){
		static bool runOnce;
		WarnOnce(runOnce,"DeltaTot >= 1.0. DeltaTot being set equal to unity.");
		a = 1.0;
	}
	double b = iontemp/etemp;
	double C = Me/Mp;

	double x1 = guess - ( (( 1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b))/((a-1)*exp(-guess) - sqrt(C/b) ) );

	while(fabs(guess-x1)>1e-2){
		guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
	}
	return guess;
}

int ChargeTest(char Element){
	clock_t begin = clock();
	// ********************************************************** //
	// FIRST, define program default behaviour

	// Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
	char EmissivityModel = 'c'; 	// Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
	char ExpansionModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
													// and (z)ero expansion
	char HeatCapacityModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
	char BoilingModel = 'y'; 	// Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 

 	// Parameters describing the heating model
	double Size=5e-8; 		// m
	double Temp=280;		// K
	double Potential = 2.5;		// Normalised Potential
	Matter *Sample;			// Define the sample matter type

	PlasmaData *Pdata = new PlasmaData;
	Pdata->NeutralDensity = 3e19;		// m^-3, Neutral density
	Pdata->ElectronDensity = 8e17;	 	// m^-3, Electron density
	Pdata->IonTemp = 10*1.16e4;	 	// K, Ion Temperature
	Pdata->NeutralTemp = 10*1.16e4; 	// K, Neutral Temperature, convert from eV
	Pdata->ElectronTemp = 10*1.16e4;	// K, Electron Temperature, convert from eV
	threevector GravityForce(0, 0, -9.81);
	Pdata->Gravity = GravityForce;
	threevector Efield(1, -2, 3);
	Pdata->ElectricField = Efield;
	threevector Bfield(0.0, 0.0, 0.0);
	Pdata->MagneticField = Bfield;

	std::array<bool,CMN> ChargeModels  = {true,false,false};

	// Models and ConstModels are placed in an array in this order:
	std::array<char, 4> ConstModels =
		{ EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};

	if	(Element == 'W'){ 
		Sample = new Tungsten(Size,Temp,ConstModels);
	}else if (Element == 'B'){ 
		Sample = new Beryllium(Size,Temp,ConstModels);
	}else if (Element == 'F'){
		Sample = new Iron(Size,Temp,ConstModels);
	}else if (Element == 'G'){
		Sample = new Graphite(Size,Temp,ConstModels);
	}else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}
	Sample->set_potential(Potential);
	ChargingModel MyModel("Tests/IntegrationTests/Data/out_ConstantChargingTest.txt",1.0,ChargeModels,Sample,Pdata);

	MyModel.UpdateTimeStep();
	MyModel.Charge();

	double AnalyticPotential(0);
	std::cout << "\nSample->get_deltatot() = " << Sample->get_deltatot();
	if( Sample->get_deltatot() < 1.0 ){ // solveOML only defined for deltatot < 1.0
		AnalyticPotential = solveOML(Sample->get_deltatot(),Potential,Pdata->IonTemp,Pdata->ElectronTemp);
	}else{ // If the grain is in fact positive ...
		AnalyticPotential = solveOML(Sample->get_deltatot(),Potential,Pdata->IonTemp,Pdata->ElectronTemp);
		if( Potential < 0.0 ){
			AnalyticPotential = solveOML(0.0,Potential,Pdata->IonTemp,Pdata->ElectronTemp)
					-Kb*Sample->get_temperature()/(echarge*Pdata->ElectronTemp);
		}
	}

	double ReturnVal = 0;
	double ModelPotential = Sample->get_potential();

	if( AnalyticPotential == ModelPotential ) 					ReturnVal = 1;
	else if( fabs((AnalyticPotential-ModelPotential)/Potential < 0.01 ) ) 		ReturnVal = 2;
	else										ReturnVal = -1;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModelPot = " << ModelPotential << "arb : AnalyticPot = " << AnalyticPotential << "arb";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelPotential/AnalyticPotential) <<"%\n*****\n\n";

	delete Sample;

	return ReturnVal;
}
