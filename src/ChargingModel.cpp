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
	ModelDataFile << std::scientific << std::setprecision(16) << std::endl;
	ModelDataFile << "Time\tCharge\tSign\tDeltatot\tPotential\n";
	ModelDataFile.close();
	ModelDataFile.clear();
	Print();
}

void ChargingModel::Print(){
	C_Debug("\tIn ChargingModel::Print()\n\n");
	ModelDataFile.open(FileName,std::ofstream::app);
	ModelDataFile << TotalTime << "\t" 
				<< -(4.0*PI*epsilon0*Sample->get_radius()*Sample->get_potential()*Kb*Pdata->ElectronTemp)
				/(echarge*echarge);
	if( Sample->is_positive() )  ModelDataFile << "\tPos";
	if( !Sample->is_positive() ) ModelDataFile << "\tNeg";
	ModelDataFile << "\t" << Sample->get_deltatot() << "\t" << Sample->get_potential() << "\n";

	//Test();

	ModelDataFile.close();
	ModelDataFile.clear();
}

void ChargingModel::Test(){
	double DSec = DeltaSec();
	double DTherm = DeltaTherm();

	std::cout << "\n"; std::cout << Sample->get_temperature(); std::cout << "\t";
	std::cout << solveOML(DSec+DTherm,Sample->get_potential()); std::cout << "\t";
	std::cout << solveMOML(); std::cout << "\t";
	std::cout << solveSOML(Sample->get_potential()); std::cout << "\t";
	std::cout << solveSMOML(Sample->get_potential()); std::cout << "\t";
	std::cout << solveCW(Sample->get_potential()); std::cout << "\t";
	std::cout << solvePHL(Sample->get_potential()); std::cout << "\t";
	
	double Potential(0.0);
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
	std::cout << Potential; std::cout << "\t";
	std::cout << solveOML(0.0,Sample->get_potential())-Sample->get_temperature()/Pdata->ElectronTemp;
	//std::cout << solveMOMLWEM(Sample->get_potential())); std::cout << "\n");
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
	if( UseModel[0] ){	// Dynamically choose charging model
		double DebyeRatio = Sample->get_radius()/sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));
		double Betae = Sample->get_radius()
			/(sqrt(PI*Kb*Pdata->ElectronTemp*Me/(2.0*echarge*echarge*Pdata->MagneticField.mag3()*Pdata->MagneticField.mag3())));
		DSec = DeltaSec();
		DTherm = ThermFluxSchottky(Potential)/OMLElectronFlux(Sample->get_potential());
		if( (DSec+DTherm) <= Accuracy ){// Emission IS NOT important
			if( Betae <= 0.01*Accuracy ){ // electrons are un-magnetised
				// Chris-Willis fit to OM Theory
				//std::cout << "\nSolving CW!";
				Potential = solveCW(Sample->get_potential());

			}else{ // electrons are magnetised
				// THS fit to POT and DiMPl results DOESN'T EXIST YET
				// So instead, we do PHL
				//std::cout << "\nSolving PHL!";
				Potential = solvePHL(Sample->get_potential());
				if( Sample->get_potential() >= 0.0 ){
					DTherm = Richardson*pow(Sample->get_temperature(),2)*
						exp(-echarge*Sample->get_workfunction()/(Kb*Sample->get_temperature()))
						/(echarge*PHLElectronFlux(Sample->get_potential()));
				}else{	// OML Electron flux, SOML Ion Flux
					DTherm = Richardson*pow(Sample->get_temperature(),2)*(1.0-Sample->get_potential()
							*(Pdata->ElectronTemp/Sample->get_temperature()))
							*exp((-echarge*Sample->get_workfunction()+Sample->get_potential()*Kb*Pdata->ElectronTemp)
							/(Kb*Sample->get_temperature()))/(echarge*PHLElectronFlux(Sample->get_potential()));
				}
			}
		}else{ // Emission IS important
			if( DebyeRatio <= 0.1*Accuracy ){ // Small dust grains wrt the debye length
				// OMLWEM like Nikoleta's theory MOML-EM, DOESN'T EXIST YET
				// So instead, we do SOMLWEM
				//std::cout << "\nSolving SOMLWEM!";
				Potential = solveSOML( Sample->get_potential());
			}else{ // large dust grains wrt the debye length
				// MOML-EM since emission is far more important than magnetic field effects
				// HASN'T BEEN IMPLEMENTED YET, solve SMOMLWEM
				//std::cout << "\nSolving SMOMLWEM!";
				Potential = solveSMOML(Sample->get_potential());
			}
		}
		
	}else if( UseModel[1] ){	// OML Charging model
		Potential = solveOML(0.0,Potential);
		if( Potential < 0 ){ // In this case we have a positive grain!
			Potential = solvePosOML( 0.0, Sample->get_potential());
		}
	}else if( UseModel[2] ){	// MOML Charging model
		Potential = solveMOML();
		if( Potential < 0 ){ // In this case we have a positive grain!
			Potential = solvePosOML( 0.0, Sample->get_potential());
		}
	}else if( UseModel[3] ){	// SOML Charging model

		Potential = solveSOML(Sample->get_potential());
		DSec = DeltaSec();
		DTherm = ThermFluxSchottky(Potential)/OMLElectronFlux(Sample->get_potential());
	}else if( UseModel[4] ){	// SMOML Charging model
		Potential = solveSMOML(Sample->get_potential());
		DSec = DeltaSec();
		DTherm = ThermFluxSchottky(Potential)/OMLElectronFlux(Sample->get_potential());
	}else if(UseModel[5] ){		// CW Charging Model, Chris willis fit to Sceptic results
		Potential = solveCW(Sample->get_potential());
		DSec = DeltaSec();
		DTherm = ThermFluxSchottky(Potential)/OMLElectronFlux(Sample->get_potential());
	}else if( UseModel[6] ){	// In this case, use Hutchinson, Patterchini and Lapenta
		Potential = solvePHL(Sample->get_potential());
		DSec = DeltaSec();
		DTherm = ThermFluxSchottky(Potential)/PHLElectronFlux(Sample->get_potential());

	}else if( UseModel[7] ){ // Original DTOKS Charging scheme
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
	}else if( UseModel[8] ){ // In this case, maintain a potential well for entire temperature range
		Potential = solveOML( 0.0, Sample->get_potential()) 
			- Kb*Sample->get_temperature()/(Pdata->ElectronTemp*Kb);
		if( Potential < 0 ){ // In this case we have a positive grain! Assume well disappears
			Potential = solvePosOML( 0.0, Sample->get_potential());
		}
	}else if( UseModel[9] ){	// MOMLWEM Charging model
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		
		if( (DSec + DTherm) >= 1.0 ){ // In this case we have a positive grain!
			Potential = solvePosOML(0.0, Sample->get_potential());
		}else{
			Potential = solveMOMLWEM(DSec + DTherm);
			if( Potential < 0.0 ){
				Potential = solvePosOML(0.0, Sample->get_potential());
			}
		}
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

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
// or
// R. D. Smirnov, A. Y. Pigarov, M. Rosenberg, S. I. Krasheninnikov, and D. a Mendis, 
// Plasma Phys. Control. Fusion 49, 347 (2007). Equation [1] & [2] 
double ChargingModel::solveSOML(double guess){
	C_Debug("\tIn ChargingModel::solveSOML(double guess)\n\n");

	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;

	
	double IonFluxDiff(0.0), ElectronFluxDiff(0.0), x1(0.0), dT(0.0), dTdiff(0.0);
	int i(0);
	double guessinit = guess;
	do{
		if( i > 0 ) // If first time in loop, don't update guess. Avoids numerical instabilities.
			guess = x1;
		if( i == 1000 ){ // Too many loops! Maybe we need to flip sign?

			std::cerr << "\nError in ChargingModel::solveSOML()! Newton Rhapson Method didn't converge to accuracy: " << Accuracy;
			return 0.0;
		}

		if( guessinit >= 0.0 ){ // OML Electron flux, SOML Ion flux
			if( uz == 0.0 ){
				IonFluxDiff = Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*Pdata->Z/Tau;
			}else{
				IonFluxDiff = (Pdata->IonDensity*IonThermalVelocity/(2.0*sqrt(2.0)*uz))*
					Pdata->Z*erf(uz)/Tau;
			}
			ElectronFluxDiff = -OMLElectronFlux(fabs(guess));
			dT = Richardson*pow(Sample->get_temperature(),2)*exp(-echarge*Sample->get_workfunction()/(Kb*Sample->get_temperature()))/echarge;
			dTdiff = 0.0;
		}else{ // OML Electron flux, SOML Ion flux
			double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
			double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
			if( uz == 0.0 ){
				IonFluxDiff = (Pdata->ElectronTemp/Pdata->IonTemp)*SOMLIonFlux(guess);
			}else{
				IonFluxDiff = Pdata->IonDensity*
					IonThermalVelocity*(1.0/(4.0*sqrt(2.0)*uz))*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm));
			}
			
			ElectronFluxDiff = -Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
			
			dTdiff = (Pdata->ElectronTemp/Sample->get_temperature())*(dT-
					Richardson*pow(Sample->get_temperature(),2)
					*exp((-echarge*Sample->get_workfunction()-guess*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge);
		}

		x1 = guess - ( (( 1.0-DeltaSec())*OMLElectronFlux(guess) - dT - SOMLIonFlux(guess))
			/((1.0-DeltaSec())*ElectronFluxDiff - dTdiff - IonFluxDiff ) );
		i ++;
	}while( fabs(guess-x1) > Accuracy );
	assert( dT != INFINITY );
	return guess;
}

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solveSMOML(double guess){
	C_Debug("\tIn ChargingModel::solveSMOML(double guess)\n\n");

	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	
	
	double IonFluxDiff(0.0), ElectronFluxDiff(0.0), x1(0.0), dT(0.0), dTdiff(0.0);
	int i(0);
	double guessinit = guess;
	do{
		if( i > 0 ) // If first time in loop, don't update guess. Avoids numerical instabilities.
			guess = x1;

		if( i == 2000 ){
			std::cerr << "\nError in ChargingModel::solveSMOML()! Newton Rhapson Method didn't converge to accuracy: " << Accuracy;
			return 0.0;
		}
		if( guessinit >= 0.0 ){  // OML Electron flux, SMOML Ion Flux
			if( uz == 0.0 ){
				IonFluxDiff = Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0))*(Pdata->Z/Tau);
			}else{
				IonFluxDiff = (Pdata->IonDensity*IonThermalVelocity/(2.0*sqrt(2.0)*uz))*
					Pdata->Z*erf(uz)/Tau;
			}
			ElectronFluxDiff = -OMLElectronFlux(guess);
			dT = Richardson*pow(Sample->get_temperature(),2)*exp(-echarge*Sample->get_workfunction()/(Kb*Sample->get_temperature()))/echarge;
			dTdiff = 0.0;
		}else{	// OML Electron flux, SOML Ion Flux
			double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
			double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
			if( uz == 0.0 ){
				IonFluxDiff = Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0))*(Pdata->Z/Tau);
			}else{
				IonFluxDiff = Pdata->IonDensity*
					IonThermalVelocity*(1.0/(4.0*sqrt(2.0)*uz))*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm));
			}
			
			ElectronFluxDiff = -Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
			
			dT = Richardson*pow(Sample->get_temperature(),2)*(1.0-guess
					*(Pdata->ElectronTemp/Sample->get_temperature()))
					*exp((-echarge*Sample->get_workfunction()-guess*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge;
			dTdiff = (Pdata->ElectronTemp/Sample->get_temperature())*(dT-
					Richardson*pow(Sample->get_temperature(),2)
					*exp((-echarge*Sample->get_workfunction()-guess*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge);
		}


		x1 = guess - ( ((1.0 - DeltaSec())*OMLElectronFlux(guess) - dT - SMOMLIonFlux(guess))
			/((1.0-DeltaSec())*ElectronFluxDiff - dTdiff - IonFluxDiff ) );
		i ++;

	}while( fabs(guess-x1) > Accuracy );
	assert( dT != -INFINITY );
	return guess;
}

// Solve the Potential for a sphere in a collisionless magnetoplasma following
// L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).
// Assuming that Lambda_s >> r_p, i.e dust grain is small.
double ChargingModel::solvePHL(double guess){
	C_Debug("\tIn ChargingModel::solvePHL(double guess)\n\n");

	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;

	double Beta = Sample->get_radius()
			/(sqrt(PI*Kb*Pdata->ElectronTemp*Me/(2.0*echarge*echarge*Pdata->MagneticField.mag3()*Pdata->MagneticField.mag3())));
	
	double AtomicNumber = Pdata->Z;	
	double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));

	double z = Beta/(1+Beta);
	double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+1.150*z*z*z*z*z;
	double eta = (guess/Beta)*(1.0+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
	double d_eta = (1.0/Beta)*(1.0+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
	
	double w = (eta/(1.0+eta));
	double d_w = d_eta/pow(1.0+eta,2.0);
	double d_A = 0.678*d_w+1.543*2.0*w*d_w-1.212*3.0*w*w*d_w;

	double IonFluxDiff(0.0), ElectronFluxDiff(0.0), x1(0.0), dT(0.0), dTdiff(0.0);
	int i(0);
	double guessinit = guess;
	do{
		if( i > 0 ) // If first time in loop, don't update guess. Avoids numerical instabilities.
			guess = x1;
		if( i == 1000 ){ // Too many loops! Maybe we need to flip sign?
			guessinit = -guessinit;
		}
		if( i == 2000 ){
			//std::cerr << "\nError in ChargingModel::solvePHL()! Newton Rhapson Method didn't converge to accuracy: " << Accuracy;
			return 0.0;
		}
		if( guess >= 0.0 ){  // PHL Electron flux, SOML Ion Flux
			eta = -(guess/Beta)*(1.0+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
			w = (eta/(1.0+eta));
			d_w = (d_eta/(1.0+eta)-eta*d_eta/pow(1.0+eta,2.0));
			d_A = 0.678*d_w+1.543*2.0*w*d_w-1.212*3.0*w*w*d_w;
			ElectronFluxDiff=Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(d_A-d_A*i_star)-PHLElectronFlux(guess);
			IonFluxDiff = (Pdata->IonDensity*IonThermalVelocity/(2.0*sqrt(2.0)*uz))*
				Pdata->Z*erf(uz)/Tau;
			dT = Richardson*pow(Sample->get_temperature(),2)*exp(-echarge*Sample->get_workfunction()/(Kb*Sample->get_temperature()))/echarge;
			dTdiff = 0.0;
			x1 = guess - ( ((1.0-DeltaSec())*PHLElectronFlux(guess) - dT - SOMLIonFlux(guess))/((1.0-DeltaSec())*ElectronFluxDiff - dTdiff - IonFluxDiff ) );
		}else{	// OML Electron flux, SOML Ion Flux
			double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
			double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
			IonFluxDiff = Pdata->IonDensity*
				IonThermalVelocity*(1.0/(4.0*sqrt(2.0)*uz))*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm));
			
			ElectronFluxDiff = -Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
			dT = Richardson*pow(Sample->get_temperature(),2)*(1.0-guess
					*(Pdata->ElectronTemp/Sample->get_temperature()))
					*exp((-echarge*Sample->get_workfunction()+guess*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge;
			dTdiff = (Pdata->ElectronTemp/Sample->get_temperature())*dT-
					Richardson*pow(Sample->get_temperature(),2)
					*(Pdata->ElectronTemp/Sample->get_temperature())
					*exp((-echarge*Sample->get_workfunction()+guess*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge;
			x1 = guess - ( ((1.0-DeltaSec())*OMLElectronFlux(guess) - dT - SOMLIonFlux(guess))
				/((1.0-DeltaSec())*ElectronFluxDiff - dTdiff - IonFluxDiff ) );
		}

		i ++;
	}while( fabs(guess-x1) > Accuracy );
	return guess;
}

// Following Semi-empirical fit to Sceptic results as detailed in Chris Willis Thesis,
// https://spiral.imperial.ac.uk/handle/10044/1/9329, pages 68-70
double ChargingModel::solveCW(double guess){
	C_Debug("\tIn ChargingModel::solveCW()\n\n");
	double A = Pdata->mi;
	double b = Pdata->IonTemp/Pdata->ElectronTemp;
	double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));
	double Rho = Sample->get_radius()/DebyeLength;
	double Rho_OML = 1.25*pow(b,0.4);
	double Rho_Upper = 50.0;
	double Potential(0.0);
	if( Rho <= Rho_OML ){ // This is the OML Limit
		if( b <=2 ){ 	// Ti <= 2.0*Te
			Potential = 0.405*log(Pdata->A)+(0.253+0.021*log(Pdata->A))*log(b)+2.454;
		}else{ 			// Ti > 2.0*Te
			Potential = 0.401*log(Pdata->A)+(-0.122+0.029*log(Pdata->A))*log(b)+2.698;
		}
	}else if( Rho <= Rho_Upper && Rho > Rho_OML ){ // This is the transition region
		double Gradient = (log(Rho/Rho_Upper)/log(Rho_Upper/Rho_OML))+1.0;
		double DeltaPhi = 0.5*log(2.0*PI*(Me/Pdata->mi)*(1.0+5.0*b/3.0))*Gradient;
		int i(0);
		do{
			if( i > 0 ) // If first time in loop, don't update guess. Avoids numerical instabilities.
				guess = Potential;
			Potential = guess - ( ( exp(-guess) - sqrt(b*(Me/Pdata->mi))*(1+guess/b-DeltaPhi/b) ) /( -exp(-guess) - sqrt((Me/Pdata->mi)/b) ) );
			i ++;
		}while( fabs(guess-Potential) > Accuracy );
		return guess;
	}else if( Rho > Rho_Upper ){ // This is the MOML limit
		if( b <=2 ){ 	// Ti <= 2.0*Te
			Potential = 0.456*log(Pdata->A)+3.179;
		}else{ 			// Ti > 2.0*Te
			Potential = 0.557*log(Pdata->A)-(0.386+0.024*log(Pdata->A))*log(b)+3.399;
		}
	}else{
		std::cerr << "\nError in ChargingModel::solveCW()! Rho is poorly defined!\n\n";
	}
	return Potential;
}

double ChargingModel::solveOML(double a, double guess){
	C_Debug("\tIn ChargingModel::solveOML(double a, double guess)\n\n");
	if( a >= 1.0 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"DeltaTot >= 1.0 in solveOML(). DeltaTot being set equal to unity.\n");
		a = 1.0;
	}
	double b = Pdata->IonTemp/Pdata->ElectronTemp;

	double C = Me/Pdata->mi;

	double x1(0.0);
	int i(0);
	do{
		if( i > 0 )
			guess = x1;
		x1 = guess - ( ( (1-a)*exp(-guess) - sqrt(b*C)*(1+guess/b) ) /( (a-1)*exp(-guess) - sqrt(C/b) ) );
		i ++;
	}while( fabs(guess-x1) > Accuracy );
	return guess;
}

double ChargingModel::solvePosOML(double a, double guess){
	C_Debug("\tIn ChargingModel::solvePosOML(double a, double guess)\n\n");
	if( a >= 1.0 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"DeltaTot >= 1.0. DeltaTot being set equal to unity.");
		a = 1.0;
	}
	double b = Pdata->IonTemp/Pdata->ElectronTemp;

	double C = Me/Pdata->mi;

	double x1 = guess - ( (( 1-a)*(1+guess) - sqrt(b*C)*exp(-guess))/((1-a) + sqrt(b*C)*exp(-guess) ) );
	
	while( fabs(guess-x1) > Accuracy ){
		guess = x1;
		x1 = guess - ( (( 1-a)*(1+guess) - sqrt(b*C)*exp(-guess))/((1-a) + sqrt(b*C)*exp(-guess) ) );
	}
	return guess;
}

// Solve the Modified orbital motion limited potential for large dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solveMOML(){
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double MassRatio = Pdata->mi/Me;
	double HeatCapacityRatio = 1.0;

	double Arg = sqrt(2*PI*TemperatureRatio*(1+HeatCapacityRatio*TemperatureRatio))
    		*exp(TemperatureRatio);
	if( std::isinf(exp(TemperatureRatio)) ){
		static bool runOnce = true;
		WarnOnce(runOnce,"LambertW exp(TemperatureRatio) is infinite in solveMOML()! Assuming Potential = 0.0");
		return 0;
	}else if( std::isinf(Arg) || Arg < 0.0 ){
		std::cout << "\nLambertW(Arg == " << Arg << ")!";
		throw LambertWFailure();
		return 0;
	}

    return -1.0*(TemperatureRatio
    		-LambertW(Arg)+log(2*PI*(1+HeatCapacityRatio*TemperatureRatio)/MassRatio)/2);
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
	double Delta_Phi_em = 0.5*log((2.0*PI/MassRatio)*(1.0+HeatCapacityRatio*TemperatureRatio)/pow(1.0-DeltaTot,2.0));
    // Uncomment following line to compare MOMLWEM results with figure (1) of paper:
    // N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
    //std::cout << DeltaTot << "\t" << Delta_Phi_em << "\n";
    double Arg = (1.0-DeltaTot)*sqrt(MassRatio*TemperatureRatio)
                    *exp(TemperatureRatio/Ionization+Delta_Phi_em/Ionization)/Ionization;
	if( std::isinf(exp(TemperatureRatio/Ionization+Delta_Phi_em/Ionization)) ){
		static bool runOnce = true;
		WarnOnce(runOnce,"LambertW exp(TemperatureRatio/Ionization+Delta_Phi_em/Ionization) is infinite in solveMOMLWEM()! Assuming Potential = 0.0");
		return 0;
	}else if( std::isinf(Arg) || Arg < 0.0 ){
		std::cout << "\nLambertW(Arg == " << Arg << ")!";
		throw LambertWFailure();
		return 0;
	}
    return -1.0*(TemperatureRatio/Ionization+Delta_Phi_em/Ionization-LambertW(Arg));
}

// Solve the Modified orbital motion limited potential for large emitting dust grains.
// See the paper by Minas and Nikoleta, equation (10) using (1), (2), (5) & (9)
// N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
double ChargingModel::solveMOMLEM(){
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double MassRatio = Pdata->mi/Me;
	double Chi = 2.0*Sample->get_temperature()/(Pdata->ElectronTemp);
	double Delta = 0.01; // Ratio of emitted electron density to sheath edge density
	double PotentialWellDiff = solveDeltaMOMLEM(Tau, MassRatio, Chi, Delta);
	double DeltaTot = 0.0;
	if( PotentialWellDiff == 0.0 ){ // No potential well has formed!
		DeltaTot = Chi*sqrt(Delta)/PotentialWellDiff;
	}else{ // Potential well has formed!
		DeltaTot = Chi*sqrt(Delta)*PotentialWellDiff;
	}

	double HeatCapacityRatio = 1.0;
	double Ionization = Pdata->Z; 	// Ionization
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double PlasmaFlowSpeed = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	double Delta_Phi_em = 0.5*log((2.0*PI/MassRatio)*(1.0+HeatCapacityRatio*Tau)/pow(1.0-DeltaTot,2.0));
    // Uncomment following line to compare MOMLWEM results with figure (1) of paper:
    // N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
    //std::cout << DeltaTot << "\t" << Delta_Phi_em << "\n";
    double Arg = (1.0-DeltaTot)*sqrt(MassRatio*Tau)
                    *exp(Tau/Ionization+Delta_Phi_em/Ionization)/Ionization;
	if( std::isinf(exp(Tau/Ionization+Delta_Phi_em/Ionization)) ){
		static bool runOnce = true;
		WarnOnce(runOnce,"LambertW exp(Tau/Ionization+Delta_Phi_em/Ionization) is infinite in solveMOMLWEM()! Assuming Potential = 0.0");
		return 0;
	}else if( std::isinf(Arg) || Arg < 0.0 ){
		std::cout << "\nLambertW(Arg == " << Arg << ")!";
		throw LambertWFailure();
		return 0;
	}
    return -1.0*(Tau/Ionization+Delta_Phi_em/Ionization-LambertW(Arg));
}

double ChargingModel::DeltaTherm()const{
	C_Debug("\tIn ChargingModel::DeltaTherm()\n\n");

	double dtherm = (Richardson*pow(Sample->get_temperature(),2)*exp(-(Sample->get_workfunction()*echarge)
					/(Kb*Sample->get_temperature())))/(echarge*OMLElectronFlux(Sample->get_potential()));
	assert(dtherm >= 0.0 && dtherm == dtherm && dtherm != INFINITY );
	return dtherm;
}

double ChargingModel::ThermFluxSchottky( double Potential )const{
	C_Debug("\tIn ChargingModel::DeltaTherm()\n\n");
	// Returns the flux of electrons due to Thermionic emission
	// Following the Richard-Dushmann formula with Schottky Correction.
	// Negative Dust grains have the normal Richard-Dushmann form, dependant on work funciton
	// Positive Dust grains are the same flux multiplied by (1-phi)e^(-e phi/(kb Td))
	if( Potential >= 0.0 ){
		return Richardson*pow(Sample->get_temperature(),2)
					*exp(-echarge*Sample->get_workfunction()/(Kb*Sample->get_temperature()))/echarge;
	}else{

		return Richardson*pow(Sample->get_temperature(),2)*(1.0-Potential
					*(Pdata->ElectronTemp/Sample->get_temperature()))
					*exp((-echarge*Sample->get_workfunction()-Potential*Kb*Pdata->ElectronTemp)
					/(Kb*Sample->get_temperature()))/echarge;
	}
}

double ChargingModel::DeltaSec()const{
	C_Debug("\tIn ChargingModel::DeltaSec()\n\n");
	double ConvertKtoev(8.6173303e-5);
	return sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem());
}
