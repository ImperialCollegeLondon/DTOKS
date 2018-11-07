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
	if( UseModel[0] ){	// SOML Charging model
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		Potential = solveSOML(DSec+DTherm, Sample->get_potential());

	}else if( UseModel[1] ){	// SMOML Charging model
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		Potential = solveSMOML(DSec+DTherm, Sample->get_potential());

	}else if( UseModel[2] ){	// In this case, use Hutchinson, Patterchini and Lapenta
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		Potential = solvePHL(DSec+DTherm, Sample->get_potential());

	}else if( UseModel[3] ){	// Dynamically choose charging model
		double DebyeRatio = Sample->get_radius()/sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));
		double Betae = Sample->get_radius()
			/(sqrt(PI*Pdata->ElectronTemp*Me)/(2.0*echarge*echarge*Pdata->MagneticField*Pdata->MagneticField));

		if( DebyeRatio <= 0.1 && Betae <= 0.01 ){ // Unmagnetised, small dust grains
			Potential = solveSOML(0.0, Sample->get_potential());
		}else if( DebyeRatio > 0.1 && Betae <= 0.01 ){ // Unmagnetised, large dust grains
			Potential = solveSMOML(0.0, Sample->get_potential());
		}else if( DebyeRatio <= 0.1 && Betae > 0.01 ){ // Magnetised, small dust grains
			Potential = solvePHL(0.0, Sample->get_potential());
		}else if( DebyeRatio <= 0.1 && Betae > 0.01 ){ // Magnetised, large dust grains!!!
			Potential = solvePHL(0.0, Sample->get_potential());
		}
	}else if( UseModel[4] ){	// OML Charging model
		Potential = solveOML(0.0,Potential);
		if( Potential < 0 ){ // In this case we have a positive grain!
			Potential = solvePosOML( 0.0, Sample->get_potential());
		}
	}else if( UseModel[5] ){	// MOML Charging model
		Potential = solveMOML();
		if( Potential < 0 ){ // In this case we have a positive grain!
			Potential = solvePosOML( 0.0, Sample->get_potential());
		}
	}else if( UseModel[6] ){ // Original DTOKS Charging scheme
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
	}else if( UseModel[7] ){ // In this case, maintain a potential well for entire temperature range
		Potential = solveOML( 0.0, Sample->get_potential()) 
			- Kb*Sample->get_temperature()/(Pdata->ElectronTemp*Kb);
		if( Potential < 0 ){ // In this case we have a positive grain! Assume well disappears
			Potential = solvePosOML( 0.0, Sample->get_potential());
		}
	}else if( UseModel[8] ){ // 
		// Assume the grain is negative and calculate potential
		DSec = DeltaSec();
		DTherm = DeltaTherm();
		if ( (DSec + DTherm) >= 1.0 ){ // If dust grain is actually positive, (negative normalised potential) recalculate it...
			Potential = solvePosSchottkyOML(); // Quasi-neutrality assumed, we take ni as ni=ne
		}else{
			Potential = solveNegSchottkyOML(0.0);
			if( Potential < 0.0 ){
				Potential = solvePosSchottkyOML();
			}
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
double ChargingModel::solveSOML(double a, double guess){
	C_Debug("\tIn ChargingModel::solveSOML(double a, double guess)\n\n");

	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;

	
	double IonFluxDiff(0.0);
	if( guess >= 0.0 ){
		IonFluxDiff = (1.0/4.0)*Pdata->IonDensity*Pdata->Z*
			(IonThermalVelocity/uz)*Pdata->Z*erf(uz)/Tau;
	}else{
		double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
		double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
		IonFluxDiff = Pdata->IonDensity*Pdata->Z*
			IonThermalVelocity*(1.0/(16.0*uz))*(2.0*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm))
				+(1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp))
				+(1.0/sqrt(PI))*sqrt(Pdata->Z/(Tau*guess))*(2.0*uz*uz+2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm)));
			//std::cout << "\n" << (1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp)) << "\t" 
			//	<< (1.0/sqrt(PI))*sqrt(Pdata->Z/(Tau*guess))*(2.0*uz*uz+2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm));
	}

	double x1 = guess - ( ((1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
	//std::cout << "\n" << guess << "\t" << x1 << "\t" << guess << "\t" << OMLElectronFlux(guess) << "\t" << SOMLIonFlux(guess) << "\t" << IonFluxDiff;
	
	while( fabs(guess-x1) > Accuracy ){
		guess = x1;

		//std::cout << "\n" << guess << "\t" << x1 << "\t" << guess << "\t" << OMLElectronFlux(guess) << "\t" << SOMLIonFlux(guess) << "\t" << IonFluxDiff;
		if( guess >= 0.0 ){
			IonFluxDiff = (1.0/4.0)*Pdata->IonDensity*Pdata->Z*
				(IonThermalVelocity/uz)*Pdata->Z*erf(uz)/Tau;
		}else{
			double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
			double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
			IonFluxDiff = Pdata->IonDensity*Pdata->Z*
				IonThermalVelocity*(1.0/(4.0*uz))*(2.0*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm))
					+(1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp))
					+sqrt(Pdata->Z/(PI*Tau*guess))*(2.0*uz*uz-2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm)));
			//std::cout << "\n" << (1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp)) << "\t" 
			//	<< (1.0/sqrt(PI))*sqrt(Pdata->Z/(Tau*guess))*(2.0*uz*uz+2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm));
		}

		//std::cin.get();
		x1 = guess - ( (( 1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
	}
	return guess;
}

// Solve the shifted orbital motion limited potential for samll dust grains.
// See drews Thesis, pg 55
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solveSMOML(double a, double guess){
	C_Debug("\tIn ChargingModel::solveSMOML(double a, double guess)\n\n");

	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	
	
	double IonFluxDiff(0.0);
	double x1(0.0);

	if( guess >= 0.0 ){ // OML Electron flux, SMOML Ion Flux
		IonFluxDiff = Pdata->IonDensity*Pdata->Z*(IonThermalVelocity/sqrt(2.0*PI))*sqrt(PI)*erf(uz)/(Tau*2.0*uz);
		x1 = guess - ( ((1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SMOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
	}else{ // OML Electron flux, SOML Ion Flux
		double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
		double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
		IonFluxDiff = Pdata->IonDensity*Pdata->Z*IonThermalVelocity*(1.0/(16.0*uz))*(2.0*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm))
				+(1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp))
				+(1.0/sqrt(PI))*sqrt(Pdata->Z/(Tau*guess))*(2.0*uz*uz+2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm)));
		x1 = guess - ( ((1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SMOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
	}
	
	while( fabs(guess-x1) > Accuracy ){
		guess = x1;

		if( guess >= 0.0 ){  // PHL Electron flux, SMOML Ion Flux
			IonFluxDiff = Pdata->IonDensity*Pdata->Z*(IonThermalVelocity/sqrt(2.0*PI))*sqrt(PI)*erf(uz)/(Tau*2.0*uz);
			x1 = guess - ( ((1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SMOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
		}else{	// OML Electron flux, SOML Ion Flux
			double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
			double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
			IonFluxDiff = Pdata->IonDensity*Pdata->Z*IonThermalVelocity*(1.0/(16.0*uz))*(2.0*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm))
					+(1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp))
					+(1.0/sqrt(PI))*sqrt(Pdata->Z/(Tau*guess))*(2.0*uz*uz+2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm)));
			x1 = guess - ( ((1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SMOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
		}
	}
	return guess;
}

// Solve the Potential for a sphere in a collisionless magnetoplasma following
// L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).
// Assuming that Lambda_s >> r_p, i.e dust grain is small.
double ChargingModel::solvePHL(double a, double guess){
	C_Debug("\tIn ChargingModel::solvePHL(double a, double guess)\n\n");
	double x1(0.0);

	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;

	double Beta = Sample->get_radius()
			/(sqrt(PI*Pdata->ElectronTemp*Me)/(2.0*echarge*echarge*Pdata->MagneticField*Pdata->MagneticField));
	
	double AtomicNumber = Pdata->Z;	
	double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));

	double z = Beta/(1+Beta);
	double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+1.150*z*z*z*z*z;
	double eta = -(guess/Beta)*(1.0+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
	double d_eta = -(1.0/Beta)*(1.0+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
	
	double w = (eta/(1.0+eta));
	double d_w = (d_eta/(1.0+eta)-eta*d_eta/pow(1.0+eta,2.0));
	double d_A = 0.678*d_w+1.543*2.0*w*d_w-1.212*3.0*w*w*d_w;

	double IonFluxDiff(0.0);


	if( guess >= 0.0 ){ // PHL Electron flux, SOML Ion Flux
		IonFluxDiff = (1.0/4.0)*Pdata->IonDensity*Pdata->Z*(IonThermalVelocity/uz)*Pdata->Z*erf(uz)/Tau;
		double ElectronFluxDiff=Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(d_A-d_A*i_star)-PHLElectronFlux(guess);
		x1 = guess - ( ((1.0-a)*PHLElectronFlux(guess) - Pdata->Z*SOMLIonFlux(guess))/((1.0-a)*ElectronFluxDiff - IonFluxDiff ) );
	}else{ // OML Electron flux, SOML Ion Flux
		double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
		double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
		IonFluxDiff = Pdata->IonDensity*Pdata->Z*
			IonThermalVelocity*(1.0/(16.0*uz))*(2.0*(Pdata->Z/Tau)*(erf(uzp)+erf(uzm))
				+(1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp))
				+(1.0/sqrt(PI))*sqrt(Pdata->Z/(Tau*guess))*(2.0*uz*uz+2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm)));
		x1 = guess - ( ((1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
	}

	while( fabs(guess-x1) > Accuracy ){
		guess = x1;

		if( guess >= 0.0 ){  // PHL Electron flux, SOML Ion Flux
			eta = -(guess/Beta)*(1.0+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));
			w = (eta/(1.0+eta));
			d_w = (d_eta/(1.0+eta)-eta*d_eta/pow(1.0+eta,2.0));
			d_A = 0.678*d_w+1.543*2.0*w*d_w-1.212*3.0*w*w*d_w;
			double ElectronFluxDiff=Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(d_A-d_A*i_star)-PHLElectronFlux(guess);
			IonFluxDiff = (1.0/4.0)*Pdata->IonDensity*Pdata->Z*(IonThermalVelocity/uz)*Pdata->Z*erf(uz)/Tau;
			x1 = guess - ( ((1.0-a)*PHLElectronFlux(guess) - Pdata->Z*SOMLIonFlux(guess))/((1.0-a)*ElectronFluxDiff - IonFluxDiff ) );
		}else{	// OML Electron flux, SOML Ion Flux
			double uzp = uz+sqrt(-Pdata->Z*guess/Tau);
			double uzm = uz-sqrt(-Pdata->Z*guess/Tau);
			IonFluxDiff = Pdata->IonDensity*Pdata->Z*
				IonThermalVelocity*(1.0/(16.0*uz))*(2.0*Pdata->Z/Tau*(erf(uzp)+erf(uzm))
					+(1.0+2.0*(uz*uz+Pdata->Z*guess/Tau))*sqrt(Pdata->Z/(Tau*PI*guess))*(exp(-uzm*uzm)-exp(-uzp*uzp))
					+(1.0/sqrt(PI))*sqrt(Pdata->Z/(Tau*guess))*(2.0*uz*uz+2.0*Pdata->Z*guess/Tau+1.0)*(exp(-uzp*uzp)-exp(-uzm*uzm)));
			x1 = guess - ( ((1.0-a)*OMLElectronFlux(guess) - Pdata->Z*SOMLIonFlux(guess))/((a-1)*OMLElectronFlux(guess) - IonFluxDiff ) );
		}
	}
	return guess;
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
		WarnOnce(runOnce,"LambertW exp(TemperatureRatio) is infinite in solveSOML()! Assuming Potential = 0.0");
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
	double Delta_Phi_em = 0.5*log((2.0*PI/MassRatio)*(1+HeatCapacityRatio*TemperatureRatio)/pow(1.0-DeltaTot,2.0));
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

// Solve the orbital motion limited potential for small dust grains accounting for electron emission.
// WARNING: This is only valid for positively charged dust grains.
// See drews Thesis, pg 52
// https://spiral.imperial.ac.uk/bitstream/10044/1/32003/1/Thomas-D-2016-PhD-thesis.pdf
double ChargingModel::solvePosSchottkyOML(){
	C_Debug("\tIn ChargingModel::solvePosSchottkyOML()\n\n");
	double Wf = echarge*Sample->get_workfunction();
	double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
	double Ionization = Pdata->Z; 		// Ionization
	double vi = sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi));
	double ve = sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
	double Arg = TemperatureRatio*exp(TemperatureRatio)*(vi+Richardson*pow(Sample->get_temperature(),2)
			*exp(-Wf/(Kb*Sample->get_temperature()))/(Ionization*echarge*Pdata->IonDensity))/(ve*(1-DeltaSec()));

	if( std::isinf(exp(TemperatureRatio)) ){
		static bool runOnce = true;
		WarnOnce(runOnce,"LambertW exp(TemperatureRatio) is infinite in solvePosSchottkyOML()! Assuming Potential = 0.0");
		return 0;
	}else if( std::isinf(Arg) || Arg < 0.0 ){
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
	C_Debug("\tIn ChargingModel::solveNegSchottkyOML()\n\n");

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

double ChargingModel::DeltaTherm()const{
	C_Debug("\tIn ChargingModel::DeltaTherm()\n\n");

	double dtherm = (Richardson*pow(Sample->get_temperature(),2)*exp(-(Sample->get_workfunction()*echarge)
					/(Kb*Sample->get_temperature())))/(echarge*OMLElectronFlux(Sample->get_potential()));

	assert(dtherm >= 0.0 && dtherm == dtherm && dtherm != INFINITY );
	return dtherm;
}

double ChargingModel::DeltaSec()const{
	C_Debug("\tIn ChargingModel::DeltaSec()\n\n");
	double ConvertKtoev(8.6173303e-5);
	return sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem());
}
