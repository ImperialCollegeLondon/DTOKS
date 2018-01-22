//#define PAUSE
//#define FORCE_DEBUG
//#define FORCE_DEEP_DEBUG

//#include <math.h>
#include "ForceModel.h"
#include "MathHeader.h" // This is weird

// Default Constructor, no arguments
ForceModel::ForceModel():Model(){
	F_Debug("\n\nIn ForceModel::ForceModel():Model()\n\n");
	UseModel = {false,false,false,false,false,false};
	CreateFile("Default_Force_filename.txt");
}

ForceModel::ForceModel(std::string filename, float accuracy, std::array<bool,NumModels> models, 
			Matter *& sample, PlasmaData *& pdata) : Model(sample,pdata,accuracy){
	F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, float accuracy, std::array<bool,3> models, Matter *& sample, PlasmaData const *& pdata) : Model(sample,pdata,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, float accuracy, std::array<bool,NumModels> models, 
			Matter *& sample, PlasmaGrid & pgrid) : Model(sample,pgrid,accuracy){
	F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, float accuracy, std::array<bool,3> models, Matter *& sample, PlasmaGrid const& pgrid) : Model(sample,pgrid,accuracy)\n\n");
	UseModel = models;
	CreateFile(filename);
}


void ForceModel::CreateFile(std::string filename){
	F_Debug("\tIn ForceModel::CreateFile(std::string filename)\n\n");
	FileName=filename;
	ModelDataFile.open(FileName);
	ModelDataFile << std::fixed << std::setprecision(16) << std::endl;
	ModelDataFile << "Time\tPosition\tVelocity\tRotationFreq";
	bool PrintGravity = true; // Lol
	if( UseModel[0] && PrintGravity ) 	ModelDataFile << "\tGravity";
	if( UseModel[1] ) 			ModelDataFile << "\tCentrifugal";
	if( UseModel[2] ) 			ModelDataFile << "\tLorentz";
	if( UseModel[3] ) 			ModelDataFile << "\tIonDrag";
	if( UseModel[4] ) 			ModelDataFile << "\tHybridIonDrag";
	if( UseModel[5] ) 			ModelDataFile << "\tNeutralDrag";
	
	ModelDataFile << "\n";
	ModelDataFile.close();
	ModelDataFile.clear();
	Print();
}

void ForceModel::Print(){
	F_Debug("\tIn ForceModel::Print()\n\n");
	ModelDataFile.open(FileName,std::ofstream::app);
	ModelDataFile << TotalTime << "\t" 
			<< Sample->get_position() << "\t" << Sample->get_velocity() << "\t" << Sample->get_rotationalfreq();
	bool PrintGravity = true; // Lol
	if( UseModel[0] && PrintGravity ) 	ModelDataFile << "\t" << Gravity(); 
	if( UseModel[1] ) 			ModelDataFile << "\t" << Centrifugal();
	if( UseModel[2] ) 			ModelDataFile << "\t" << LorentzForce();
	if( UseModel[3] ) 			ModelDataFile << "\t" << DTOKSIonDrag();
	if( UseModel[4] ) 			ModelDataFile << "\t" << HybridIonDrag();
	if( UseModel[5] ) 			ModelDataFile << "\t" << NeutralDrag();
	ModelDataFile << "\n";
	ModelDataFile.close();
	ModelDataFile.clear();
}

double ForceModel::ProbeTimeStep()const{
	F_Debug( "\tIn ForceModel::ProbeTimeStep()const\n\n" );
	// Deal with case where power/time step causes large temperature change.
	double timestep(0);
	threevector Acceleration = CalculateAcceleration();
	// For Accuracy = 1.0, requires change in velocity less than 0.01m or 10cm/s
	if( Acceleration.mag3() == 0 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"Zero Acceleration!\ntimestep being set to unity");
		timestep = 1;	// Set arbitarily large time step (Should this be float max?)
	}else{
		timestep = (0.01*Accuracy)*(1.0/Acceleration.mag3());
	}

	// Check if the timestep should be shortened such that particles don't cross many grid cells in a single step
	// (This is often the case without this condition.)
	if( Sample->get_velocity().mag3() != 0.0 && ( get_dlx()*Accuracy/(2*Sample->get_velocity().mag3()) ) < timestep ){
		F_Debug("\ntimestep limited by grid size!");
		timestep = get_dlx()*Accuracy/(2*Sample->get_velocity().mag3());
	}
	
	// Check if the timestep is limited by the gyration of the particle in a magnetic field.
	double GyromotionTimeStep = 0.1*pow(Sample->get_radius(),2)*Sample->get_density()
			*sqrt(1-(Pdata->MagneticField.getunit()*Sample->get_velocity().getunit()))
			/(3.0*epsilon0*Pdata->ElectronTemp*fabs(Sample->get_potential())*Pdata->MagneticField.mag3());

/*	if( GyromotionTimeStep < timestep && GyromotionTimeStep > 0.0 ){
		std::cout << "\ntimestep limited by magnetic field (Gyromotion)\n";
		timestep = GyromotionTimeStep;
	}*/

	F1_Debug( "\t\tAcceleration = " << Acceleration << "\n\t\ttimestep = " << timestep << "\n");
	assert(timestep == timestep);
	assert(timestep > 0);
	return timestep;
}

double ForceModel::UpdateTimeStep(){
	F_Debug( "\tIn ForceModel::UpdateTimeStep()\n\n" );
	TimeStep = ProbeTimeStep();
	
	return TimeStep;
}

// Move the dust grain by calculating the forces acting on it
void ForceModel::Force(double timestep){
	F_Debug("\tIn ForceModel::Force(double timestep)\n\n");

	// Make sure timestep input time is valid. Shouldn't exceed the timescale of the process.
	assert(timestep > 0 && timestep <= TimeStep );

	threevector Acceleration = CalculateAcceleration();

	threevector ChangeInPosition(
					Sample->get_velocity().getx()*timestep,
					(Sample->get_velocity().gety()*timestep)/Sample->get_position().getx(),
					Sample->get_velocity().getz()*timestep);
	// WARNING! If Sample->get_position().getx() == 0, then we will get nans in the position element
	if( Sample->get_position().getx() == 0.0 ){
		ChangeInPosition.sety(0.0);
	}

	threevector ChangeInVelocity = Acceleration*timestep;

	// Assert change in absolute vel less than ten times accuracy
	assert( ChangeInVelocity.mag3() < 0.1*Accuracy ); // Assert change in velocity less than 10* TimeStep accuracy

	// Krasheninnikov, S. I. (2006). On dust spin up in uniform magnetized plasma. Physics of Plasmas, 13(11), 2004–2007.
//	double TimeOfSpinUp = Sample->get_radius()*sqrt(Mp/(Kb*Pdata->IonTemp))*Sample->get_density()/(Mp*Pdata->IonDensity);
//	TimeOfSpinUp = 1;
//	double RotationalSpeedUp = timestep*sqrt(Kb*Pdata->IonTemp/Mp)/(TimeOfSpinUp*Sample->get_radius());
//	if( Sample->get_radius() < sqrt(Kb*Pdata->IonTemp*Mp)/(echarge*Pdata->MagneticField.mag3()) ){
//	      	RotationalSpeedUp = (timestep*echarge*Pdata->MagneticField.mag3()/(TimeOfSpinUp*Mp));
//		F1_Debug( "\nREGIME TWO!" );
//	}
//	F1_Debug( "\nRho_{Ti} = " << sqrt(Kb*Pdata->IonTemp*Mp)/(echarge*Pdata->MagneticField.mag3()) <<
//		 "\nV_{Ti} = " << sqrt(Kb*Pdata->IonTemp/Mp) << "\nOmega_{i} = " 
//		<< echarge*Pdata->MagneticField.mag3()/Mp << "\ntimestep = " << timestep );
//	F1_Debug( "\nBfield.mag = " << Pdata->MagneticField.mag3() << "\nTimeOfSpinUp = " << TimeOfSpinUp
//		 << "\nSpeedUp = " << RotationalSpeedUp << "\n" );


//	Introduced on 11/10/17, this is informed by over two months work on the theory and dynamics of dust rotation
//	due to ion collection. Even with the full theory, we find values of B that are too small for tokamak conditions
//	to lead to dust breakup. We need to find a way to increase the theoretical breakup speed.
//	double B = (5*sqrt(2*PI)*Pdata->IonDensity*Mp*sqrt((Kb*Pdata->IonTemp)/Mp)*pow(Sample->get_radius(),2))
//			/(2*Sample->get_mass()); 
	double B = 1.0;
	double RotationalSpeedUp = timestep*B*(2*(echarge*(Pdata->MagneticField.mag3())/Mp)-Sample->get_rotationalfreq());

	Sample->update_motion(ChangeInPosition,ChangeInVelocity,RotationalSpeedUp);

	F1_Debug( "\nChangeInPosition : " << ChangeInPosition << "\nChangeInVelocity : " << ChangeInVelocity << "\nAcceleration : " << Acceleration << "\nTimeStep : " << TimeStep << "\n");
	F_Debug("\t"); Print();
	TotalTime += timestep;
}


// Move the dust grain by calculating the forces acting on it
void ForceModel::Force(){
	F_Debug("\tIn ForceModel::Force()\n\n");
	Force(TimeStep);
}

threevector ForceModel::CalculateAcceleration()const{
	F_Debug("\tIn ForceModel::CalculateAcceleration()const\n\n");
	// Forces: Lorentz + ion drag + gravity
	threevector Accel(0.0,0.0,0.0);
	if( UseModel[0] ) Accel += Gravity();
	if( UseModel[1] ) Accel += Centrifugal(); // Centrifugal terms. CHECK THESE TWO LINES AT SOME POINT
	if( UseModel[2] ) Accel += LorentzForce();
	if( UseModel[3] ) Accel += DTOKSIonDrag();
	if( UseModel[4] ) Accel += HybridIonDrag();
	if( UseModel[5] ) Accel += NeutralDrag();

	F1_Debug( "\n\t\tg = " << Gravity() );
	F1_Debug( "\n\t\tcentrifugal = " << Centrifugal() );
	F1_Debug( "\n\t\tlorentzforce = " << LorentzForce() );
	F1_Debug( "\n\t\tDTOKSIonDrag = " << DTOKSIonDrag() );
	F1_Debug( "\n\t\tHybridIonDrag = " << HybridIonDrag() );
	F1_Debug( "\n\t\tNeutralDrag = " << NeutralDrag() );
	F1_Debug( "\n\t\tAccel = " << Accel << "\n\n" );

	return Accel;
}

threevector ForceModel::Gravity()const{
	threevector gravity(0.0,0.0,0.0);

	double Theta = Sample->get_position().gety();
	// Setup for Magnum-PSI
	gravity.setx(-Pdata->Gravity.mag3()*cos(Theta));
	gravity.sety(Pdata->Gravity.mag3()*sin(Theta));
	return gravity;
}

threevector ForceModel::HybridIonDrag()const{
	double IonThermalVelocity = sqrt(Kb*Pdata->IonTemp/Mp);
	double u = Pdata->PlasmaVel.mag3()*(1.0/IonThermalVelocity); 	// Normalised ion flow velocity
	double Tau = Pdata->ElectronTemp/Pdata->IonTemp;

	if( Pdata->IonTemp == 0.0 || u == 0.0 ){
		threevector Zero(0.0,0.0,0.0);
		return Zero;
	}

	double z = Sample->get_potential()*4.0*PI*epsilon0;
	double CoulombLogarithm = 17.0;		// Approximation of coulomb logarithm	


	double Coefficient = sqrt(2*PI)*pow(Sample->get_radius(),2.0)*Pdata->IonDensity*Mp*Pdata->PlasmaVel.square();
	double term1 = sqrt(PI/2.0)*erf(u/sqrt(2))*(1.0+u*u+(1.0-(1.0/(u*u)))*(1.0+2*Tau*z)+4*z*z*Tau*Tau*CoulombLogarithm/(u*u));
	double term2 = (1.0/u)*exp(-u*u/2.0)*(1.0+2.0*Tau*z+u*u-4*z*z*Tau*Tau*CoulombLogarithm);	
	
	threevector HybridDrag = Coefficient*(term1+term2)*(Pdata->PlasmaVel-Sample->get_velocity()).getunit();
	return HybridDrag*(1.0/Sample->get_mass());
}

// Calculations for ion drag: Mach number, shielding length with fitting function and thermal scattering parameter
threevector ForceModel::DTOKSIonDrag()const{
	F_Debug("\tIn ForceModel::DTOKSIonDrag()\n\n");
	threevector Fid(0.0,0.0,0.0);

// ION TEMPERATURE IN THIS FUNCTION IS IN ev.
	double ConvertKelvsToeV(8.621738e-5);
	threevector Mt(0.0,0.0,0.0);
	if( Pdata->IonTemp != 0 ) Mt = (Pdata->PlasmaVel-Sample->get_velocity())*sqrt(Mp/(Kb*Pdata->IonTemp)); 
	F_Debug("\nMt = " << Mt);

        if( Pdata->IonDensity == 0 || Pdata->IonTemp == 0 || Pdata->ElectronTemp == 0  || Mt.mag3() == 0 ){ 
		// || Mt.mag3() == 0 ){
                F_Debug("\nWarning! IonDensity = " << Pdata->IonDensity << "\nIonTemp = " 
			<< Pdata->IonTemp << "\nElectronTemp = " 
			<< Pdata->ElectronTemp << "!\nThis blows up calculations! Setting Fid=0.");
                Fid = threevector(0.0,0.0,0.0);
        }else{
		if(Mt.mag3()<2.0){ // Relative speed less than twice mach number, use Fortov et al theory with screening length 'Lambda'.

			double lambda = sqrt(epsilon0/(Pdata->IonDensity*echarge*exp(-Mt.mag3()*Mt.mag3()/2)
				*(1.0/(Pdata->IonTemp*ConvertKelvsToeV))+1.0/(Pdata->ElectronTemp*ConvertKelvsToeV)));
			double beta = Pdata->ElectronTemp*ConvertKelvsToeV*Sample->get_radius()
				*fabs(Sample->get_potential())/(lambda*Pdata->IonTemp*ConvertKelvsToeV);
			F_Debug("\nlambda = " << lambda << "\nbeta = " << beta << "\nPot = " << Sample->get_potential());

			if(beta>13.0) std::cout << "nonlinear drag parameter" << std::endl;

			threevector FidS(0.0,0.0,0.0);
			if( Sample->get_potential() == 0 || beta == 0 ){
				FidS = threevector(0.0,0.0,0.0);
			}else{
				double Lambda = -exp(beta/2.0)*Exponential_Integral_Ei(-beta/2.0); 
				FidS = Mt*(sqrt(32*PI)/3.0*epsilon0*pow(Pdata->IonTemp*ConvertKelvsToeV,2)*Lambda*pow(beta,2));
			}

			threevector FidC =(Pdata->PlasmaVel-Sample->get_velocity())*4.0*PI*pow(Sample->get_radius(),2)
				*Pdata->IonDensity*Mp*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*exp(-Sample->get_potential()); 
			//for John's ion drag... I assume here and in other places in the 
			//calculation that the given potential is normalised to kTe/e
	    
			//Do The same but with the ion current instead of the electron current
	    
			//.....
	    
			Fid=FidS+FidC;
			F_Debug("\nFidS = " << FidS << "\nFidC = " << FidC);
		}else{	// Relative speed greater than twice the mach number, use just plain collection area
//			double lambdadi = sqrt(epsilon0*Pdata->IonTemp*ConvertKelvsToeV)/sqrt(Pdata->IonDensity*echarge);
			//F_Debug("\nlambdadi = " << lambdadi);
			Fid = Mt*Mt.mag3()*PI*Pdata->IonTemp*ConvertKelvsToeV*pow(Sample->get_radius(),2)*Pdata->IonDensity*echarge;
		}
	}
	return Fid*(3/(4*PI*pow(Sample->get_radius(),3)*Sample->get_density()));
}

// Calculations for Neutral Drag
threevector ForceModel::NeutralDrag()const{
	F_Debug("\tIn ForceModel::DTOKSIonDrag()\n\n");
//	return (Pdata->PlasmaVel-Sample->get_velocity())*Pdata->NeutralDensity*sqrt(2*Kb*Pdata->IonTemp*Mp)*PI
//			*pow(Sample->get_radius(),2);
	return (Pdata->PlasmaVel-Sample->get_velocity())*Mp*sqrt(4*PI)*NeutralFlux()*PI*pow(Sample->get_radius(),2)*(1.0/Sample->get_mass());
}

threevector ForceModel::LorentzForce()const{
	F_Debug("\tIn ForceModel::LorentzForce()\n\n");
	// Dust grain charge to mass ratio
	double Charge = 4.0*PI*epsilon0*Sample->get_radius()*Kb*Pdata->ElectronTemp*Sample->get_potential()/echarge;
	double qtom = Charge/Sample->get_mass();

//	double ConvertKelvsToeV(8.621738e-5);
	// Estimate charge from potential difference
//	if( Sample->is_positive() ) (qtom = 3.0*epsilon0*Kb*Sample->get_temperature())
//					/(echarge*pow(Sample->get_radius(),2)*Sample->get_density());
//	else qtom = -3.0*epsilon0*Pdata->ElectronTemp*ConvertKelvsToeV*Sample->get_potential()
//			/(pow(Sample->get_radius(),2)*Sample->get_density());	
	//else qtom = -1000.0*3.0*epsilon0*V/(a*a*rho);//why it had Te???????
	//edo to eixa allaksei se ola ta runs sto After 28_Feb ara prepei na ksana ginoun
	// Google Translate: Here I had changed it to all the brides in After 28_Feb so they have to be done again
	
	threevector returnvec = (Pdata->ElectricField+(Sample->get_velocity()^Pdata->MagneticField))*qtom;
	return returnvec;
}

threevector ForceModel::Centrifugal()const{
	F_Debug("\tIn ForceModel::Centrifugal()\n\n");
	threevector returnval(
		Sample->get_velocity().gety()*Sample->get_velocity().gety()/Sample->get_position().getx(),
		-Sample->get_velocity().getx()*Sample->get_velocity().gety()/Sample->get_position().getx(),
		0.0);
	//assert( returnval == returnval );
	return returnval;
}
