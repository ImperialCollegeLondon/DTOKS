//#define PAUSE
//#define FORCE_DEBUG
//#define FORCE_DEEP_DEBUG

#include "ForceModel.h"
#include "MathHeader.h" // This is weird

// Default Constructor, no arguments
ForceModel::ForceModel():Model(){
	F_Debug("\n\nIn ForceModel::ForceModel():Model()\n\n");
	UseModel = {false,false,false};
	TimeStep = 0;
	TotalTime = 0;
	CreateFile("Default_Force_filename.txt");
}

ForceModel::ForceModel(std::string filename, double accuracy, std::array<bool,4> models, 
			Matter *& sample, PlasmaData & pdata) : Model(sample,pdata,accuracy){
	F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, std::array<bool,3> models, Matter *& sample, PlasmaData const& pdata) : Model(sample,pdata,accuracy)\n\n");
	UseModel = models;
	TimeStep = 0;
	TotalTime = 0;
	CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, double accuracy, std::array<bool,4> models, 
			Matter *& sample, PlasmaGrid & pgrid) : Model(sample,pgrid,accuracy){
	F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, std::array<bool,3> models, Matter *& sample, PlasmaGrid const& pgrid) : Model(sample,pgrid,accuracy)\n\n");
	UseModel = models;
	TimeStep = 0;
	TotalTime = 0;
	CreateFile(filename);
}


void ForceModel::CreateFile(std::string filename){
	F_Debug("\tIn ForceModel::CreateFile(std::string filename)\n\n");
	ForceFile.open(filename);
	ForceFile << "Position\tVelocity";
	if( UseModel[0] ) ForceFile << "\tGravity";
       	if( UseModel[1] ) ForceFile << "\tCentrifugal";
       	if( UseModel[2] ) ForceFile << "\tLorentz";
       	if( UseModel[3] ) ForceFile << "\tIonDrag";

	ForceFile << "\n";
}

void ForceModel::Print(){
	F_Debug("\tIn ForceModel::Print()\n\n");
	ForceFile << Sample->get_position() << "\t" << Sample->get_velocity();
	if( UseModel[0] ) ForceFile << "\t(0.0,0.0,-9.81)"; // Maybe this should be coded better...
	if( UseModel[1] ) ForceFile << "\t" << Centrifugal();
	if( UseModel[2] ) ForceFile << "\t" << LorentzForce();
	if( UseModel[3] ) ForceFile << "\t" << DTOKSIonDrag();
	ForceFile << "\n";
}

double ForceModel::CheckTimeStep(){
	F_Debug( "\tIn ForceModel::CheckTimeStep()\n\n" );
	// Deal with case where power/time step causes large temperature change.
	
	threevector Acceleration = CalculateAcceleration();

	// For Accuracy = 1.0, requires change in velocity less than 0.01m or 10cm/s
	if( Acceleration.mag3() == 0 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"Warning! Zero Acceleration!\nTimeStep being set to unity");
		TimeStep = 1;	// Set arbitarily large time step (Should this be float max?)
	}else{
		TimeStep = (0.01*Accuracy)*(1.0/Acceleration.mag3());
	}
	F1_Debug( "\t\tAcceleration = " << Acceleration << "\n\t\tTimeStep = " << TimeStep << "\n");
	assert(TimeStep == TimeStep);
	assert(TimeStep > 0);
	return TimeStep;
}

// Move the dust grain by calculating the forces acting on it
void ForceModel::Force(){
	F_Debug("\tIn ForceModel::Force()\n\n");

	// Forces: Lorentz + ion drag + gravity
	threevector Acceleration = CalculateAcceleration();

	threevector ChangeInPosition(
			Sample->get_velocity().getx()*TimeStep,
			(Sample->get_velocity().gety()*TimeStep)/(Sample->get_position().getx()),
			Sample->get_velocity().getz()*TimeStep);

	threevector ChangeInVelocity = Acceleration*TimeStep;
	
	Sample->update_motion(ChangeInPosition,ChangeInVelocity);

	F_Debug("\t"); 
	Print();
	TotalTime += TimeStep;
}

threevector ForceModel::CalculateAcceleration()const{
	F_Debug("\tIn ForceModel::CalculateAcceleration()\n\n");
	// Forces: Lorentz + ion drag + gravity

	threevector g(0.0,0.0,-9.81);
	F_Debug("\t"); threevector Fid = DTOKSIonDrag();
        F_Debug("\t"); threevector centrifugal = Centrifugal(); // Centrifugal terms. CHECK THESE TWO LINES AT SOME POINT
        F_Debug("\t"); threevector lorentzforce = LorentzForce();
	g 		= g*UseModel[0];
	Fid 		= Fid*UseModel[1];
	centrifugal 	= centrifugal*UseModel[2];
	lorentzforce 	= lorentzforce*UseModel[3];
	F1_Debug( "\n\t\tFid = " << Fid );
	F1_Debug( "\n\t\tg = " << g );
	F1_Debug( "\n\t\tcentrifugal = " << centrifugal );
	F1_Debug( "\n\t\tlorentzforce = " << lorentzforce << "\n\n" );
        return (lorentzforce + Fid + g + centrifugal);
}

// Calculations for ion drag: Mach number, shielding length with fitting function and thermal scattering parameter
threevector ForceModel::DTOKSIonDrag()const{
	F_Debug("\tIn ForceModel::DTOKSIonDrag()\n\n");
	threevector Fid;

// ION TEMPERATURE IN THIS FUNCTION IS IN ev.
	double ConvertKelvsToeV(8.621738e-5);

	threevector Mt = (Pdata.PlasmaVel-Sample->get_velocity())*sqrt(Mp/(echarge*Pdata.IonTemp*ConvertKelvsToeV)); 
	F_Debug("\nMt = " << Mt);

        if( Pdata.IonDensity == 0 || Pdata.IonTemp*ConvertKelvsToeV == 0 || Pdata.ElectronTemp*ConvertKelvsToeV == 0 ){ 
		// || Mt.mag3() == 0 ){
                std::cerr << "\nError! IonDensity = " << Pdata.IonDensity << "\nIonTemp*ConvertKelvsToeV = " 
			<< Pdata.IonTemp*ConvertKelvsToeV << "\nElectronTemp*ConvertKelvsToeV = " 
			<< Pdata.ElectronTemp*ConvertKelvsToeV << "!\nThis blows up calculations! Setting Fid=0.";
                Fid = threevector(0.0,0.0,0.0);
        }else{
		double lambda = sqrt(epsilon0/(Pdata.IonDensity*echarge*exp(-Mt.mag3()*Mt.mag3()/2)
				*(1.0/(Pdata.IonTemp*ConvertKelvsToeV))+1.0/(Pdata.ElectronTemp*ConvertKelvsToeV)));
		double beta = Pdata.ElectronTemp*ConvertKelvsToeV*Sample->get_radius()
				*fabs(Sample->get_potential())/(lambda*Pdata.IonTemp*ConvertKelvsToeV);
		F_Debug("\nlambda = " << lambda << "\nbeta = " << beta);

		if(beta>13.0) std::cout << "nonlinear drag parameter" << std::endl;
		if(Mt.mag3()<2.0){ // Relative speed less than twice mach number, use Fortov et al theory with screening length 'Lambda'.

			assert( beta != 0 );
			double Lambda = -exp(beta/2.0)*Exponential_Integral_Ei(-beta/2.0); 
			threevector FidC,FidS;
			FidS = Mt*(sqrt(32*PI)/3.0*epsilon0*pow(Pdata.IonTemp*ConvertKelvsToeV,2)*Lambda*pow(beta,2));
			FidC =(Pdata.PlasmaVel-Sample->get_velocity())*4.0*PI*pow(Sample->get_radius(),2)*Pdata.IonDensity*Mp
				*sqrt(echarge*Pdata.ElectronTemp*ConvertKelvsToeV/(2.0*PI*Me))*exp(-Sample->get_potential()); 
			//for John's ion drag... I assume here and in other places in the 
			//calculation that the given potential is normalised to kTe/e
	    
			//Do The same but with the ion current instead of the electron current
	    
			//.....
	    
			Fid=FidS+FidC;
			F_Debug("\nLambda = " << Lambda << "FidS = " << FidS << "\nFidC = " << FidC);
		}else{	// Relative speed greater than twice the mach number, use just plain collection area
			double lambdadi = sqrt(epsilon0*Pdata.IonTemp*ConvertKelvsToeV)/sqrt(Pdata.IonDensity*echarge);
			F_Debug("\nlambdadi = " << lambdadi);
			Fid = Mt*(Mt.mag3()*4*PI*epsilon0*pow(Pdata.IonTemp*ConvertKelvsToeV,2)*(pow(Sample->get_radius(),2)
					/(4*lambdadi*lambdadi)));	
		}
	}
	return Fid*(3/(4*PI*pow(Sample->get_radius(),3)*Sample->get_density()));
}

threevector ForceModel::LorentzForce()const{
	F_Debug("\tIn ForceModel::LorentzForce()\n\n");
	// Dust grain charge to mass ratio
	double qtom(0);
	double ConvertKelvsToeV(8.621738e-5);
	// Estimate charge from potential difference

	if( Sample->is_positive() ) (qtom = 3.0*epsilon0*Kb*Sample->get_temperature())
					/(echarge*pow(Sample->get_radius(),2)*Sample->get_density());
	else qtom = -3.0*epsilon0*Pdata.ElectronTemp*ConvertKelvsToeV*Sample->get_potential()
			/(pow(Sample->get_radius(),2)*Sample->get_density());	
	//else qtom = -1000.0*3.0*epsilon0*V/(a*a*rho);//why it had Te???????
	//edo to eixa allaksei se ola ta runs sto After 28_Feb ara prepei na ksana ginoun
	// Google Translate: Here I had changed it to all the brides in After 28_Feb so they have to be done again
	threevector returnvec = (Pdata.ElectricField+(Sample->get_velocity()^Pdata.MagneticField))*qtom;
	
	
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
