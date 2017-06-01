//#define PAUSE
#define FORCE_DEBUG

#include "ForceModel.h"
#include "MathHeader.h" // This is weird

// Default Constructor, no arguments
ForceModel::ForceModel():Model(){
	F_Debug("\n\nIn ForceModel::ForceModel():Model()\n\n");
	CreateFile("Default_Force_filename.txt");
	UseModel = {false,false,false};
	TimeStep = 0;
	TotalTime = 0;
}

ForceModel::ForceModel(std::string filename, double accuracy, std::array<bool,4> models, 
			Matter *& sample, PlasmaData const& pdata) : Model(sample,pdata,accuracy){
	F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, std::array<bool,3> models, Matter *& sample, PlasmaData const& pdata) : Model(sample,pdata)\n\n");
	CreateFile(filename);
	UseModel = models;
	TimeStep = 0;
	TotalTime = 0;
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
	
	assert(TimeStep > 0);
	TotalTime += TimeStep;
}

// Move the dust grain by calculating the forces acting on it
void ForceModel::Force(){
	F_Debug("\tIn ForceModel::Force()\n\n");

//	std::cout << "\n(3)Temp = " << Sample->get_temperature() << "\nVapPressure = " << Sample->get_vapourpressure() << "\nCv = " << Sample->get_heatcapacity() << "\n";
/*
	// Forces: Lorentz + ion drag + gravity
	threevector Fid = DTOKSIonDrag();
	// Calculations for ion drag: Mach number, shielding length with fitting function and thermal scattering parameter
	// MAKE SURE THE ION DENSITY IS NOT ZERO! OTHERWISE Fid = threevector(0.0,0.0,0.0);
	if( Pdata.IonDensity == 0 ){
		std::cout << "\nError! IonDensity = " << Pdata.IonDensity << "!\nThis blows up calculations! Setting Fid=0";
		Fid = threevector(0.0,0.0,0.0);
	}

	threevector g(0.0,0.0,-9.81);
	// Centrifugal terms. MAKE SURE TO CHECK THESE TWO LINES AT SOME POINT
	F_Debug("\t"); threevector centrifugal = Centrifugal();

	threevector ChangeInPosition(
			Sample->get_position().getx()+Sample->get_velocity().getx()*TimeStep,
			Sample->get_position().gety()+(Sample->get_velocity().gety()*TimeStep)/(Sample->get_position().getx()),
			Sample->get_position().getz()+Sample->get_velocity().getz()*TimeStep);
	F_Debug("\t"); threevector LorentzForce = LorentzForce();
	threevector ChangeInVelocity = (LorentzForce + Fid + g + centrifugal)*TimeStep;
	//vd += (LorentzForce() + Fid + g)*TimeStep;
	
	//vd += centrifugal*TimeStep;

	Sample->update_motion(ChangeInPosition,ChangeInVelocity);
*/
	F_Debug("\t"); Print();
}

threevector ForceModel::DTOKSIonDrag()const{
	F_Debug("\tIn ForceModel::DTOKSIonDrag()\n\n");
	threevector Fid;

	threevector Mt = (Pdata.PlasmaVel-Sample->get_velocity());//*sqrt(Mp/echarge/Pdata.IonTemp); // FIND OUT WHAT THIS OPERATION MEANS!
	double lambda = sqrt(epsilon0/(Pdata.IonDensity*echarge))/sqrt(exp(-Mt.mag3()*Mt.mag3()/2)/Pdata.IonTemp+1.0/Pdata.ElectronTemp);
	double beta = Pdata.ElectronTemp*Sample->get_radius()*fabs(Sample->get_potential())/(lambda*Pdata.IonTemp);
	if(beta>13.0) std::cout << "nonlinear drag parameter" << std::endl;
	if(Mt.mag3()<2.0){ // Relative speed less than twice mach number, use Fortov et al theory with screening length 'Lambda'.
		double Lambda = -exp(beta/2.0)*Exponential_Integral_Ei(-beta/2.0); 
		threevector FidC,FidS;
		FidS = Mt*(sqrt(32*PI)/3.0*epsilon0*pow(Pdata.IonTemp,2)*Lambda*pow(beta,2));
		FidC =(Pdata.PlasmaVel-Sample->get_velocity())*4.0*PI*pow(Sample->get_radius(),2)*Pdata.IonDensity*Mp
			*sqrt(echarge*Pdata.ElectronTemp/2.0/PI/Me)*exp(-Sample->get_potential()); 
		//for John's ion drag... I assume here and in other places in the 
		//calculation that the given potential is normalised to kTe/e
    
		//Do The same but with the ion current instead of the electron current
    
		//.....
    
		Fid=FidS+FidC;
	}else{	// Relative speed greater than twice the mach number, use just plain collection area
		double lambdadi = sqrt(epsilon0*Pdata.IonTemp)/sqrt(Pdata.IonDensity*echarge);
		Fid = Mt*(Mt.mag3()*4*PI*epsilon0*pow(Pdata.IonTemp,2)*(pow(Sample->get_radius(),2)/(4*lambdadi*lambdadi)));	
	}
	return Fid*(3/(4*PI*pow(Sample->get_radius(),3)*Sample->get_density()));
}

threevector ForceModel::LorentzForce()const{
	F_Debug("\tIn ForceModel::LorentzForce()\n\n");
	// Dust grain charge to mass ratio
	double qtom(0);
	// Estimate charge from potential difference
	if( Sample->is_positive() ) (qtom = 3.0*epsilon0*Kb*Sample->get_temperature())
					/(echarge*pow(Sample->get_radius(),2)*Sample->get_density());
	else qtom = -3.0*epsilon0*Pdata.ElectronTemp*Sample->get_potential()/(pow(Sample->get_radius(),2)*Sample->get_density());	
	//else qtom = -1000.0*3.0*epsilon0*V/(a*a*rho);//why it had Te???????
	//edo to eixa allaksei se ola ta runs sto After 28_Feb ara prepei na ksana ginoun
	// Google Translate: Here I had changed it to all the brides in After 28_Feb so they have to be done again


	return ((Pdata.ElectricField+(Sample->get_velocity()^Pdata.MagneticField))*qtom);
}

threevector ForceModel::Centrifugal()const{
	F_Debug("\tIn ForceModel::Centrifugal()\n\n");
	threevector returnval(
		Sample->get_velocity().gety()*Sample->get_velocity().gety()/Sample->get_position().getx(),
		-Sample->get_velocity().getx()*Sample->get_velocity().gety()/Sample->get_position().getx(),
		0.0);
	return returnval;
}
