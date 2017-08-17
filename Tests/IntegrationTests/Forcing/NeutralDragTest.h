#include "ForceModel.h"

int NeutralDragTest(char Element){
	clock_t begin = clock();
	// ********************************************************** //
	// FIRST, define program default behaviour

	// Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
	char EmissivityModel = 'c'; 	// Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
	char ExpansionModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
													// and (z)ero expansion
	char HeatCapacityModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
	char BoilingModel = 'y'; 	// Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
													// and (t)homson

 	// Parameters describing the heating model
	double Size=1e-6; 		// m
	double Temp=280;		// K
	double Potential = 1;		// Normalised Potential
	Matter *Sample;			// Define the sample matter type

	// Set to true all Force models that are wanted
	bool Gravity = false;		// IS OFF!
	bool Centrifugal = false;	// IS OFF!
	bool Lorentz = false;		// IS OFF!
	bool IonDrag = false;		// IS OFF!
	bool NeutralDrag = true;	// IS ON!

	PlasmaData *Pdata = new PlasmaData();
	Pdata->NeutralDensity	= 10e18;
	Pdata->NeutralTemp	= 10*1.16e4;	// K, Neutral Temperature, convert from eV
	threevector PlasmaVel(0.0,0.0,10.0);
	Pdata->PlasmaVel	= PlasmaVel;

	std::array<bool,5> ForceModels  = {Gravity,Centrifugal,Lorentz,IonDrag,NeutralDrag};

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
	threevector vinit(2.0,-4.0,5.0);
	threevector zeros(0.0,0.0,0.0);
	Sample->update_motion(zeros,vinit,0.0);
	// START NUMERICAL MODEL
	ForceModel MyModel("out_ConstantForcingTest.txt",1.0,ForceModels,Sample,Pdata);
	double Mass = Sample->get_mass();
	MyModel.UpdateTimeStep();
	size_t imax(100);

	for( size_t i(0); i < imax; i ++)
		MyModel.Force();

	threevector ModelVelocity = Sample->get_velocity();
	// END NUMERICAL MODEL

	// START ANALYTICAL MODEL
	// Factor of sqrt(4*PI) exists to match definitions between DTOKS and DUSTT of velocity
	threevector NeutralDragAccel = (Pdata->PlasmaVel-vinit)*Mp*sqrt(4*PI)*Pdata->NeutralDensity*sqrt(Kb*Pdata->NeutralTemp/(2*PI*Mp)) *PI*pow(Sample->get_radius(),2)*(1/Sample->get_mass());
	threevector AnalyticVelocity = vinit + NeutralDragAccel*(imax*MyModel.get_timestep());
	delete Pdata;
	// END ANALYTICAL MODEL

	double ReturnVal = 0;

	if( (ModelVelocity - AnalyticVelocity).mag3() == 0 ) 				ReturnVal = 1;
	else if( (( AnalyticVelocity.getx() - ModelVelocity.getx() ) / AnalyticVelocity.getx() ) > 0.01 )	ReturnVal = -1;
	else if( (( AnalyticVelocity.gety() - ModelVelocity.gety() ) / AnalyticVelocity.gety() ) > 0.01 )	ReturnVal = -1;
	else if( (( AnalyticVelocity.getz() - ModelVelocity.getz() ) / AnalyticVelocity.getz() ) > 0.01 )	ReturnVal = -1;
	else										ReturnVal = 2;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModelVel = " << ModelVelocity << "m s^-1 : AnalyticVel = " << AnalyticVelocity << "m s^-1";
	std::cout << "\nPercentage Deviation: Mag = " << fabs(100-100*ModelVelocity.mag3()/AnalyticVelocity.mag3()) <<"%\n*****\n\n";
	std::cout << "\nPercentage Deviation: Xdir = " << fabs(100-100*ModelVelocity.getx()/AnalyticVelocity.getx()) <<"%\n*****\n\n";
	std::cout << "\nPercentage Deviation: Ydir = " << fabs(100-100*ModelVelocity.gety()/AnalyticVelocity.gety()) <<"%\n*****\n\n";
	std::cout << "\nPercentage Deviation: Zdir = " << fabs(100-100*ModelVelocity.getz()/AnalyticVelocity.getz()) <<"%\n*****\n\n";
	delete Sample;

	return ReturnVal;
}
