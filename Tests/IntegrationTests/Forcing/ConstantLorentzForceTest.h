#include "ForceModel.h"

int ConstantLorentzForceTest(char Element){
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
	double Size=5e-8; 		// m
	double Temp=280;		// K
	double Potential = 1;		// Normalised Potential
	Matter *Sample;			// Define the sample matter type

	// Set to true all Force models that are wanted
	bool Gravity = false;		// IS OFF!
	bool Centrifugal = false;	// IS OFF!
	bool Lorentz = true;		// IS ON!
	bool IonDrag = false;		// IS OFF!
	bool NeutralDrag = false;	// IS OFF!

	PlasmaData *Pdata = new PlasmaData();
	Pdata->ElectronTemp = 10*1.16e4;	// K, Electron Temperature, convert from eV
	threevector Efield(0.5, 2.0, 0.0);
	Pdata->ElectricField = Efield;
	threevector Bfield(0.0, 1.0, 0.0);
	Pdata->MagneticField = Bfield;

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
	threevector vinit(0.0,0.0,1.0); // Initial velocity with component parallel and perpendicular
	threevector zeros(0.0,0.0,0.0);
	Sample->update_motion(zeros,vinit,0.0);
	// START NUMERICAL MODEL
	ForceModel MyModel("out_ConstantForcingTest.txt",1.0,ForceModels,Sample,Pdata);

	double Mass = Sample->get_mass();
	MyModel.UpdateTimeStep();
	size_t imax(800000);

	for( size_t i(0); i < imax; i ++)
		MyModel.Force();
	double ModelTimeStep = MyModel.get_timestep();
	threevector ModelVelocity = Sample->get_velocity();
	// END NUMERICAL MODEL

	// START ANALYTICAL MODEL

	threevector Eparr(Bfield.getx()*Efield.getx(),Bfield.gety()*Efield.gety(),Bfield.getz()*Efield.getz());
	threevector Eperp = Efield - Eparr;

	threevector vE = (Eperp^Bfield)*(1/(pow(Bfield.mag3(),2)));

	threevector Vparr(Bfield.getx()*vinit.getx(),Bfield.gety()*vinit.gety(),Bfield.getz()*vinit.getz());
	threevector Vperp = vinit - Vparr - vE;
	double ConvertKelvsToeV(8.621738e-5);
	// NOTE: this is the charge to mass ratio for a negative grain only...
	double qtom = -3.0*epsilon0*Pdata->ElectronTemp*ConvertKelvsToeV*Sample->get_potential()
			/(pow(Sample->get_radius(),2)*Sample->get_density());
	double AngularVel = qtom*Bfield.mag3(); // Direction of this as vector is parallel to B

	double rc = 0.1; 
	double Vperpprime = rc*AngularVel;

	double vx = -Vperp.mag3()*sin(AngularVel*imax*MyModel.get_timestep());
	double vz = Vperp.mag3()*cos(AngularVel*imax*MyModel.get_timestep());

	std::cout << "\nVret = (" << vx << ", 0.0, " << vz << ")\n";
	std::cout << "\nVperpprime = (" << Vperpprime << ")\n";
	std::cout << "\nAngularVel = (" << AngularVel << ")\n";
	std::cout << "\nimax*ModelTimeStep = " << imax*ModelTimeStep;
	std::cout << "\nqtom = " << qtom;
	std::cout << "\nvE = " << vE;
	threevector AnalyticVelocity = Vparr + Eparr*qtom*imax*ModelTimeStep + vE + threevector(vx,0.0,vz); // + vE;
	// END ANALYTICAL MODEL

	double ReturnVal = 0;

	if( (ModelVelocity - AnalyticVelocity).mag3() == 0 ) 				ReturnVal = 1;
	else if( (( AnalyticVelocity.getx() - ModelVelocity.getx() ) / AnalyticVelocity.getx() ) > 0.01 )	ReturnVal = -1;
	else if( (( AnalyticVelocity.gety() - ModelVelocity.gety() ) / AnalyticVelocity.gety() ) > 0.01 )	ReturnVal = -1;
	else if( (( AnalyticVelocity.getz() - ModelVelocity.getz() ) / AnalyticVelocity.getz() ) > 0.01 )	ReturnVal = -1;
	else										ReturnVal = 2;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 2 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModelVel = " << ModelVelocity << "m s^-1 : AnalyticVel = " << AnalyticVelocity << "m s^-1";
	std::cout << "\nPercentage Deviation: Mag = " << fabs(100-100*ModelVelocity.mag3()/AnalyticVelocity.mag3()) <<"%\n*****\n\n";
	std::cout << "\nPercentage Deviation: Xdir = " << fabs(100-100*ModelVelocity.getx()/AnalyticVelocity.getx()) <<"%\n*****\n\n";
	std::cout << "\nPercentage Deviation: Ydir = " << fabs(100-100*ModelVelocity.gety()/AnalyticVelocity.gety()) <<"%\n*****\n\n";
	std::cout << "\nPercentage Deviation: Zdir = " << fabs(100-100*ModelVelocity.getz()/AnalyticVelocity.getz()) <<"%\n*****\n\n";
	delete Sample;

	return ReturnVal;
}
