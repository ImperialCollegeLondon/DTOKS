#include "HeatingModel.h"

int CompareConstEmissivTest(char Element,double Emissiv){
	clock_t begin = clock();
	// ********************************************************** //
	// FIRST, define program default behaviour

	// Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
	char EmissivityModel = 'f'; 	// Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
	char ExpansionModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
													// and (z)ero expansion
	char HeatCapacityModel = 'c'; 	// Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
	char BoilingModel = 'y'; 	// Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
													// and (t)homson
	char TimeStepType = 'f';	// Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
													// and (f)ixed steps
	std::string Name="constant";	// Describes heating model

 	// Parameters describing the heating model
	double Power=1e-9;		// Kilo-Watts power in addition to heating model powers
	double Size=1e-6; 		// m
	double Temp=280;		// K
	double TimeStep=1e-9;		// s
	Matter *Sample;			// Define the sample matter type

	// Set to true all heating models that are wanted
	bool RadiativeCooling = true;
	bool EvaporativeCooling = false;
	bool NewtonCooling = false;		// This model is equivalent to Electron and Ion heat flux terms
	// Plasma heating terms
	bool NeutralHeatFlux = false;
	bool ElectronHeatFlux = false;
	bool IonHeatFlux = false;
	bool NeutralRecomb = false;
	// Electron Emission terms
	bool TEE = false;
	bool SEE = false;

	PlasmaData *Pdata = new PlasmaData;
	Pdata->AmbientTemp = 0;

	bool PlasmaHeating = false; 		// If we want plasma heating terms turned off
	if( !PlasmaHeating ){
		NeutralHeatFlux = false;
		ElectronHeatFlux = false;
		IonHeatFlux = false;
		NeutralRecomb = false;
	}

	// Models and ConstModels are placed in an array in this order:
	std::array<bool, 9> Models = 
		{RadiativeCooling, EvaporativeCooling, NewtonCooling, IonHeatFlux, ElectronHeatFlux, NeutralHeatFlux, 
		NeutralRecomb, SEE, TEE };
	std::array<char, 4> ConstModels =
		{ EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};

	if	(Element == 'W'){ 
		Sample = new Tungsten(Size,Temp,ConstModels);
		TimeStep=1e-10;
	}else if (Element == 'B'){ 
		Sample = new Beryllium(Size,Temp,ConstModels);
		TimeStep=1e-10;
	}else if (Element == 'F'){
		Sample = new Iron(Size,Temp,ConstModels);
		TimeStep=1e-11;
	}else if (Element == 'G'){
		Sample = new Graphite(Size,Temp,ConstModels);
		TimeStep=1e-10;
	}else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}
	threevector xinit(1.15,0.0,-1.99);// default injection right hand side
	threevector vinit(0.0,0.0,0.0);
	Sample->update_motion(xinit,vinit,0.0);
	HeatingModel MyModel("out_ConstantHeatingTest.txt",1.0,Models,Sample,Pdata);
	MyModel.set_PowerIncident(Power);
	MyModel.UpdateTimeStep();
	MyModel.Vapourise();

	double ModelTime = MyModel.get_totaltime();
	

	double a = Power;
	std::cout << "Emiss = " << Sample->get_emissivity();
	double b = Emissiv*Sample->get_surfacearea()*Sigma;
	double ti(0), tf(0), t1(0), t2(0), t3(0), t4(0);
	double FinalTemp = pow(a/b,0.25)-0.000000001;
	if( FinalTemp < Sample->get_meltingtemp() ){ 
		tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
	}else{
		tf=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
		t2 = Sample->get_latentfusion()*Sample->get_mass()/(Power-b*Sample->get_meltingtemp()); 
		tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		std::cout << "\ntf is " << tf;
		t3 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
		if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
			t4 = Sample->get_latentvapour()*Sample->get_mass()/(Power-b*Sample->get_boilingtemp());
		}
	}

	if( Element == 'G' ){
		
		// Only one phase transition
                double p1 = atan(FinalTemp*pow(b/a,0.25));
                double p2 = atanh(FinalTemp*pow(b/a,0.25));
                double p3 = 2*pow(pow(a,3)*b,0.25);
                tf = (p1 + p2)/p3;


		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.24)))
					/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
		if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
                        t4 = Sample->get_latentvapour()*Sample->get_mass()/(Power-b*Sample->get_boilingtemp());
                }
		t2 = 0;
		t3 = 0;
	}

	if( t1 != t1 && t3 != t3 ){
		std::cout << "\nt1 AND t3 are nan!";
		return -1;
	}else if( t1 != t1 ){
		std::cout << "\nt1 is a nan!";
		return -1;
	}else if( t3 != t3 ){
		std::cout << "\nt3 is a nan!";
		return -1;
	}

	double AnalyticTime = t1 + t2 + t3 + t4;
	double ReturnVal = 0;

	if( ModelTime == AnalyticTime ) 			ReturnVal = 1;
	else if( fabs(1-ModelTime/AnalyticTime) < 0.01 ) 	ReturnVal = 2;
	else							ReturnVal = -1;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModel Time = " << ModelTime << "s : Analytic Time = " << AnalyticTime << "s";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelTime/AnalyticTime) <<"%\n*****\n\n";
	std::cout << "\n\n*****\nModel Temp = " << Sample->get_temperature() << "K : Analytic Temp = " << FinalTemp << "K";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*Sample->get_temperature()/FinalTemp) <<"%\n*****\n\n";


	return ReturnVal;
}
