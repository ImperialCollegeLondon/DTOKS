#include "HeatingModel.h"

int ConstantHeatEmissivTest(char Element){
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
	char BreakupModel = 'n';	// Possible values 'r', 'e', 'b'  and 'n': Corresponding to (r)otational, (e)lectrostatic, (b)oth and (n)o
	char TimeStepType = 'f';	// Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
													// and (f)ixed steps
	std::string Name="constant";	// Describes heating model

 	// Parameters describing the heating model
	double Power=1e-4;		// Kilo-Watts power in addition to heating model powers
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
	std::array<bool, HMN> Models = 
		{RadiativeCooling, EvaporativeCooling, NewtonCooling, IonHeatFlux, ElectronHeatFlux, NeutralHeatFlux, 
		NeutralRecomb, SEE, TEE };
	std::array<char, CM> ConstModels =
		{ EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel,BreakupModel};

	if	(Element == 'W'){ 
		Sample = new Tungsten(Size,Temp,ConstModels);
		TimeStep=1e-12;
	}else if (Element == 'B'){ 
		Sample = new Beryllium(Size,Temp,ConstModels);
		TimeStep=1e-12;
	}else if (Element == 'F'){
		Sample = new Iron(Size,Temp,ConstModels);
		TimeStep=1e-14;
	}else if (Element == 'G'){
		Sample = new Graphite(Size,Temp,ConstModels);
		TimeStep=1e-12;
	}else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}
	threevector xinit(1.15,0.0,-1.99);// default injection right hand side
	threevector vinit(0.0,0.0,0.0);
	Sample->update_motion(xinit,vinit,0.0);
	HeatingModel MyModel("Tests/IntegrationTests/Data/out_ConstantHeatingTest.txt",1.0,Models,Sample,Pdata);
//	MyModel.Vapourise("out_ConstantHeatingTest.txt",ConstModels,TimeStepType);
	
	double Mass = Sample->get_mass();
	double SA = Sample->get_surfacearea();
	MyModel.set_PowerIncident(Power);
	MyModel.UpdateTimeStep();
	MyModel.Vapourise();
	double ModelTime = MyModel.get_totaltime();
	
	double a = Power;
	double b = Sample->get_emissivity()*SA*Sigma/1000;

	double FinalTemp = pow(a/b,1.0/4.0)-0.000000001;

	double ti(0), tf(0), t1(0), t2(0), t3(0), t4(0);
	if( FinalTemp >= Sample->get_boilingtemp() ){
		FinalTemp = Sample->get_boilingtemp();
		t4 = Sample->get_latentvapour()*Mass/(a-b*pow(Sample->get_boilingtemp(),4));
	}
	if( Element != 'G' ){
		if( FinalTemp < Sample->get_meltingtemp() ){
			tf=(atan(FinalTemp*pow(b/a,1.0/4.0))+atanh(FinalTemp*pow(b/a,1.0/4.0)))
						/(2*pow(pow(a,3)*b,1.0/4.0));
			ti=(atan(Temp*pow(b/a,1.0/4.0))+atanh(Temp*pow(b/a,1.0/4.0)))
						/(2*pow(pow(a,3)*b,1.0/4.0));
			t1 = Mass*Sample->get_heatcapacity()*(tf-ti);
		}else{
			tf=(atan(Sample->get_meltingtemp()*pow(b/a,1.0/4.0))+atanh(Sample->get_meltingtemp()*pow(b/a,1.0/4.0)))
						/(2*pow(pow(a,3)*b,1.0/4.0));
			ti=(atan(Temp*pow(b/a,1.0/4.0))+atanh(Temp*pow(b/a,1.0/4.0)))
						/(2*pow(pow(a,3)*b,1.0/4.0));
			t1 = Mass*Sample->get_heatcapacity()*(tf-ti);
	
			t2 = Sample->get_latentfusion()*Mass/(a-b*Sample->get_meltingtemp());

			tf=(atan(FinalTemp*pow(b/a,1.0/4.0))+atanh(FinalTemp*pow(b/a,1.0/4.0)))
						/(2*pow(pow(a,3)*b,1.0/4.0));
			ti=(atan(Sample->get_meltingtemp()*pow(b/a,1.0/4.0))+atanh(Sample->get_meltingtemp()*pow(b/a,1.0/4.0)))
						/(2*pow(pow(a,3)*b,1.0/4.0));	
	
			t3 = Mass*Sample->get_heatcapacity()*(tf-ti);
		}
	}else{
		tf=(atan(FinalTemp*pow(b/a,1.0/4.0))+atanh(FinalTemp*pow(b/a,1.0/4.0)))
					/(2*pow(pow(a,3)*b,1.0/4.0));
		ti=(atan(Temp*pow(b/a,1.0/4.0))+atanh(Temp*pow(b/a,1.0/4.0)))
					/(2*pow(pow(a,3)*b,1.0/4.0));
		t1 = Mass*Sample->get_heatcapacity()*(tf-ti);
	}

	if( t1 != t1 || t3 != t3 ){
		std::cout << "\nThermal Equilibrium reached before end of test!";
		return -1;
	}

//	std::cout << "\nt1 = " << t1 << "\nt2 = " << t2 << "\nt3 = " << t3 << "\nt4 = " << t4;
	
	double AnalyticTime = t1 + t2 + t3 + t4;	

	double ReturnVal = 0;

	if( ModelTime == AnalyticTime ) 			ReturnVal = 1;
	else if( fabs(1-ModelTime/AnalyticTime) < 0.01 ) 	ReturnVal = 2;
	else							ReturnVal = -1;

	clock_t end = clock();
	double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
		
	std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
	std::cout << "\n\n*****\nModelTime = " << ModelTime << "s : AnalyticTime = " << AnalyticTime << "s";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelTime/AnalyticTime) <<"%\n*****\n\n";
	std::cout << "\n\n*****\nModel Temp = " << Sample->get_temperature() << "K : Analytic Temp = " << FinalTemp << "K";
	std::cout << "\nPercentage Deviation = " << fabs(100-100*Sample->get_temperature()/FinalTemp) <<"%\n*****\n\n";


	return ReturnVal;
}

