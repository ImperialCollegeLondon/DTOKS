#include "HeatingModel.h"

int ConstantPlasmaHeatingTest(char Element){
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
	char TimeStepType = 'f';	// Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
													// and (f)ixed steps
	std::string Name="constant";	// Describes heating model

 	// Parameters describing the heating model
	double Power=1e-8;		// Kilo-Watts power in addition to heating model powers
	double Size=5e-8; 		// m
	double Temp=280;		// K
	double TimeStep=1e-12;		// s
	Matter *Sample;			// Define the sample matter type

	// Set to true all heating models that are wanted
	bool RadiativeCooling = true;
	bool EvaporativeCooling = false;
	bool NewtonCooling = false;		// This model is equivalent to Electron and Ion heat flux terms
	// Plasma heating terms
	bool NeutralHeatFlux = true;
	bool ElectronHeatFlux = true;
	bool IonHeatFlux = true;
	bool NeutralRecomb = false;
	// Electron Emission terms
	bool TEE = false;
	bool SEE = false;

	PlasmaData Pdata;
	Pdata.NeutralDensity = 3e19;		// m^-3, Neutral density
	Pdata.ElectronDensity = 8e17;	 	// m^-3, Electron density
	double Potential = 1;			// arb, assumed negative, potential normalised to dust temperature, (-e*phi)/(Kb*Td)
	Pdata.IonTemp = 100*1.16e5;	 	// K, Ion Temperature
	Pdata.ElectronTemp = 100*1.16e5;	// K, Electron Temperature, convert from eV
	Pdata.NeutralTemp = 100*1.16e5; 	// K, Neutral Temperature, convert from eV
	Pdata.AmbientTemp = 0;

	// Models and ConstModels are placed in an array in this order:
	std::array<bool, 9> Models = 
		{RadiativeCooling, EvaporativeCooling, NewtonCooling, IonHeatFlux, ElectronHeatFlux, NeutralHeatFlux, 
		NeutralRecomb, SEE, TEE };
	std::array<char, 4> ConstModels =
		{ EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel};

	if	(Element == 'W'){ 
		Sample = new Tungsten(Size,Temp,ConstModels);
		TimeStep=1e-12;
	}else if (Element == 'B'){ 
		Sample = new Beryllium(Size,Temp,ConstModels);
		TimeStep=1e-11;
	}else if (Element == 'F'){
		Sample = new Iron(Size,Temp,ConstModels);
		TimeStep=1e-13;
	}else if (Element == 'G'){
		Sample = new Graphite(Size,Temp,ConstModels);
		TimeStep=1e-11;
	}else{ 
		std::cerr << "\nInvalid Option entered";
		return -1;
	}

	HeatingModel MyModel("out_ConstantHeatingTest.txt",1.0,Models,Sample,Pdata);
//	MyModel.Vapourise(("out_ConstantHeatingTest_" + Element + ".txt").c_str(),ConstModels,TimeStepType);
//	MyModel.Vapourise("out_ConstantHeatingTest.txt",ConstModels,TimeStepType);
	MyModel.Vapourise();

	double ModelTime = MyModel.get_totaltime();
	

	// *********************** BEGIN ANALYTICAL MODEL ************************************ //

	double ElectronFlux = Pdata.ElectronDensity*exp(-Potential)*sqrt(Kb*Pdata.ElectronTemp/(2*PI*Me));
	double NeutralFlux = Pdata.NeutralDensity*sqrt(Kb*Pdata.NeutralTemp/(2*PI*Mp));
	double IonFlux = ElectronFlux;

	double ElectronFluxPower = Sample->get_surfacearea()*2*ElectronFlux*Pdata.ElectronTemp*Kb/1000; // Convert from Joules to KJ
	double NeutralFluxPower = Sample->get_surfacearea()*2*NeutralFlux*Pdata.NeutralTemp*Kb/1000; // Convert from Joules to KJ
	double IonFluxPower = (Sample->get_surfacearea()*IonFlux*Pdata.IonTemp*Kb/1000) // Convert from Joules to KJ
	*(2+2*Potential*(Pdata.ElectronTemp/Pdata.IonTemp)+pow(Potential*(Pdata.ElectronTemp/Pdata.IonTemp),2))/(1+Potential*(Pdata.ElectronTemp/Pdata.IonTemp));

	double a = Power+ElectronFluxPower+NeutralFluxPower+IonFluxPower;

//	std::cout << "Emiss = " << Sample->get_emissivity();
	double b = Sample->get_emissivity()*Sample->get_surfacearea()*Sigma/1000;
	double ti(0), tf(0), t1(0), t2(0), t3(0), t4(0);
	double FinalTemp = pow(a/b,0.25)-0.000000001;
//	std::cout << "\nFinalTemp = " << FinalTemp; std::cin.get();
	if( FinalTemp < Sample->get_meltingtemp() ){ 
		tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
	}else{
		if( FinalTemp > Sample->get_boilingtemp() )
			FinalTemp = Sample->get_boilingtemp();
		tf=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
		t2 = Sample->get_latentfusion()*Sample->get_mass()/(a-b*Sample->get_meltingtemp()); 
		tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
//		std::cout << "\ntf is " << tf;
		t3 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
		if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
			t4 = Sample->get_latentvapour()*Sample->get_mass()/(a-b*Sample->get_boilingtemp());
		}

	}

	if( Element == 'G' ){
		if( FinalTemp > Sample->get_boilingtemp() )
		FinalTemp = Sample->get_boilingtemp();
	
		// Only one phase transition
//double p1 = atan(FinalTemp*pow(b/a,0.25));
//double p2 = atanh(FinalTemp*pow(b/a,0.25));
//double p3 = 2*pow(pow(a,3)*b,0.25);
                tf = (atan(FinalTemp*pow(b/a,0.25)) + atanh(FinalTemp*pow(b/a,0.25)))/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);

		if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
                        t4 = Sample->get_latentvapour()*Sample->get_mass()/(a-b*Sample->get_boilingtemp());
                }
		t2 = 0;
		t3 = 0;
		std::cout << "\nt1 = " << t1 << "\nt2 = " << t2 << "\nt3 = " << t3 << "\nt4 = " << t4; 

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
