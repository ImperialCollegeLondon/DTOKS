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
	char TimeStepType = 'f';	// Possible values 'o', 's' and 'f': Corresponding to (o)ne degree steps, (s)mall steps 
													// and (f)ixed steps
	std::string Name="constant";	// Describes heating model

 	// Parameters describing the heating model
	double Power=1e-6;		// Kilo-Watts power in addition to heating model powers
	double Size=5e-8; 		// m
	double Temp=270;		// K
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
	bool NeutralRecomb = true;
	// Electron Emission terms
	bool TEE = false;
	bool SEE = false;

	PlasmaData Pdata;
	Pdata.AmbientTemp = 0;

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

	HeatingModel MyModel("out_ConstantHeatingTest.txt",1.0,Models,Sample,Pdata);
//	MyModel.Vapourise("out_ConstantHeatingTest.txt",ConstModels,TimeStepType);
	MyModel.Vapourise();
	double ModelTime = MyModel.get_totaltime();
	

	double a = Power;
	double b = Sample->get_emissivity()*Sample->get_surfacearea()*Sigma;
	double ti(0), tf(0), t1(0), t3(0);
	if( Element != 'G' ){
		tf=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
					/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);

	

		double p1 = atan(Sample->get_boilingtemp()*pow(b/a,0.25));
		double p2 = atanh(Sample->get_boilingtemp()*pow(b/a,0.25));
		double p3 = 2*pow(pow(a,3)*b,0.25);

		tf = (p1 + p2)/p3;
		p1 = atan(Sample->get_meltingtemp()*pow(b/a,0.25));
		p2 = atanh(Sample->get_meltingtemp()*pow(b/a,0.25));
		p3 = 2*pow(pow(a,3)*b,0.25);
		ti= (p1 + p2)/p3;

		t3 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);

	}else{
		
		// Only one phase transition
                double p1 = atan(Sample->get_boilingtemp()*pow(b/a,0.25));
                double p2 = atanh(Sample->get_boilingtemp()*pow(b/a,0.25));
                double p3 = 2*pow(pow(a,3)*b,0.25);
                tf = (p1 + p2)/p3;


		ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.24)))
					/(2*pow(pow(a,3)*b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);


	}

	if( t1 != t1 || t3 != t3 ){
		std::cout << "\nThermal Equilibrium reached before end of test!";
		return -1;
	}

	double t2 = Sample->get_latentfusion()*Sample->get_mass()/(Power-b*Sample->get_meltingtemp()); 
	double t4 = Sample->get_latentvapour()*Sample->get_mass()/(Power-b*Sample->get_boilingtemp());
	
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

	return ReturnVal;
}


/*
//		std::cout << "\n\np1 = " << atan(Sample->get_meltingtemp()*pow(b/a,0.25)) << "\np2 = " << atanh(Sample->get_meltingtemp()*pow(b/a,0.25)) << "\np3 = " << (2*pow(pow(a,3)*b,0.25));
//		std::cout << "\n\nSample->get_boilingtemp() = " << Sample->get_boilingtemp();

//		std::cout << "\nARG = " << Sample->get_boilingtemp()*pow(b/a,0.25);
//		std::cout << "\n\np1 = " << p1 << "\np2 = " << p2 << "\np3 = " << p3;
//		std::cout << "\n\nSample->get_boilingtemp() = " << Sample->get_boilingtemp();


//		std::cout << "\n\ntf = " << tf;
//		std::cout << "\n\nti = " << ti;
	
//		std::cout << "\n\np1 = " << p1 << "\np2 = " << p2 << "\np3 = " << p3;
//		std::cout << "\n\nSample->get_boilingtemp() = " << Sample->get_meltingtemp();



		tf = (log10((pow(a,0.25)+pow(b,0.25)*Sample->get_meltingtemp())/(pow(a,0.25)-pow(b,0.25)*Sample->get_meltingtemp()))+2*atan(pow(b/a,0.25)*Sample->get_meltingtemp()))/(4*pow(a,0.75)*pow(b,0.25));
		
		ti = (log10((pow(a,0.25)+pow(b,0.25)*Temp)/(pow(a,0.25)-pow(b,0.25)*Temp))+2*atan(pow(b/a,0.25)*Temp))/(4*pow(a,0.75)*pow(b,0.25));
		t1 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
		tf = (log10((pow(a,0.25)+pow(b,0.25)*Sample->get_boilingtemp()))+2*atan(pow(b/a,0.25)*Sample->get_boilingtemp()))/(4*pow(a,0.75)*pow(b,0.25));
		ti = (log10((pow(a,0.25)+pow(b,0.25)*Sample->get_meltingtemp()))-log10((pow(a,0.25)-pow(b,0.25)*Sample->get_meltingtemp()))+2*atan(pow(b/a,0.25)*Sample->get_meltingtemp()))/(4*pow(a,0.75)*pow(b,0.25));
		t3 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
*/


/*
		ti = (log10((pow(a,0.25)+pow(b,0.25)*Temp)/(pow(a,0.25)-pow(b,0.25)*Temp))+2*atan(pow(b/a,0.25)*Temp))/(4*pow(a,0.75)*pow(b,0.25));
		tf = (log10((pow(a,0.25)+pow(b,0.25)*Sample->get_boilingtemp()))+2*atan(pow(b/a,0.25)*Sample->get_boilingtemp()))/(4*pow(a,0.75)*pow(b,0.25));
		t3 = Sample->get_mass()*Sample->get_heatcapacity()*(tf-ti);
*/

