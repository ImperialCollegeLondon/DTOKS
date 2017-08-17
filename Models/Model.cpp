//#define PAUSE
//#define MODEL_DEBUG

#include "Model.h"

struct PlasmaData PlasmaDefaults = {
	1e20,		// m^-3, Neutral Density
	1e20,		// m^-3, Electron Density
	1e20,		// m^-3, Ion Density
	116045.25,	// K, Ion Temperature
	116045.25,	// K, Electron Temperature
	116045.25,	// K, Neutral Temperature
	300,		// K, Ambient Temperature
	threevector(),	// m s^-1, Plasma Velocity (Should eventually be normalised to sound speed cs)
	threevector(),	// m s^-2, Acceleration due to gravity
	threevector(),	// V m^-1, Electric field at dust location (Normalised later) 
	threevector(),	// T, Magnetic field at dust location (Normalised later)
};
//PlasmaGrid *DefaultGrid = new PlasmaGrid('h','m',0.01);

Model::Model():Sample(new Tungsten),Pgrid(new PlasmaGrid('h','m',0.01)),Pdata(&PlasmaDefaults),Accuracy(1.0),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model():Sample(new Tungsten),Pgrid('h','m',0.01)Pdata(PlasmaDefaults),Accuracy(1.0),ContinuousPlasma(true)\n\n");
	i = 0; k = 0;
	PlasmaDataFile.open("pd.txt");
	update_plasmadata();
}

// Constructor for Matter sample sitting in a constant plasma background given by PlasmaData (pdata) with a Default grid
Model::Model( Matter *&sample, PlasmaData *&pdata, double accuracy )
		:Sample(sample),Pgrid(new PlasmaGrid('h','m',0.01)),Pdata(pdata),Accuracy(accuracy),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model( Matter *& sample, PlasmaData *&pdata, double accuracy ):Sample(sample),Pgrid(new PlasmaGrid('h','m',0.01)),Pdata(pdata),Accuracy(accuracy),ContinuousPlasma(true)\n\n");
	assert(Accuracy > 0);
	i = 0; k = 0;
	PlasmaDataFile.open("pd.txt");
	update_plasmadata(pdata);
//	std::cout << "\nAccuracy = " << Accuracy;
}

// Constructor for Matter sample moving in a varying plasma background given by PlasmaGrid (pgrid) with the current information 
// stored in pdata.
Model::Model( Matter *&sample, PlasmaGrid &pgrid, double accuracy )
		:Sample(sample),Pgrid(&pgrid),Pdata(&PlasmaDefaults),Accuracy(accuracy),ContinuousPlasma(false),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model( Matter *& sample, PlasmaGrid &pgrid, double accuracy ):Sample(sample),Pgrid(pgrid),Pdata(PlasmaDefaults),Accuracy(accuracy), ContinuousPlasma(false)\n\n");
	assert(Accuracy > 0);
	i = 0; k = 0;
	PlasmaDataFile.open("pd.txt");
	update_plasmadata();
	//	std::cout << "\nAccuracy = " << Accuracy;
}

bool Model::new_cell()const{
	int j(0), p(0);
	Pgrid->locate(j,p,Sample->get_position());
	if( j != i || p != k )	return true;
	return false;
}

void Model::update_plasmadata(PlasmaData *&pdata){
	Mo_Debug( "\tIn Model::update_plasmadata(PlasmaData *&pdata)\n\n");
	Pdata = pdata;
}

bool Model::update_plasmadata(){
	Mo_Debug( "\tIn Model::update_plasmadata()\n\n");
	double ConvertJtoK(7.24297166e22);		// Conversion factor from ev to K
	bool InGrid = Pgrid->locate(i,k,Sample->get_position());
	if( !InGrid ) return InGrid;			// Particle has escaped simulation domain
	update_fields(i,k);
	Pdata->NeutralDensity 	= Pgrid->getna0(i,k);  	// NEUTRAL DENSITY EQUALS ION DENSITY
	Pdata->ElectronDensity 	= Pgrid->getna1(i,k);  
	Pdata->IonDensity 	= Pgrid->getna0(i,k);
	Pdata->IonTemp		= Pgrid->getTi(i,k)*ConvertJtoK;
	Pdata->ElectronTemp 	= Pgrid->getTe(i,k)*ConvertJtoK;
	Pdata->NeutralTemp 	= Pgrid->getTi(i,k)*ConvertJtoK; 	// NEUTRAL TEMP EQUAL TO ION TEMP
//	std::cout << "\nTi = " << Pdata->IonTemp;
//	std::cout << "\nTe = " << Pdata->ElectronTemp;
//	std::cout << "\nTn = " << Pdata->NeutralTemp; std::cin.get();
	Pdata->AmbientTemp 	= 300; 			// NOTE THIS IS HARD CODED OHMEINGOD

	return true;
//	std::cout << "\t\tPdata->MagneticField = " << Pdata->MagneticField << "\n";
}

// CHECK THIS FUNCTION IS THE SAME AS IT WAS BEFORE!
void Model::update_fields(int i, int k){
	Mo_Debug( "\tIn Model::update_fields(int " << i << ", int " << k << ")\n\n");
	threevector vp, E, B;

	// Get Average plasma velocity
	double aveu(0.0);
	if( (Pgrid->getna0(i,k)>0.0) || (Pgrid->getna1(i,k)>0.0) ){
		aveu = (Pgrid->getna0(i,k)*Pgrid->getua0(i,k)+Pgrid->getna1(i,k)*Pgrid->getua1(i,k))
			/(Pgrid->getna0(i,k)+Pgrid->getna1(i,k));
	}
	else aveu = 0.0;

	// Read magnetic field in
	B.setx(Pgrid->getbx(i,k));
	B.sety(Pgrid->getby(i,k));
	B.setz(Pgrid->getbz(i,k));

	// Plasma velocity is parallel to the B field
	vp.setx(aveu*(B.getunit().getx()));
	vp.sety(aveu*(B.getunit().gety()));
	vp.setz(aveu*(B.getunit().getz()));

	if(Pgrid->getgridflag(i,k)==1){
		if( (Pgrid->getgridflag(i+1,k)==1) && (Pgrid->getgridflag(i-1,k)==1) ){
			E.setx(-(Pgrid->getpo(i+1,k)-Pgrid->getpo(i-1,k))/(2.0*Pgrid->getdl()));
		}else if(Pgrid->getgridflag(i+1,k)==1){
			E.setx(-(Pgrid->getpo(i+1,k)-Pgrid->getpo(i,k))/Pgrid->getdl());
		}else if(Pgrid->getgridflag(i-1,k)==1){
			E.setx(-(Pgrid->getpo(i,k)-Pgrid->getpo(i-1,k))/Pgrid->getdl());
		}else E.setx(0.0);
		if((Pgrid->getgridflag(i,k+1)==1)&&(Pgrid->getgridflag(i,k-1)==1)){
			E.setz(-(Pgrid->getpo(i,k+1)-Pgrid->getpo(i,k-1))/(2.0*Pgrid->getdl()));
		}else if(Pgrid->getgridflag(i,k+1)==1){
			E.setz(-(Pgrid->getpo(i,k+1)-Pgrid->getpo(i,k))/Pgrid->getdl());
		}else if(Pgrid->getgridflag(i,k-1)==1){
			E.setz(-(Pgrid->getpo(i,k)-Pgrid->getpo(i,k-1))/Pgrid->getdl());
		}else E.sety(0.0);
	}else{
		E.setx(0.0);
		E.setz(0.0);
	}

        Pdata->PlasmaVel         = vp;
	Pdata->Gravity 		= threevector(0.0,0.0,-9.8);
        Pdata->ElectricField     = E;
        Pdata->MagneticField     = B;
}

void Model::RecordPlasmadata(){
	PlasmaDataFile << "\n" << Pdata->NeutralDensity << "\t" << Pdata->ElectronDensity << "\t" << Pdata->IonDensity
			<< "\t" << Pdata->IonTemp << "\t" << Pdata->ElectronTemp << "\t" << Pdata->NeutralTemp
			<< "\t" << Pdata->AmbientTemp << "\t" << Pdata->PlasmaVel << "\t" << Pdata->Gravity 
			<< "\t" << Pdata->ElectricField << "\t" << Pdata->MagneticField;
}

const double Model::ElectronFlux(double DustTemperature)const{
	return DTOKSElectronFlux(DustTemperature);
//	return OMLElectronFlux(DustTemperature);
}

const double Model::IonFlux(double DustTemperature)const{
	return DTOKSIonFlux(DustTemperature);
//	return OMLIonFlux(DustTemperature);
}

const double Model::DTOKSIonFlux(double DustTemperature)const{
	H_Debug("\n\tIn Model::IonFlux():");

	double IonFlux=0;

	if( Sample->is_positive() ) IonFlux = DTOKSElectronFlux(DustTemperature); //Positive grain, DeltaTot() > 1
	else	IonFlux = DTOKSElectronFlux(DustTemperature)*(1-Sample->get_deltatot());
	assert(IonFlux >= 0);
	return IonFlux;
}

const double Model::DTOKSElectronFlux(double DustTemperature)const{
	H_Debug("\n\tIn Model::ElectronFlux():\n\n");

	return Pdata->ElectronDensity*exp(-Sample->get_potential())*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
}

const double Model::OMLIonFlux(double DustTemperature)const{
	H_Debug("\n\tIn Model::IonFlux():");

	double IonFlux=0;

	if( Pdata->IonTemp <= 0 ) return 0.0;
	if( Sample->is_positive() ) IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Mp))*exp(Sample->get_potential());
	else	IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Mp))*(1+Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp));
	assert(IonFlux >= 0);
	return IonFlux;
}

const double Model::OMLElectronFlux(double DustTemperature)const{
	H_Debug("\n\tIn Model::ElectronFlux():\n\n");
	if( Pdata->IonTemp <= 0 ) return 0.0;
	if(Sample->is_positive()) return Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me))*(1-Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp));

	else 	return Pdata->ElectronDensity*exp(-Sample->get_potential())*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
}

const double Model::NeutralFlux()const{
	H_Debug("\n\tIn Model::NeutralFlux():\n\n");
	return Pdata->NeutralDensity*sqrt(Kb*Pdata->NeutralTemp/(2*PI*Mp));
}


