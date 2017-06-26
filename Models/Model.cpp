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
PlasmaGrid *DefaultGrid = new PlasmaGrid('h','m',0.01);


Model::Model():Sample(new Tungsten),Pgrid(DefaultGrid),Pdata(&PlasmaDefaults),Accuracy(1.0),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model():Sample(new Tungsten),Pgrid('h','m',0.01)Pdata(PlasmaDefaults),Accuracy(1.0),ContinuousPlasma(true)\n\n");
	update_plasmadata(Sample->get_position());
}

// Constructor for Matter sample sitting in a constant plasma background given by PlasmaData (pdata) with a Default grid
Model::Model( Matter *&sample, PlasmaData *&pdata, double accuracy )
		:Sample(sample),Pgrid(DefaultGrid),Pdata(pdata),Accuracy(accuracy),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model( Matter *& sample, PlasmaData *&pdata, double accuracy ):Sample(sample),Pgrid(DefaultGrid),Pdata(pdata),Accuracy(accuracy),ContinuousPlasma(true)\n\n");
	assert(Accuracy > 0);
	update_plasmadata(pdata);
//	std::cout << "\nAccuracy = " << Accuracy;
}

// Constructor for Matter sample moving in a varying plasma background given by PlasmaGrid (pgrid) with the current information 
// stored in pdata.
Model::Model( Matter *&sample, PlasmaGrid &pgrid, double accuracy )
		:Sample(sample),Pgrid(&pgrid),Pdata(&PlasmaDefaults),Accuracy(accuracy),ContinuousPlasma(false),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model( Matter *& sample, PlasmaGrid &pgrid, double accuracy ):Sample(sample),Pgrid(pgrid),Pdata(PlasmaDefaults),Accuracy(accuracy), ContinuousPlasma(false)\n\n");
	assert(Accuracy > 0);
	update_plasmadata(sample->get_position());
	//	std::cout << "\nAccuracy = " << Accuracy;
}

void Model::update_plasmadata(PlasmaData *&pdata){
	Mo_Debug( "\tIn Model::update_plasmadata(PlasmaData *&pdata->\n\n");
	Pdata->NeutralDensity 	= pdata->NeutralDensity;
	Pdata->ElectronDensity 	= pdata->ElectronDensity;
	Pdata->IonDensity 	= pdata->IonDensity;
	Pdata->IonTemp		= pdata->IonTemp;
	Pdata->ElectronTemp 	= pdata->ElectronTemp;
	Pdata->NeutralTemp 	= pdata->NeutralTemp;
	Pdata->AmbientTemp 	= pdata->AmbientTemp;
	Pdata->PlasmaVel         = pdata->PlasmaVel;
	Pdata->Gravity 		= pdata->Gravity;
	Pdata->ElectricField     = pdata->ElectricField;
	Pdata->MagneticField     = pdata->MagneticField;
}

bool Model::update_plasmadata(threevector pos){
	Mo_Debug( "\tIn Model::update_plasmadata(" << pos << ")\n\n");
	double ConvertevtoK(1.16e4);		// Conversion factor from ev to K
	int i(0), k(0);
	bool InGrid = Pgrid->locate(i,k,pos);
	if( !InGrid ) return InGrid;			// Particle has escaped simulation domain
	update_fields(i,k);
	Pdata->NeutralDensity 	= Pgrid->getna0(i,k);  	// NEUTRAL DENSITY EQUALS ION DENSITY
	Pdata->ElectronDensity 	= Pgrid->getna1(i,k);  
	Pdata->IonDensity 	= Pgrid->getna0(i,k);
	Pdata->IonTemp		= Pgrid->getTi(i,k)*ConvertevtoK;
	Pdata->ElectronTemp 	= Pgrid->getTe(i,k)*ConvertevtoK;
	Pdata->NeutralTemp 	= Pgrid->getTi(i,k)*ConvertevtoK; 	// NEUTRAL TEMP EQUAL TO ION TEMP
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
	double aveu;
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

