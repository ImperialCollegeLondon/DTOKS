//#define PAUSE
//#define MODEL_DEBUG

#include "Model.h"

//PlasmaGrid *DefaultGrid = new PlasmaGrid('h','m',0.01);
struct PlasmaData PlasmaDataDefaults = {
	1e20,		// m^-3, Neutral Density
	1e20,		// m^-3, Electron Density
	1e20,		// m^-3, Ion Density
	116045.25,	// K, Ion Temperature
	116045.25,	// K, Electron Temperature
	116045.25,	// K, Neutral Temperature
	300,		// K, Ambient Temperature
	1.66054e-27,// kg, Mass of ions
	threevector(),	// m s^-1, Plasma Velocity (Should eventually be normalised to sound speed cs)
	threevector(),	// m s^-2, Acceleration due to gravity
	threevector(),	// V m^-1, Electric field at dust location (Normalised later) 
	threevector(),	// T, Magnetic field at dust location (Normalised later)
};

struct PlasmaGrid_Data PlasmaGrid_DataDefaults = {
	// Plasma Parameters
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<double> >(),
	std::vector< std::vector<int> >	(),

	251,
	401,
	0,
	1.5,
	4.0,
	-2.0,
	2.0,

	1.66054e-27,
	'h',
	'j',
};

Model::Model():Sample(new Tungsten),PG_data(&PlasmaGrid_DataDefaults),Pdata(&PlasmaDataDefaults),Accuracy(1.0),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model():Sample(new Tungsten),PG_data(PlasmaGridDefaults),Pdata(&PlasmaDefaults),Accuracy(1.0),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0)\n\n");
	i = 0; k = 0;
	PlasmaDataFile.open("Data/pd.txt");
	PlasmaDataFile << "#t\ti\tk\tNn\tNe\tNi\tTi\tTe\tTn\tT0\tPvel\tgravity\tE\tB";
	PlasmaDataFile.close();
	PlasmaDataFile.clear();
	update_plasmadata();
}

// Constructor for Matter sample sitting in a constant plasma background given by PlasmaData (pdata) with a Default grid
Model::Model( Matter *&sample, PlasmaData &pdata, float accuracy )
		:Sample(sample),PG_data(&PlasmaGrid_DataDefaults),Pdata(std::make_shared<PlasmaData>(pdata)),Accuracy(accuracy),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model( Matter *&sample, PlasmaData *&pdata, float accuracy ):Sample(sample),PG_data(PlasmaGridDefaults),Pdata(pdata),Accuracy(accuracy),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0)\n\n");
	assert(Accuracy > 0);
	i = 0; k = 0;
	PlasmaDataFile.open("Data/pd.txt");
	PlasmaDataFile << "#t\ti\tk\tNn\tNe\tNi\tTi\tTe\tTn\tT0\tPvel\tgravity\tE\tB";
	PlasmaDataFile.close();
	PlasmaDataFile.clear();
	update_plasmadata(pdata);
}

// Constructor for Matter sample moving in a varying plasma background given by PlasmaGrid (pgrid) with the current information 
// stored in pdata.
Model::Model( Matter *&sample, PlasmaGrid_Data &pgrid, float accuracy )
		:Sample(sample),PG_data(std::make_shared<PlasmaGrid_Data>(pgrid)),Pdata(&PlasmaDataDefaults),Accuracy(accuracy),ContinuousPlasma(false),TimeStep(0.0),TotalTime(0.0){
	Mo_Debug("\n\nIn Model::Model( Matter *& sample, PlasmaGrid_Data &pgrid, float accuracy ):Sample(sample),PG_data(pgrid),Pdata(PlasmaDefaults),Accuracy(accuracy), ContinuousPlasma(false)\n\n");
	assert(Accuracy > 0);
	i = 0; k = 0;
	PlasmaDataFile.open("Data/pd.txt");
	PlasmaDataFile << "#t\ti\tk\tNn\tNe\tNi\tTi\tTe\tTn\tT0\tPvel\tgravity\tE\tB\n\n";
	PlasmaDataFile.close();
	PlasmaDataFile.clear();
	update_plasmadata();
}

// *************************************** MUTATION FUNCTIONS *************************************** //

// Locate dust particle in the plasma grid
const bool Model::locate(int &i, int &k, const threevector xd)const{
	P_Debug("\tIn Model::locate(int &" << i << ", int &" << k << ", " << xd << ")\n\n");
	// Adding 0.5 makes the rounding work properly
	i = int(0.5+(xd.getx()-get_gridxmin())/get_dlx());
	k = int(0.5+(xd.getz()-get_gridzmin())/get_dlz());
	if( xd.getx() < get_gridxmin() )
		i = -1;
	if( xd.getz() < get_gridzmin() )
		k = -1;
	return checkingrid(i,k);

}

const bool Model::checkingrid(const int i, const int k)const{
	P_Debug("\tIn Model::checkingrid(int " << i << ", int " << k  << ")\n\n");
	bool returnval(true);
	if( i >= get_gridx() || i < 0) returnval = false;
	if( k >= get_gridz() || k < 0) returnval = false;
	return returnval;
}


bool Model::new_cell()const{
	int j(0), p(0);
	locate(j,p,Sample->get_position());
	if( (j == i) && (p == k) )	return false;
	return true;
}

void Model::close_file(){
	ModelDataFile.close();
}

void Model::update_plasmadata(PlasmaData &pdata){
	Mo_Debug( "\tIn Model::update_plasmadata(PlasmaData &pdata)\n\n");
	Pdata = std::make_shared<PlasmaData>(pdata);
}

const bool Model::update_plasmadata(){
	Mo_Debug( "\tIn Model::update_plasmadata()\n\n");
	double ConvertJtoK(7.24297166e22);		// Conversion factor from ev to K
	bool InGrid = locate(i,k,Sample->get_position());
	if( !InGrid ) return InGrid;			// Particle has escaped simulation domain
	update_fields(i,k);
	Pdata->NeutralDensity 	= 1e19;  	// NEUTRAL DENSITY EQUALS 10e19
	Pdata->ElectronDensity 	= get_na1(i,k);  
	Pdata->IonDensity 		= get_na0(i,k);
	Pdata->IonTemp			= get_Ti(i,k)*ConvertJtoK;
	Pdata->ElectronTemp 	= get_Te(i,k)*ConvertJtoK;
	Pdata->NeutralTemp 		= 0.025*echarge*ConvertJtoK; 	// NEUTRAL TEMP EQUAL TO ION TEMP
	Pdata->AmbientTemp 		= 300; 						// NOTE THIS IS HARD CODED OHMEINGOD
	
	return true;
}

// CHECK THIS FUNCTION IS THE SAME AS IT WAS BEFORE!
void Model::update_fields(int i, int k){
	Mo_Debug( "\tIn Model::update_fields(int " << i << ", int " << k << ")\n\n");
	threevector vp, E, B;

	// Get Average plasma velocity
	double aveu(0.0);
	if( (get_na0(i,k)>0.0) || (get_na1(i,k)>0.0) ){
		aveu = (get_na0(i,k)*get_ua0(i,k)+get_na1(i,k)*get_ua1(i,k))
			/(get_na0(i,k)+get_na1(i,k));
	}
	else aveu = 0.0;

	// Read magnetic field in
	B.setx(get_bx(i,k));
	B.sety(get_by(i,k));
	B.setz(get_bz(i,k));

	// Plasma velocity is parallel to the B field
	vp.setx(aveu*(B.getunit().getx()));
	vp.sety(aveu*(B.getunit().gety()));
	vp.setz(aveu*(B.getunit().getz()));

	if(get_gridflag(i,k)==1){
		if( (get_gridflag(i+1,k)==1) && (get_gridflag(i-1,k)==1) ){
			E.setx(-(get_po(i+1,k)-get_po(i-1,k))/(2.0*get_dlx()));
		}else if(get_gridflag(i+1,k)==1){
			E.setx(-(get_po(i+1,k)-get_po(i,k))/get_dlx());
		}else if(get_gridflag(i-1,k)==1){
			E.setx(-(get_po(i,k)-get_po(i-1,k))/get_dlx());
		}else E.setx(0.0);
		if((get_gridflag(i,k+1)==1)&&(get_gridflag(i,k-1)==1)){
			E.setz(-(get_po(i,k+1)-get_po(i,k-1))/(2.0*get_dlz()));
		}else if(get_gridflag(i,k+1)==1){
			E.setz(-(get_po(i,k+1)-get_po(i,k))/get_dlz());
		}else if(get_gridflag(i,k-1)==1){
			E.setz(-(get_po(i,k)-get_po(i,k-1))/get_dlz());
		}else E.sety(0.0);
	}else{
		E.setx(0.0);
		E.setz(0.0);
	}

	Pdata->PlasmaVel		= vp;
	Pdata->Gravity 			= threevector(0.0,0.0,-9.8);
	Pdata->ElectricField	= E;
	Pdata->MagneticField	= B;
}

void Model::RecordPlasmadata(){
	PlasmaDataFile.open("Data/pd.txt",std::ofstream::app);
	PlasmaDataFile << "\n" << TotalTime 
			<< "\t" << i << "\t" << k << "\t" << Pdata->NeutralDensity << "\t" << Pdata->ElectronDensity 
			<< "\t" << Pdata->IonDensity << "\t" << Pdata->IonTemp << "\t" << Pdata->ElectronTemp 
			<< "\t" << Pdata->NeutralTemp << "\t" << Pdata->AmbientTemp << "\t" << Pdata->PlasmaVel 
			<< "\t" << Pdata->Gravity << "\t" << Pdata->ElectricField << "\t" << Pdata->MagneticField;
	PlasmaDataFile.close();
	PlasmaDataFile.clear();
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

// *************************************** PRINTING FUNCTIONS *************************************** //

// Print the inside and the outside of the tokamak
void Model::vtkcircle(double r, std::ofstream &fout){
	P_Debug("\tModel::vtkcircle(double r, std::ofstream &fout)\n\n");
	int i,imax;
	double phi, dphi = PI/100.0;
	imax = 201;
	fout << "# vtk DataFile Version 2.0" << std::endl << "plasma" << std::endl << "ASCII" << std::endl << std::endl
		<< "DATASET UNSTRUCTURED_GRID" << std::endl << "POINTS " << imax << " float" << std::endl;
	for(phi=0.0;phi<=2*PI+dphi;phi+=dphi)
	{
		fout << r*cos(phi) << '\t' << r*sin(phi) << '\t' << 0.0 << std::endl;
	}
	fout << "CELLS " << imax-1 << " " << 3*imax << std::endl;
	for(i=0;i<=imax-2;i++)
	{
		fout << 2 << " " << i << " " << i+1 << std::endl;
	}
	fout << "CELL_TYPES " << imax-1 << std::endl;
	for(i=0;i<=imax-2;i++)
	{
		fout << "3" << std::endl;
	}
}

void Model::vtktoroid(){
	P_Debug("\tModel::vtktoroid()\n\n");
	std::ofstream inner("output/innerplasma.vtk"),outer("output/outerplasma.vtk");
	vtkcircle(get_gridxmin(),inner);
	vtkcircle(get_gridxmax(),outer);
	inner.close();
	outer.close();
}

void Model::datadump(){
	P_Debug("\tModel::datadump()\n\n");
	if(get_device()=='p'){
		std::cout << "#i\tk\tr\tz\tna0\tna1\tpo\tua0\tua1\tbz\n";
		for(unsigned int i=0; i< get_gridx(); i++){
			for(unsigned int k=0; k< get_gridz(); k++){
				std::cout << i << "\t" << k << "\t" << i*0.15/get_gridx() << "\t" << k*1.0/get_gridz() << "\t" << get_na0(i,k) 
					<< "\t" << get_na1(i,k) << "\t" << get_po(i,k) << "\t" << get_ua0(i,k) << "\t" << get_ua1(i,k)
					<< "\t" << get_bz(i,k) << "\n";
			}
		}
	}else{
		vtktoroid();
	}
	//impurityprint(totalmass);
}

/*
void Model::impurityprint(double totalmass)
{
	int i,k;
	impurity << "# vtk DataFile Version 2.0" << std::endl << "//impurity" << std::endl << "ASCII" << std::endl << std::endl
		<< "DATASET STRUCTURED_GRID" << std::endl << "DIMENSIONS " << gridx << " " << gridz << " " << 1 << std::endl
		<< "POINTS " << gridx*gridz << " float" << std::endl;
	for(k=0;k<=gridz-1;k++)
	{
		for(i=0;i<=gridx-1;i++)
		{
			impurity << x[i][k] << " " << z[i][k] << " " << 0.0 << std::endl;
		}
	}

	// Print cell data
	impurity << std::endl << "POINT_DATA " << gridx*gridz << std::endl << "SCALARS " << "//impurity " << "float" << std::endl
		<< "LOOKUP_TABLE " << "default" << std::endl;
	for(k=0;k<=gridz-1;k++)
	{
		for(i=0;i<=gridx-1;i++)
		{
			impurity << mevap[i][k]/totalmass << std::endl;
		}
	}
}
*/
// *************************************** PRINTING FUNCTIONS *************************************** //

/*
void Model::impurityprint(double totalmass)
{
	int i,k;
	impurity << "# vtk DataFile Version 2.0" << std::endl << "//impurity" << std::endl << "ASCII" << std::endl << std::endl
		<< "DATASET STRUCTURED_GRID" << std::endl << "DIMENSIONS " << gridx << " " << gridz << " " << 1 << std::endl
		<< "POINTS " << gridx*gridz << " float" << std::endl;
	for(k=0;k<=gridz-1;k++)
	{
		for(i=0;i<=gridx-1;i++)
		{
			impurity << x[i][k] << " " << z[i][k] << " " << 0.0 << std::endl;
		}
	}

	// Print cell data
	impurity << std::endl << "POINT_DATA " << gridx*gridz << std::endl << "SCALARS " << "//impurity " << "float" << std::endl
		<< "LOOKUP_TABLE " << "default" << std::endl;
	for(k=0;k<=gridz-1;k++)
	{
		for(i=0;i<=gridx-1;i++)
		{
			impurity << mevap[i][k]/totalmass << std::endl;
		}
	}
}
*/