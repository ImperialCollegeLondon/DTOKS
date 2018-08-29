//#define PAUSE
//#define MODEL_DEBUG

#include "Model.h"

Model::Model():Sample(new Tungsten),PG_data(std::make_shared<PlasmaGrid_Data>(PlasmaGrid_DataDefaults)),Pdata(&PlasmaDataDefaults),Accuracy(1.0),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
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
		:Sample(sample),PG_data(std::make_shared<PlasmaGrid_Data>(PlasmaGrid_DataDefaults)),Pdata(std::make_shared<PlasmaData>(pdata)),Accuracy(accuracy),ContinuousPlasma(true),TimeStep(0.0),TotalTime(0.0){
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
	static bool runOnce = true;
	WarnOnce(runOnce,"Default values being taken: Tn = 0.025*116045.25K, Nn = 1e19m^-3, Ta = 300K & Mi = 1.66054e-27Kg!");
	PlasmaDataFile.close();
	PlasmaDataFile.clear();
	update_plasmadata();
}

// Constructor for Matter sample moving in a varying plasma background given by PlasmaGrid (pgrid) with the current information 
// stored in pdata.
Model::Model( Matter *&sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata, float accuracy )
		:Sample(sample),PG_data(std::make_shared<PlasmaGrid_Data>(pgrid)),Pdata(std::make_shared<PlasmaData>(pdata)),Accuracy(accuracy),ContinuousPlasma(false),TimeStep(0.0),TotalTime(0.0){
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
	i = int(0.5+(xd.getx()-PG_data->gridxmin)/PG_data->dlx);
	k = int(0.5+(xd.getz()-PG_data->gridzmin)/PG_data->dlz);
	if( xd.getx() < PG_data->gridxmin )
		i = -1;
	if( xd.getz() < PG_data->gridzmin )
		k = -1;
	return checkingrid(i,k);

}

const bool Model::checkingrid(const int i, const int k)const{
	P_Debug("\tIn Model::checkingrid(int " << i << ", int " << k  << ")\n\n");
	bool returnval(true);
	if( i >= PG_data->gridx || i < 0) returnval = false;
	if( k >= PG_data->gridz || k < 0) returnval = false;
	if( !ContinuousPlasma )
		return returnval;

	return true;
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
	
	bool InGrid = locate(i,k,Sample->get_position());
	if( !InGrid ) return InGrid;			// Particle has escaped simulation domain
	if( ContinuousPlasma ) return InGrid;
	update_fields(i,k);

	Pdata->NeutralDensity 	= PG_data->na2[i][k];  
	Pdata->ElectronDensity 	= PG_data->na1[i][k];  
	Pdata->IonDensity 		= PG_data->na0[i][k];
	Pdata->IonTemp			= PG_data->Ti[i][k];
	Pdata->ElectronTemp 	= PG_data->Te[i][k];
	Pdata->NeutralTemp		= PG_data->Tn[i][k];
	Pdata->AmbientTemp		= PG_data->Ta[i][k];
	
	return true;
}

// CHECK THIS FUNCTION IS THE SAME AS IT WAS BEFORE!
void Model::update_fields(int i, int k){
	Mo_Debug( "\tIn Model::update_fields(int " << i << ", int " << k << ")\n\n");
	threevector vp, E, B, gravity(0.0,0.0,-9.81);

	// Get Average plasma velocity
	double aveu(0.0);
	if( (PG_data->na0[i][k]>0.0) || (PG_data->na1[i][k]>0.0) ){
		aveu = (PG_data->na0[i][k]*PG_data->ua0[i][k]+PG_data->na1[i][k]*PG_data->ua1[i][k])
			/(PG_data->na0[i][k]+PG_data->na1[i][k]);
	}
	else aveu = 0.0;

	// Read magnetic field in
	B.setx(PG_data->bx[i][k]);
	B.sety(PG_data->by[i][k]);
	B.setz(PG_data->bz[i][k]);

	// Plasma velocity is parallel to the B field
	vp.setx(aveu*(B.getunit().getx()));
	vp.sety(aveu*(B.getunit().gety()));
	vp.setz(aveu*(B.getunit().getz()));

	if(PG_data->gridflag[i][k]==1){
		if( (PG_data->gridflag[i+1][k]==1) && (PG_data->gridflag[i-1][k]==1) ){
			E.setx(-(PG_data->po[i+1][k]-PG_data->po[i-1][k])/(2.0*PG_data->dlx));
		}else if(PG_data->gridflag[i+1][k]==1){
			E.setx(-(PG_data->po[i+1][k]-PG_data->po[i][k])/PG_data->dlx);
		}else if(PG_data->gridflag[i-1][k]==1){
			E.setx(-(PG_data->po[i][k]-PG_data->po[i-1][k])/PG_data->dlx);
		}else E.setx(0.0);
		if((PG_data->gridflag[i][k+1]==1)&&(PG_data->gridflag[i][k-1]==1)){
			E.setz(-(PG_data->po[i][k+1]-PG_data->po[i][k-1])/(2.0*PG_data->dlz));
		}else if(PG_data->gridflag[i][k+1]==1){
			E.setz(-(PG_data->po[i][k+1]-PG_data->po[i][k])/PG_data->dlz);
		}else if(PG_data->gridflag[i][k-1]==1){
			E.setz(-(PG_data->po[i][k]-PG_data->po[i][k-1])/PG_data->dlz);
		}else E.sety(0.0);
	}else{
		E.setx(0.0);
		E.setz(0.0);
	}

	double Theta = Sample->get_position().gety();
	// Setup for Magnum-PSIq
	
	if( PG_data->device == 'p' ){ // For Magnum PSI, Gravity is not in -z direction but is radial & Azimuthal
		gravity.setx(-1.0*gravity.mag3()*cos(Theta));
		gravity.sety(gravity.mag3()*sin(Theta));
		gravity.setz(0.0);
	}else{
		gravity = Pdata->Gravity;
	}

	Pdata->PlasmaVel		= vp;
	Pdata->ElectricField	= E;
	Pdata->MagneticField	= B;
	Pdata->Gravity			= gravity;
}

void Model::RecordPlasmadata(std::string filename){
	PlasmaDataFile.open("Data/" + filename,std::ofstream::app);
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

	if( Pdata->IonTemp <= 0 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"\nWarning!  Pdata->IonTemp <= 0!");
		return 0.0;
	}
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
	vtkcircle(PG_data->gridxmin,inner);
	vtkcircle(PG_data->gridxmax,outer);
	inner.close();
	outer.close();
}

void Model::datadump(){
	P_Debug("\tModel::datadump()\n\n");
	if(PG_data->device=='p'){
		std::cout << "#i\tk\tr\tz\tna0\tna1\tpo\tua0\tua1\tbz\n";
		for(unsigned int i=0; i< PG_data->gridx; i++){
			for(unsigned int k=0; k< PG_data->gridz; k++){
				std::cout << i << "\t" << k << "\t" << i*0.15/PG_data->gridx << "\t" << k*1.0/PG_data->gridz << "\t" << PG_data->na0[i][k] 
					<< "\t" << PG_data->na1[i][k] << "\t" << PG_data->na2[i][k] << "\t" << PG_data->po[i][k] << "\t" << PG_data->ua0[i][k] 
					<< "\t" << PG_data->ua1[i][k] << "\t" << PG_data->bz[i][k] << "\n";
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
