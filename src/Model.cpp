//#define PAUSE
//#define MODEL_DEBUG

#include "Model.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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
	Mo_Debug("\n\nIn Model::Model( Matter *&sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata, float accuracy ):Sample(sample),PG_data(pgrid),Pdata(PlasmaDefaults),Accuracy(accuracy), ContinuousPlasma(false)\n\n");
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

// Interpolate between grid points to determine plasma data
const double Model::interpolatepdata(const int i,const int k)const{
	Mo_Debug( "\tIn Model::interpolatepdata(const int i,const int k)const\n\n");
	threevector vp1, vp2, E, B, gravity(0.0,0.0,-9.81);

	if( PG_data->dlx != PG_data->dlz ){
		static bool runOnce = true;
		WarnOnce(runOnce,"PlasmaGrid Interpolation only valid for square Grid! PG_data->dlx != PG_data->dlz!");
	}

	// Get the position of the dust as a decimal number of grid cells
	double i_pos = 0.5+(Sample->get_position().getx()-PG_data->gridxmin)/PG_data->dlx;
	double k_pos = 0.5+(Sample->get_position().getz()-PG_data->gridzmin)/PG_data->dlz;

	// Calculate the distance to the relative (0,0) of the current square
	double i_diff = i_pos - i;
	double k_diff = k_pos - k;

	// Calculate the SI distance to each edge (m)
	double dx_1 = PG_data->dlx*(1-i_diff);
	double dz_1 = PG_data->dlz*(1-k_diff);
	double dx_2 = PG_data->dlx*i_diff;
	double dz_2 = PG_data->dlz*k_diff;

	// Calculate the normalisation
	double Coeff = 1.0/sqrt(4.0*PG_data->dlx*PG_data->dlx+4.0*PG_data->dlz*PG_data->dlz);
	
	// If it is not edge, in which case we aren't in a square
	if( PG_data->gridx >= i_pos && PG_data->gridz >= k_pos 
		&& 0 <= i_pos && 0 <= k_pos ){
		/*std::cout << "\n\nCoeff = " << Coeff;
		std::cout << "\ni = " << i;
		std::cout << "\nk = " << k;
		std::cout << "\ni_pos = " << i_pos;
		std::cout << "\nk_pos = " << k_pos; 
		std::cout << "\ni_diff = " << i_diff;
		std::cout << "\nk_diff = " << k_diff;
		std::cout << "\ni+sgn(i_diff) = " << i+1*sgn(i_diff);
		std::cout << "\nk+sgn(k_diff) = " << k+1*sgn(k_diff);
		std::cout << "\nPG_data->na0[" << i << "][" << k << "] = " << PG_data->na0[i][k]; 
		std::cout << "\nPG_data->na0[" << i+sgn(i_diff) << "][" << k << "] = " << PG_data->na0[i+1*sgn(i_diff)][k]; 
		std::cout << "\nPG_data->na0[" << i << "][" << k+sgn(k_diff) << "] = " << PG_data->na0[i][k+1*sgn(k_diff)]; 
		std::cout << "\nPG_data->na0[" << i+sgn(i_diff) << "][" << k+sgn(k_diff) << "] = " << PG_data->na0[i+1*sgn(i_diff)][k+1*sgn(k_diff)];
		std::cout << "\nsqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->na0[" << i << "][" << k << "] = " << sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->na0[i][k]; 
		std::cout << "\nsqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->na0[" << i+sgn(i_diff) << "][" << k << "] = " << sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->na0[i+1*sgn(i_diff)][k]; 
		std::cout << "\nsqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->na0[" << i << "][" << k+sgn(k_diff) << "] = " << sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->na0[i][k+1*sgn(k_diff)]; 
		std::cout << "\nsqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->na0[" << i+sgn(i_diff) << "][" << k+sgn(k_diff) << "] = " << sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->na0[i+1*sgn(i_diff)][k+1*sgn(k_diff)]; std::cin.get();
		std::cout << "\nPdata->IonDensity = " << Pdata->IonDensity;*/
		Pdata->NeutralDensity 	= Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->na2[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->na2[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->na2[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->na2[i+1*sgn(i_diff)][k+1*sgn(k_diff)]);
		
		Pdata->ElectronDensity 	= Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->na1[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->na1[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->na1[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->na1[i+1*sgn(i_diff)][k+1*sgn(k_diff)]);

		Pdata->IonDensity	 	= Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->na0[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->na0[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->na0[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->na0[i+1*sgn(i_diff)][k+1*sgn(k_diff)]);
		

		Pdata->IonTemp	 		= Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->Ti[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->Ti[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->Ti[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->Ti[i+1*sgn(i_diff)][k+1*sgn(k_diff)]);

		Pdata->ElectronTemp 	= Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->Te[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->Te[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->Te[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->Te[i+1*sgn(i_diff)][k+1*sgn(k_diff)]);

		Pdata->NeutralTemp 		= Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->Tn[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->Tn[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->Tn[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->Tn[i+1*sgn(i_diff)][k+1*sgn(k_diff)]);

		Pdata->AmbientTemp	 	= Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->Ta[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->Ta[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->Ta[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->Ta[i+1*sgn(i_diff)][k+1*sgn(k_diff)]);
		B.setx(	Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->bx[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->bx[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->bx[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->bx[i+1*sgn(i_diff)][k+1*sgn(k_diff)]));

		B.sety(	Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->by[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->by[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->by[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->by[i+1*sgn(i_diff)][k+1*sgn(k_diff)]));

		B.setz(	Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->bz[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->bz[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->bz[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->bz[i+1*sgn(i_diff)][k+1*sgn(k_diff)]));

		// Ion Plasma velocity parallel to the B field
		vp1.setx(Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->ua0[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->ua0[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->ua0[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->ua0[i+1*sgn(i_diff)][k+1*sgn(k_diff)])*(B.getunit().getx()));
		vp1.sety(Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->ua0[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->ua0[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->ua0[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->ua0[i+1*sgn(i_diff)][k+1*sgn(k_diff)])*(B.getunit().gety()));
		vp1.setz(Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->ua0[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->ua0[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->ua0[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->ua0[i+1*sgn(i_diff)][k+1*sgn(k_diff)])*(B.getunit().getz()));

		// Electron Plasma velocity parallel to the B field
		vp2.setx(Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->ua1[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->ua1[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->ua1[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->ua1[i+1*sgn(i_diff)][k+1*sgn(k_diff)])*(B.getunit().getx()));
		vp2.sety(Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->ua1[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->ua1[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->ua1[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->ua1[i+1*sgn(i_diff)][k+1*sgn(k_diff)])*(B.getunit().gety()));
		vp2.setz(Coeff*(sqrt(dx_1*dx_1+dz_1*dz_1)*PG_data->ua1[i][k]
								+sqrt(dx_2*dx_2+dz_1*dz_1)*PG_data->ua1[i+1*sgn(i_diff)][k]
								+sqrt(dx_1*dx_1+dz_2*dz_2)*PG_data->ua1[i][k+1*sgn(k_diff)]
								+sqrt(dx_2*dx_2+dz_2*dz_2)*PG_data->ua1[i+1*sgn(i_diff)][k+1*sgn(k_diff)])*(B.getunit().getz()));

	}else{
		Pdata->NeutralDensity 	= PG_data->na2[i][k];  
		Pdata->ElectronDensity 	= PG_data->na1[i][k];  
		Pdata->IonDensity 		= PG_data->na0[i][k];
		Pdata->IonTemp			= PG_data->Ti[i][k];
		Pdata->ElectronTemp 	= PG_data->Te[i][k];
		Pdata->NeutralTemp		= PG_data->Tn[i][k];
		Pdata->AmbientTemp		= PG_data->Ta[i][k];

		B.setx(PG_data->bx[i][k]);
		B.sety(PG_data->by[i][k]);
		B.setz(PG_data->bz[i][k]);

		vp1.setx(PG_data->ua0[i][k]*(B.getunit().getx()));
		vp1.sety(PG_data->ua0[i][k]*(B.getunit().gety()));
		vp1.setz(PG_data->ua0[i][k]*(B.getunit().getz()));

		vp2.setx(PG_data->ua1[i][k]*(B.getunit().getx()));
		vp2.sety(PG_data->ua1[i][k]*(B.getunit().gety()));
		vp2.setz(PG_data->ua1[i][k]*(B.getunit().getz()));
	}

	E.setx(-(PG_data->po[i+1][k]-PG_data->po[i-1][k])/(2.0*PG_data->dlx));
	E.sety(0.0);
	E.setz(-(PG_data->po[i][k+1]-PG_data->po[i][k-1])/(2.0*PG_data->dlz));

	
	// Setup for Magnum-PSI
	// For Magnum PSI, Gravity is not in -z direction but has a radial & Azimuthal
	if( PG_data->device == 'p' ){ 
		double Theta = Sample->get_position().gety();
		gravity.setx(-1.0*gravity.mag3()*cos(Theta));
		gravity.sety(gravity.mag3()*sin(Theta));
		gravity.setz(0.0);
	}else{
		gravity = Pdata->Gravity;
	}
	Pdata->PlasmaVel		= vp1;
	Pdata->ElectricField	= E;
	Pdata->MagneticField	= B;
	Pdata->Gravity			= gravity;
	return true;
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
	//interpolatepdata(i,k); //RecordPlasmadata("pd.txt");
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

const double Model::SOMLIonFlux(double Potential)const{
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	if( uz == 0.0 ){
		return OMLIonFlux(Potential);
	}
	if( Potential >= 0.0 ){
		double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+exp(-uz*uz)/2.0;
		double s2 = sqrt(PI)*erf(uz)/(2.0*uz);
		return Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*
				(s1+s2*Pdata->Z*Potential/Tau);
	}else{
		double uzp = uz+sqrt(-Pdata->Z*Potential/Tau);
		double uzm = uz-sqrt(-Pdata->Z*Potential/Tau);
		return Pdata->IonDensity*IonThermalVelocity*(1.0/(8.0*sqrt(2.0)*uz))*
				((1.0+2.0*(uz*uz+Pdata->Z*Potential/Tau))*(erf(uzp)+erf(uzm))
				+(2.0/sqrt(PI))*(uzp*exp(-uzm*uzm)+uzm*exp(-uzp*uzp)));
	}
}

const double Model::SMOMLIonFlux(double Potential)const{
	double HeatCapacityRatio = 5.0/3.0;
	double MassRatio = Pdata->mi/Me;
	double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
	double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
	double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
	if( uz == 0.0 ){ // Handle zero-relative velocity case
		return Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0))*(1.0+(1.0/Tau)*
			(Potential*Pdata->Z-0.5*log(2.0*PI*(1.0+HeatCapacityRatio*Tau)/MassRatio)));
	}
	if( Potential >= 0.0 ){ // For negative dust, do SMOML
		double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+exp(-uz*uz)/2.0;
		double s2 = sqrt(PI)*erf(uz)/(2.0*uz);
		return Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*(s1+(s2/Tau)*
			(Potential*Pdata->Z-0.5*log(2.0*PI*(1.0+HeatCapacityRatio*Tau)/MassRatio)));
	}else{	// For Positive dust, do SOML
		double uzp = uz+sqrt(-Pdata->Z*Potential/Tau);
		double uzm = uz-sqrt(-Pdata->Z*Potential/Tau);
		return Pdata->IonDensity*IonThermalVelocity*(1.0/(8.0*sqrt(2.0)*uz))*
				((1.0+2.0*(uz*uz+Pdata->Z*Potential/Tau))*(erf(uzp)+erf(uzm))
				+(2.0/sqrt(PI))*(uzp*exp(-uzm*uzm)+uzm*exp(-uzp*uzp)));
	}
}

const double Model::PHLElectronFlux(double Potential)const{
	double Tau = Pdata->ElectronTemp/Pdata->IonTemp;
	double Beta = Sample->get_radius()
			/(sqrt(PI*Pdata->ElectronTemp*Me)/(2.0*echarge*echarge*Pdata->MagneticField*Pdata->MagneticField));
	double MassRatio = Pdata->mi/Me;

	if( Beta/MassRatio > 0.01 ){
		static bool runOnce = true;
		WarnOnce(runOnce,"Beta/MassRatio > 0.01 in solvePHL! Model may not be valid in this range! see Fig 11. of  L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).");
	}
	double AtomicNumber = Pdata->Z;	
	double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/(Pdata->ElectronDensity*pow(echarge,2)));

	double z = Beta/(1.0+Beta);
	double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+1.150*z*z*z*z*z;
	double eta = (Potential/Beta)*(1.0+(Beta/4.0)*(1-exp(-4.0/(DebyeLength*Beta))));

	double w(1.0);
	if( Beta == 0.0 || eta == -1.0 ){
		w = 1.0;
	}else if( std::isnan(eta) ){
		std::cout << "\nWarning! w being set to 1.0 (Assuming high B field limit) but Phi/Beta is nan while Beta != 0.";
		w = 1.0;
	}else{
		w = eta/(1+eta);
	} 
	
	double A = 0.678*w+1.543*w*w-1.212*w*w*w;

	if( Potential >= 0.0 ){ // For negative dust, do PHL
		return Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(A+(1.0-A)*i_star)*exp(-Potential);
	}else{ // For positive dust, do OML
		return Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(1.0-Potential);
	}
}

const double Model::DTOKSIonFlux(double Potential)const{
	H_Debug("\n\tIn Model::IonFlux():");

	double IonFlux=0;

	if( Sample->is_positive() ) IonFlux = DTOKSElectronFlux(Potential); //Positive grain, DeltaTot() > 1
	else	IonFlux = DTOKSElectronFlux(Potential)*(1-Sample->get_deltatot());
	assert(IonFlux >= 0);
	return IonFlux;
}

const double Model::DTOKSElectronFlux(double Potential)const{
	H_Debug("\n\tIn Model::ElectronFlux():\n\n");
	return Pdata->ElectronDensity*exp(-Potential)*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
}

const double Model::OMLIonFlux(double Potential)const{
	H_Debug("\n\tIn Model::IonFlux():");

	double IonFlux=0;

	if( Potential >= 0 ) 
		IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
				*(1+Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
	else	
		IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
				*exp(Potential*(Pdata->ElectronTemp/Pdata->IonTemp));

	assert(IonFlux >= 0);
	return IonFlux;
}

const double Model::OMLElectronFlux(double Potential)const{
	H_Debug("\n\tIn Model::ElectronFlux():\n\n");
	if( Potential < 0.0 ){
		return Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me))*(1-Potential);
	}else{ 	
		return Pdata->ElectronDensity*exp(-Potential)*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
}	}

const double Model::NeutralFlux()const{
	H_Debug("\n\tIn Model::NeutralFlux():\n\n");

	return Pdata->NeutralDensity*sqrt(Kb*Pdata->NeutralTemp/(2*PI*Pdata->mi));
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
