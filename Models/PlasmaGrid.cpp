#define PLASMAGRID_Debug

#include "PlasmaGrid.h"

PlasmaGrid::PlasmaGrid(std::string filename, char element, char machine, double xspacing, double zspacing):mi(Mp),gamma(Me/mi),gas(element),device(machine),dlx(xspacing),dlz(zspacing){
	P_Debug( "\n\nIn PlasmaGrid::PlasmaGrid(std::string filename, char element, char machine, double xspacing, double zspacing):mi(Mp),gamma(Me/mi),gas(element),device(machine),dlx(xspacing),dlz(zspacing)\n\n");

	// Plasma parameters
	if(device=='m'){
		gridx = 121;
		gridz = 281;
		gridtheta=0;
		gridxmin = 0.2;
		gridzmin = -2.0;
		gridxmax = 1.4;
		gridzmax = 0.80;
		std::cout << "\n\n\tCalculation for MAST" << std::endl;
	}else if(device=='i'){
		gridx = 451;
		gridz = 951;
		gridtheta=0;
		gridxmin = 4.0;
		gridzmin = -4.7;
		gridxmax = 8.25;
		gridzmax = 4.80;
		std::cout << "\n\n\tCalculation for ITER" << std::endl;
	}else if(device=='d'){
		gridx = 180;
		gridz = 400;
		gridtheta=0;
		gridxmin = 0.2;
		gridzmin = -2.0;
		gridxmax = 2.0;
		gridzmax = 4.0;

		std::cout << "\n\n\tCalculation for Double Null MAST 17839files shot" << std::endl;
	}else if(device=='p'){
		gridx = 64;
		gridz = 20;
		gridtheta=64;
		gridxmin = 0.0;
		gridzmin = 0.0;
		gridxmax = 0.15;
//		gridzmax = 1.0;
		gridzmax = 1.9;
		std::cout << "\n\n\t* Calculation for Magnum-PSI *\n";
	}else std::cout << "Invalid tokamak" << std::endl;

	Te	= std::vector<std::vector<double>>(gridx,std::vector<double>(gridz));
	Ti 	= Te;
	na0	= Te;
	na1	= Te;
	po 	= Te;
	ua0	= Te;
	ua1	= Te;
	bx 	= Te;
	by 	= Te;
	bz 	= Te;
	x  	= Te;
	z 	= Te;
	mevap 	= Te;
	gridflag= std::vector<std::vector<int>>(gridx,std::vector<int>(gridz));
	std::cout << "\n\t* Begin File read *";
	//impurity.open("output///impurity///impurity.vtk");
	try{ // Read data
		int readstatus = readdata(filename);
		if(readstatus != 0)
			throw readstatus;
	}	
	catch( int e ){
		std::cerr << "Exception opening/reading/closing file\nError status: " << e;
	}
	// Everything should now be initialised correctly...
}

// *************************************** READING FUNCTIONS *************************************** //

// Read an input file
void PlasmaGrid::readscalars(std::ifstream &input){
	P_Debug("\tPlasmaGrid::readscalars(std::ifstream &input)\n\n");
	int i,k;
	char dummy;
	// Ignore first line of file
	for(i=0;i<=19;i++){
		input >> dummy;
	}
	for(k=0;k<=gridz-1;k++){
		for(i=0;i<=gridx-1;i++){
			input >> x[i][k] >> z[i][k] >> Te[i][k] >> Ti[i][k] >> na0[i][k] >> na1[i][k] >> po[i][k] >> ua0[i][k]
				>> ua1[i][k];
			// For some reason, very small non-zero values are being assigned to the 'zero' values being read in.
			// This should correct for this...
			if( Te[i][k] != 0 && Te[i][k] < 1e-100 ) {	Te[i][k]  = 0.0; } 
			if( Ti[i][k] != 0 && Ti[i][k] < 1e-100 ) {	Ti[i][k]  = 0.0; } 
			if( na0[i][k]!= 0 && na0[i][k] < 1e-100 ){ 	na0[i][k] = 0.0; }
			if( na1[i][k]!= 0 && na1[i][k] < 1e-100 ){ 	na1[i][k] = 0.0; }
			if( po[i][k] != 0 && po[i][k] < 1e-100 ) {	po[i][k]  = 0.0; } 
			if( ua0[i][k]!= 0 && ua0[i][k] < 1e-100 ){ 	ua0[i][k] = 0.0; }
			if( ua1[i][k]!= 0 && ua1[i][k] < 1e-100 ){ 	ua1[i][k] = 0.0; }
			mevap[i][k] = 0.0;
		}
	}
}

void PlasmaGrid::readthreevectors(std::ifstream &input){
	P_Debug("\tPlasmaGrid::readthreevectors(std::ifstream &input)\n\n");
	int i,k;
	char dummy;
	double dummy1, dummy2;
	// Ignore first line of file
	for(i=0;i<=19;i++){
		input >> dummy;
	}
	for(k=0;k<=gridz-1;k++){
		for(i=0;i<=gridx-1;i++){
			input >> dummy1 >> dummy2 >> bx[i][k] >> bz[i][k] >> by[i][k];
			bz[i][k] = -bz[i][k];
			if( bx[i][k] != 0 && bx[i][k] < 1e-100 ) {	bx[i][k] = 0.0; } 
			if( by[i][k] != 0 && by[i][k] < 1e-100 ) {	by[i][k] = 0.0; } 
			if( bz[i][k] != 0 && bz[i][k] < 1e-100 ) { 	bz[i][k] = 0.0; }
		}
	}
}

void PlasmaGrid::readgridflag(std::ifstream &input){
	P_Debug("\tPlasmaGrid::readgridflag(std::ifstream &input)\n\n");
	int i,k;
	double dummy1,dummy2;
	for(k=0;k<=gridz-1;k++){
		for(i=0;i<=gridx-1;i++){
			input >> dummy1 >> dummy2 >> gridflag[i][k];
		}
	}
}


int PlasmaGrid::readMPSIdata(std::string filename){
	P_Debug("\tPlasmaGrid::readMPSIdata(std::string filename)\n\n");
	
	const int NC_ERR = 2;
	std::cout << "\n\t\t* Reading data from: " << filename << " *";

	float electron_dens_mat[gridx][gridz][gridtheta];
	float electron_temp_mat[gridx][gridz][gridtheta];
	float electron_Vele_mat[gridx][gridz][gridtheta];
	float ion_dens_mat[gridx][gridz][gridtheta];
	float ion_Veli_mat[gridx][gridz][gridtheta];
	float Potential_mat[gridx][gridz][gridtheta];
	float Bxy_mat[gridx][gridz];

	// Change the error behavior of the netCDF C++ API by creating an
	// NcError object. Until it is destroyed, this NcError object will
	// ensure that the netCDF C++ API silently returns error codes on
	// any failure, and leaves any other error handling to the calling
	// program. In the case of this example, we just exit with an
	// NC_ERR error code.
	NcError err(NcError::silent_nonfatal);

	
	NcFile dataFile(filename.c_str(), NcFile::ReadOnly);
	if(!dataFile.is_valid())
		return NC_ERR;
	
	if (dataFile.num_dims() != 3 || dataFile.num_vars() != 8 ||
		dataFile.num_atts() != 0 || dataFile.rec_dim() != 0)
		return NC_ERR;
	// We get back a pointer to each NcVar we request. Get the
	// latitude and longitude coordinate variables.

	NcVar *Ne, *e_Temp, *Vel_e, *Ni, *Vel_i, *Potential; //, *B_xy;
	if (!(Ne = dataFile.get_var("Ne")))
		return NC_ERR;
	if (!(e_Temp = dataFile.get_var("Te")))
		return NC_ERR;
	if (!(Vel_e = dataFile.get_var("Vel_e")))
		return NC_ERR;
	if (!(Ni = dataFile.get_var("Ni")))
		return NC_ERR;
	if (!(Vel_i = dataFile.get_var("Vel_i")))
		return NC_ERR;
	if (!(Potential = dataFile.get_var("Potential_Phi")))
		return NC_ERR;
//	if (!(B_xy = dataFile.get_var("B_xy")))
//		return NC_ERR;

	if (!Ne->get(&electron_dens_mat[0][0][0], gridx, gridz, gridtheta))
		return NC_ERR;
	if (!e_Temp->get(&electron_temp_mat[0][0][0], gridx, gridz, gridtheta))
		return NC_ERR;
	if (!Vel_e->get(&electron_Vele_mat[0][0][0], gridx, gridz, gridtheta))
		return NC_ERR;
	if (!Ni->get(&ion_dens_mat[0][0][0], gridx, gridz, gridtheta))
		return NC_ERR;
	if (!Vel_i->get(&ion_Veli_mat[0][0][0], gridx, gridz, gridtheta))
		return NC_ERR;
	if (!Potential->get(&Potential_mat[0][0][0], gridx, gridz, gridtheta))
		return NC_ERR;
//	if (!B_xy->get(&Bxy_mat[0][0], gridx, gridz))
//		return NC_ERR;
	
	for(unsigned int i=0; i< gridx; i++){
		for(unsigned int k=0; k< gridz; k++){
			Ti[i][k] = 1.0*echarge;
			Te[i][k] = electron_temp_mat[i][k][0];
			na0[i][k] = ion_dens_mat[i][k][0];
			na1[i][k] = electron_dens_mat[i][k][0];
			po[i][k] = Potential_mat[i][k][0];
			ua0[i][k] = ion_Veli_mat[i][k][0];
			ua1[i][k] = electron_Vele_mat[i][k][0];
			bx[i][k] = 0.0;
			by[i][k] = 0.0;
			bz[i][k] = 0.4;
//			bz[i][k] = Bxy_mat[i][k];
		}
	}
	return 0;

}

int PlasmaGrid::readdata(std::string filename){
	P_Debug("\tPlasmaGrid::readdatastd::string filename()\n\n");
	// Input files
	if(device=='p'){ // Note, grid flags will be empty 
		return readMPSIdata(filename);
	}else{
		std::cout << "\nWarning! Filename: " << filename << " is unused parameter in function PlasmaGrid::readdata(std::string filename)";
		std::ifstream scalars,threevectors,gridflagfile;
		if(device=='m'){
			scalars.open("Models/PlasmaData/MAST/b2processed.dat");
			threevectors.open("Models/PlasmaData/MAST/b2processed2.dat");
			gridflagfile.open("Models/PlasmaData/MAST/locate.dat");
		}else if(device=='i'){
			scalars.open("Models/PlasmaData/regulardataITER/b2processed.dat");
			threevectors.open("Models/PlasmaData/regulardataITER/b2processed2.dat");
			gridflagfile.open("Models/PlasmaData/regulardataITER/locate.dat");
		}else if(device=='d'){
			scalars.open("Models/PlasmaData/MASTDexp/b2processed.dat");
			threevectors.open("Models/PlasmaData/MASTDexp/b2processed2.dat");
			gridflagfile.open("Models/PlasmaData/MASTDexp/locate.dat");
		}
		assert( scalars.is_open() );
		assert( threevectors.is_open() );
		assert( gridflagfile.is_open() );
		readscalars(scalars);
		readthreevectors(threevectors);
		readgridflag(gridflagfile);
		scalars.close();
		threevectors.close();
		gridflagfile.close();
	}
	return 0;
}

// *************************************** MUTATION FUNCTIONS *************************************** //

// Locate dust particle in the plasma grid
bool PlasmaGrid::locate(int &i, int &k, const threevector xd)const{
	P_Debug("\tIn PlasmaGrid::locate(int &" << i << ", int &" << k << ", " << xd << ")\n\n");
	// Adding 0.5 makes the rounding work properly
	i = int(0.5+(xd.getx()-gridxmin)/dlx);
	k = int(0.5+(xd.getz()-gridzmin)/dlz);
	if( xd.getx() < gridxmin )
		i = -1;
	if( xd.getz() < gridzmin )
		k = -1;
	return checkingrid(i,k);

}

bool PlasmaGrid::checkingrid(const int i, const int k)const{
	P_Debug("\tIn PlasmaGrid::checkingrid(int " << i << ", int " << k  << ")\n\n");
	bool returnval(true);
	if( i >= gridx || i < 0) returnval = false;
	if( k >= gridz || k < 0) returnval = false;
	return returnval;
}

// *************************************** PRINTING FUNCTIONS *************************************** //

// Print the inside and the outside of the tokamak
void PlasmaGrid::vtkcircle(double r, std::ofstream &fout){
	P_Debug("\tPlasmaGrid::vtkcircle(double r, std::ofstream &fout)\n\n");
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

void PlasmaGrid::vtktoroid(){
	P_Debug("\tPlasmaGrid::vtktoroid()\n\n");
	std::ofstream inner("output/innerplasma.vtk"),outer("output/outerplasma.vtk");
	vtkcircle(gridxmin,inner);
	vtkcircle(gridxmax,outer);
	inner.close();
	outer.close();
}
/*
void PlasmaGrid::impurityprint(double totalmass)
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
void PlasmaGrid::datadump(){
	P_Debug("\tPlasmaGrid::datadump()\n\n");
	if(device=='p'){
		std::cout << "#i\tk\tr\tz\tna0\tna1\tpo\tua0\tua1\tbz\n";
		for(unsigned int i=0; i< gridx; i++){
			for(unsigned int k=0; k< gridz; k++){
				std::cout << i << "\t" << k << "\t" << i*0.15/gridx << "\t" << k*1.0/gridz << "\t" << na0[i][k] 
					<< "\t" << na1[i][k] << "\t" << po[i][k] << "\t" << ua0[i][k] << "\t" << ua1[i][k] 
					<< "\t" << bz[i][k] << "\n";
			}
		}
	}else{
		vtktoroid();
	}
	//impurityprint(totalmass);

}
