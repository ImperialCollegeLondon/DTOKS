#include "PlasmaGrid.h"

PlasmaGrid::PlasmaGrid(char element, char machine, double spacing):mi(Mp),gamma(Me/mi),gas(element),device(machine),dl(spacing){
	std::cout << "\n\nIn PlasmaGrid::PlasmaGrid(char element, char machine, double spacing):mi(Mp),gamma(Me/mi),gas(element),device(machine),dl(spacing)\n\n";
	
	// Plasma parameters
	if(device=='m'){
		gridx = 121;
		gridz = 281;
		gridxmin = 0.2;
		gridzmin = -2.0;
		rmin = 0.2;
		rmax = 1.4;
		std::cout << "\tCalculation for MAST" << std::endl;
	}else if(device=='i'){
		gridx = 451;
		gridz = 951;
		gridxmin = 4.0;
		gridzmin = -4.7;
		rmin = 4.1;
		rmax = 8.25;
		std::cout << "\tCalculation for ITER" << std::endl;
	}else if(device=='d'){
		gridx = 180;
		gridz = 400;
		gridxmin = 0.2;
		gridzmin = -2.0;
		rmin = 0.2;
		rmax = 2.0;
		std::cout << "\tCalculation for Double Null MAST 17839files shot" << std::endl;
	}
	else std::cout << "Invalid tokamak" << std::endl;
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

	//impurity.open("output///impurity///impurity.vtk");
	readdata();	// Read data
	// Everything should now be initialised correctly...
}

// *************************************** READING FUNCTIONS *************************************** //

// Read an input file
void PlasmaGrid::readscalars(std::ifstream &input){
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
			mevap[i][k] = 0.0;
		}
	}
}

void PlasmaGrid::readthreevectors(std::ifstream &input){
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
		}
	}
}

void PlasmaGrid::readgridflag(std::ifstream &input){
	int i,k;
	double dummy1,dummy2;
	for(k=0;k<=gridz-1;k++){
		for(i=0;i<=gridx-1;i++){
			input >> dummy1 >> dummy2 >> gridflag[i][k];
		}
	}
}

void PlasmaGrid::readdata(){
	// Input files
	std::ifstream scalars,threevectors,gridflagfile;
	if(device=='m'){
		scalars.open("Models/PlasmaData/b2processed.dat");
		threevectors.open("Models/PlasmaData/b2processed2.dat");
		gridflagfile.open("Models/PlasmaData/locate.dat");
	}
	else if(device=='i')
	{
		scalars.open("Models/PlasmaData/regulardataITER/b2processed.dat");
		threevectors.open("Models/PlasmaData/regulardataITER/b2processed2.dat");
		gridflagfile.open("Models/PlasmaData/regulardataITER/locate.dat");
	}
	else if(device=='d')
	{
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

// *************************************** MUTATION FUNCTIONS *************************************** //

// Locate dust particle in the plasma grid
void PlasmaGrid::locate(int &i, int &k, const threevector xd)const{
	std::cout << "\tIn PlasmaGrid::locate(int &" << i << ", int &" << k << ", " << xd << ")\n\n";
	// Adding 0.5 makes the rounding work properly
	i = int(0.5+(xd.getx()-gridxmin)/dl);
	k = int(0.5+(xd.getz()-gridzmin)/dl);
	checkingrid(i,k);
//	std::cout << "\ni = " << i << "\nk = " << k;
}

void PlasmaGrid::checkingrid(const int i, const int k)const{
//	std::cout << "\tIn PlasmaGrid::checkingrid(int " << i << ", int " << k  << ")\n\n";
	assert( i < gridx && i > 0);
	assert( k < gridz && k > 0);
}

// *************************************** PRINTING FUNCTIONS *************************************** //

// Print the inside and the outside of the tokamak
void PlasmaGrid::vtkcircle(double r, std::ofstream &fout){
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

void PlasmaGrid::vtktoroid()
{
	std::ofstream inner("output/innerplasma.vtk"),outer("output/outerplasma.vtk");
	vtkcircle(rmin,inner);
	vtkcircle(rmax,outer);
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
void PlasmaGrid::datadump(double totalmass){
	vtktoroid();
	//impurityprint(totalmass);
}
