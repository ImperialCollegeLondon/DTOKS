#include "plasmagrid.h"

plasmagrid::plasmagrid(char element, char machine, double spacing)
{
	int i;
	// Plasma parameters
	mi=Mp;
	gamma=Me/mi;
	gas = element;
	device = machine;
	g = threevector(0.0,0.0,-9.8);
	if(device=='m')
	{
		gridx = 121;
		gridz = 281;
		gridxmin = 0.2;
		gridzmin = -2.0;
		rmin = 0.2;
		rmax = 1.4;
		std::cout << "Calculation for MAST" << std::endl;
	}
	else if(device=='i')
	{
		gridx = 451;
		gridz = 951;
		gridxmin = 4.0;
		gridzmin = -4.7;
		rmin = 4.1;
		rmax = 8.25;
		std::cout << "Calculation for ITER" << std::endl;
	}
	else if(device=='d')
	{
		gridx = 180;
		gridz = 400;
		gridxmin = 0.2;
		gridzmin = -2.0;
		rmin = 0.2;
		rmax = 2.0;
		std::cout << "Calculation for Double Null MAST 17839files shot" << std::endl;
	}
	else std::cout << "Invalid tokamak" << std::endl;
	dl = spacing;
	Te = new double*[gridx];
	Ti = new double*[gridx];
	na0 = new double*[gridx];
	na1 = new double*[gridx];
	po = new double*[gridx];
	ua0 = new double*[gridx];
	ua1 = new double*[gridx];
	bx = new double*[gridx];
	by = new double*[gridx];
	bz = new double*[gridx];
	x = new double*[gridx];
	z = new double*[gridx];
	mevap = new double*[gridx];
	gridflag = new int*[gridx];

	for(i=0;i<gridx;i++)
	{
		Te[i] = new double[gridz];
		Ti[i] = new double[gridz];
		na0[i] = new double[gridz];
		na1[i] = new double[gridz];
		po[i] = new double[gridz];
		ua0[i] = new double[gridz];
		ua1[i] = new double[gridz];
		bx[i] = new double[gridz];
		by[i] = new double[gridz];
		bz[i] = new double[gridz];
		x[i] = new double[gridz];
		z[i] = new double[gridz];
		mevap[i] = new double[gridz];
		gridflag[i] = new int[gridz];
	}
	//impurity.open("output///impurity///impurity.vtk");
}

// Destructor - frees up the memory
plasmagrid::~plasmagrid()
{
	int i;
	for(i=0;i<gridx;i++)
	{
		delete [] Te[i];
		delete [] Ti[i];
		delete [] na0[i];
		delete [] na1[i];
		delete [] po[i];
		delete [] ua0[i];
		delete [] ua1[i];
		delete [] bx[i];
		delete [] by[i];
		delete [] bz[i];
		delete [] x[i];
		delete [] z[i];
		delete [] mevap[i];
		delete [] gridflag[i];
	}
	delete [] Te;
	delete [] Ti;
	delete [] na0;
	delete [] na1;
	delete [] po;
	delete [] ua0;
	delete [] ua1;
	delete [] bx;
	delete [] by;
	delete [] bz;
	delete [] x;
	delete [] z;
	delete [] mevap;
	delete [] gridflag;
	//impurity.close();
}

// Print the inside and the outside of the tokamak
void plasmagrid::vtkcircle(double r, std::ofstream &fout)
{
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

void plasmagrid::vtktoroid()
{
	std::ofstream inner("output/innerplasma.vtk"),outer("output/outerplasma.vtk");
	vtkcircle(rmin,inner);
	vtkcircle(rmax,outer);
	inner.close();
	outer.close();
}
/*
void plasmagrid::impurityprint(double totalmass)
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
void plasmagrid::datadump(double totalmass)
{
	vtktoroid();
	//impurityprint(totalmass);
}

// Read an input file
void plasmagrid::readscalars(std::ifstream &input)
{
	int i,k;
	char dummy;
	// Ignore first line of file
	for(i=0;i<=19;i++)
	{
		input >> dummy;
	}
	for(k=0;k<=gridz-1;k++)
	{
		for(i=0;i<=gridx-1;i++)
		{
			input >> x[i][k] >> z[i][k] >> Te[i][k] >> Ti[i][k] >> na0[i][k] >> na1[i][k] >> po[i][k] >> ua0[i][k] 
				>> ua1[i][k];
			mevap[i][k] = 0.0;
		}
	}
}

void plasmagrid::readthreevectors(std::ifstream &input)
{
	int i,k;
	char dummy;
	double dummy1, dummy2;
	// Ignore first line of file
	for(i=0;i<=19;i++)
	{
		input >> dummy;
	}
	for(k=0;k<=gridz-1;k++)
	{
		for(i=0;i<=gridx-1;i++)
		{
			input >> dummy1 >> dummy2 >> bx[i][k] >> bz[i][k] >> by[i][k];
			bz[i][k] = -bz[i][k];
		}
	}
}

void plasmagrid::readgridflag(std::ifstream &input)
{
	int i,k;
	double dummy1,dummy2;
	for(k=0;k<=gridz-1;k++)
	{
		for(i=0;i<=gridx-1;i++)
		{
			input >> dummy1 >> dummy2 >> gridflag[i][k];
		}
	}
}

void plasmagrid::readdata()
{
	// Input files
	std::ifstream scalars,threevectors,gridflag;
	if(device=='m')
	{
		scalars.open("regulardata/b2processed.dat");
		threevectors.open("regulardata/b2processed2.dat");
		gridflag.open("regulardata/locate.dat");
	}
	else if(device=='i')
	{
		scalars.open("regulardataITER/b2processed.dat");
		threevectors.open("regulardataITER/b2processed2.dat");
		gridflag.open("regulardataITER/locate.dat");
	}
	else if(device=='d')
	{
		scalars.open("MASTDexp/b2processed.dat");
		threevectors.open("MASTDexp/b2processed2.dat");
		gridflag.open("MASTDexp/locate.dat");
	}
	readscalars(scalars);
	readthreevectors(threevectors);
	readgridflag(gridflag);
	scalars.close();
	threevectors.close();
	gridflag.close();
}

void plasmagrid::setfields(int i, int k)
{
	double aveu;
	if( (na0[i][k]>0.0) || (na1[i][k]>0.0) )
	{
		aveu = (na0[i][k]*ua0[i][k]+na1[i][k]*ua1[i][k])/(na0[i][k]+na1[i][k]);
	}
	else aveu = 0.0;

	// Read magnetic field in
	B.setx(bx[i][k]);
	B.sety(by[i][k]);
	B.setz(bz[i][k]);

	// Plasma velocity is parallel to the B field
	vp.setx(aveu*(B.getunit().getx()));
	vp.sety(aveu*(B.getunit().gety()));
	vp.setz(aveu*(B.getunit().getz()));

	if(gridflag[i][k]==1)
	{
		if( (gridflag[i+1][k]==1) && (gridflag[i-1][k]==1) )
		{
			E.setx(-(po[i+1][k]-po[i-1][k])/(2.0*dl));
		}
		else if(gridflag[i+1][k]==1)
		{
			E.setx(-(po[i+1][k]-po[i][k])/dl);
		}
		else if(gridflag[i-1][k]==1)
		{
			E.setx(-(po[i][k]-po[i-1][k])/dl);
		}
		else E.setx(0.0);
		if((gridflag[i][k+1]==1)&&(gridflag[i][k-1]==1))
		{
			E.setz(-(po[i][k+1]-po[i][k-1])/(2.0*dl));
		}
		else if(gridflag[i][k+1]==1)
		{
			E.setz(-(po[i][k+1]-po[i][k])/dl);
		}
		else if(gridflag[i][k-1]==1)
		{
			E.setz(-(po[i][k]-po[i][k-1])/dl);
		}
		else E.sety(0.0);
	}
	else
	{
		E.setx(0.0);
		E.setz(0.0);
	}
}

// Locate dust particle in the plasma grid
void plasmagrid::locate(int &i, int &k, threevector xd)
{
	// Adding 0.5 makes the rounding work properly
	i = int(0.5+(xd.getx()-gridxmin)/dl);
	k = int(0.5+(xd.getz()-gridzmin)/dl);
}
bool plasmagrid::withingrid(int i, int k)
{
	bool result;
	if( (i>0) && (i<gridx-1) && (k>0) && (k<gridz-1) ) result = true;
	else result = false;
	return result;
}

PlasmaData plasmagrid::get_plasmadata(threevector pos){
	PlasmaData ReturnPdata;
	int i(0), k(0);
	locate(i,k,pos);
	ReturnPdata.NeutralDensity 	= getna0(i,k); 
	ReturnPdata.ElectronDensity 	= getna1(i,k);  // ELECTRON DENSITY EQUALS ION DENSITY
	ReturnPdata.IonDensity 		= getna1(i,k);
	ReturnPdata.IonTemp		= getTi(i,k); 
	ReturnPdata.ElectronTemp 	= getTe(i,k); 
	ReturnPdata.NeutralTemp 	= getTi(i,k); // NEUTRAL TEMP EQUAL TO ION TEMP
	ReturnPdata.AmbientTemp 	= 300; // NOTE THIS IS HARD CODED OHMEINGOD
	ReturnPdata.PlasmaVel	 	= vp; 
	ReturnPdata.ElectricField 	= E; 
	ReturnPdata.MagneticField 	= B; 
	return ReturnPdata;
}

