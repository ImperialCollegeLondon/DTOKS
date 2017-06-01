#ifndef __PLASMAGRID_H_INCLUDED__   // if plasmagrid.h hasn't been included yet...
#define __PLASMAGRID_H_INCLUDED__

#include "PlasmaData.h"
#include "threevector.h"
#include "Constants.h"
#include <iostream>

class plasmagrid{
	// Plasma variables
	protected:
		double **Te,**Ti,**na0,**na1,**po,**ua0,**ua1,**bx,**by,**bz,**x,**z,**mevap;
		threevector E,B,vp,g;
		int **gridflag;
		double mi,gamma,gridxmin,gridzmin,rmin,rmax,dl;
		int gridx,gridz;  

		// Material type
		char gas, device;

		// Files
//		std::ofstream impurity;
	public:
		plasmagrid(char element, char machine, double spacing);	// Constructor: input the plasma gas and the tokamak with the grid spacing
		~plasmagrid();		// Destructor - frees up the memory

		// Methods to get variables
		PlasmaData get_plasmadata(threevector pos);
		double getTe		(int i, int k)	{return Te[i][k]/echarge;}
		double getTi		(int i, int k)	{return Ti[i][k]/echarge;}
		double getna0		(int i, int k)	{return na0[i][k];}
		double getna1		(int i, int k)	{return na1[i][k];}
		double getna1mi		(int i, int k)	{return na1[i][k]*mi;}
		double getmevap		(int i, int k)	{return mevap[i][k];}
		double getmi		()		{return mi;}
		double getgamma		()		{return gamma;}
		double getdl		()		{return dl;}
		threevector getE	()		{return E;}
		threevector getB	()		{return B;}
		threevector getvp	()		{return vp;}
		threevector getg	()		{return g;}
		void summevap(int i, int k, double dm){mevap[i][k] += dm;}

		void vtkcircle(double r, std::ofstream &fout); // Print the inside and the outside of the tokamak 
		void vtktoroid();
		void impurityprint(double totalmass);
		void datadump(double totalmass);
		void readscalars(std::ifstream &input); // Read an input file
		void readthreevectors(std::ifstream &input);
		void readgridflag(std::ifstream &input);
		void readdata();
		void setfields(int i, int k);
		void locate(int &i, int &k, threevector xd); // Locate dust particle in the plasma grid		
		bool withingrid(int i, int k);
};

#endif
