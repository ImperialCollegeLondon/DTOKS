#ifndef __PLASMAGRID_H_INCLUDED__   // if plasmagrid.h hasn't been included yet...
#define __PLASMAGRID_H_INCLUDED__

#include <assert.h>
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
		PlasmaData get_plasmadata(const threevector pos)const;
		double getTe		(int i, int k)	const{
			checkingrid(i,k);
			if(Te[i][k] == Te[i][k] ){
				return Te[i][k]/echarge;
			}else{
				return 0;
			}
		}
		double getTi		(int i, int k)	const{
			checkingrid(i,k);
			if(Ti[i][k] == Ti[i][k] ){
				return Ti[i][k]/echarge;
			}else{
				return 0;
			}
		}
		double getna0		(int i, int k)	const{
			checkingrid(i,k);
                        if(na0[i][k] == na0[i][k] ){
                                return na0[i][k];
                        }else{
                                return 0;
                        }
		}
		double getna1		(int i, int k)	const{
			checkingrid(i,k);
                        if(na1[i][k] == na1[i][k] ){
                                return na1[i][k];
                        }else{
                                return 0;
                        }
		}
		double getna1mi		(int i, int k)	const{
			checkingrid(i,k);
                        if(na1[i][k]*mi == na1[i][k]*mi ){
                                return na1[i][k]*mi;
                        }else{
                                return 0;
                        }
		}
		double getmevap		(int i, int k)	const{
			checkingrid(i,k);
                        if(mevap[i][k] == mevap[i][k] ){
                                return mevap[i][k];
                        }else{
                                return 0;
                        }
		}
		double getmi		()		const{return mi;}
		double getgamma		()		const{return gamma;}
		double getdl		()		const{return dl;}
		threevector getE	()		const{return E;}
		threevector getB	()		const{return B;}
		threevector getvp	()		const{return vp;}
		threevector getg	()		const{return g;}
		void summevap(int i, int k, double dm){mevap[i][k] += dm;}
		void checkingrid(int i, int k)const;

		void vtkcircle(double r, std::ofstream &fout); // Print the inside and the outside of the tokamak 
		void vtktoroid();
		void impurityprint(double totalmass);
		void datadump(double totalmass);
		void readscalars(std::ifstream &input); // Read an input file
		void readthreevectors(std::ifstream &input);
		void readgridflag(std::ifstream &input);
		void readdata();
		void setfields(int i, int k);
		void locate(int &i, int &k, threevector xd)const; // Locate dust particle in the plasma grid		
		bool withingrid(int i, int k);
};

#endif
