#ifndef __PLASMAGRID_H_INCLUDED__   // if PlasmaGrid.h hasn't been included yet...
#define __PLASMAGRID_H_INCLUDED__

#include "PlasmaData.h"
#include "threevector.h"
#include "Constants.h"

#include <iostream>
#include <assert.h>
#include <vector>

class PlasmaGrid{
	// Plasma variables
	protected:

		// Oh by the way, na0 is the ion density and na1 is the electron density,
		// ua0 is the ion drift velocity and ua1 is the electron drift velocity! po is the potential. You're welcome.
		std::vector<std::vector<double>> Te, Ti, na0, na1, po, ua0, ua1, bx, by, bz, x, z, mevap;
		std::vector<std::vector<int>> gridflag;
		int gridx, gridz;  
		double gridxmin, gridzmin, rmin, rmax;

		// Material type
		const double mi, gamma, dl;
		const char gas, device;

		// Files
//		std::ofstream impurity;
	public:

		// Constructor: input the plasma gas and the tokamak with the grid spacing
		PlasmaGrid(char element, char machine, double spacing);
		~PlasmaGrid(){};		// Destructor - frees up the memory

		// TO ALLOW FOR THIS CLASS TO BE CONST LIKE IT SHOULD BE:
		// YOU COULD CONSIDER MAKING THESE FOUR FUNCTIONS STATIC AND RETURN VECTORS OF EACH THING
		// THESE COULD THEN BE USED IN THE INITIALISER LIST OF THE CONSTRUCTOR.
		// Functions for reading in the data
		void readscalars(std::ifstream &input); // Read an input file
		void readthreevectors(std::ifstream &input);
		void readgridflag(std::ifstream &input);
		void readdata();

		// functions to locate dust, update electromagnetic fields, calculate mass loss and check position
		void setfields(int i, int k);
		void summevap(int i, int k, double dm){mevap[i][k] += dm;}
		bool locate(int &i, int &k, threevector xd)const; // Locate dust particle in the plasma grid		
		bool checkingrid(int i, int k)const;

		// Methods to get variables
		const PlasmaData &get_plasmadata(const threevector pos)const;
		double getTe		(int i, int k)const{
			checkingrid(i,k);
			if(Te[i][k] == Te[i][k] ){
				return Te[i][k]/echarge;
			}else{
				return 0;
			}
		}
		double getTi		(int i, int k)const{
			checkingrid(i,k);
			if(Ti[i][k] == Ti[i][k] ){
				return Ti[i][k]/echarge;
			}else{
				return 0;
			}
		}
		double getna0		(int i, int k)const{
			checkingrid(i,k);
                        if(na0[i][k] == na0[i][k] ){
                                return na0[i][k];
                        }else{
                                return 0;
                        }
		}
		double getpo		(int i, int k)const{
			checkingrid(i,k);
                        if(po[i][k] == po[i][k] ){
                                return po[i][k];
                        }else{
                                return 0;
                        }
		}
		double getua0		(int i, int k)const{
			checkingrid(i,k);
                        if(ua0[i][k] == ua0[i][k] ){
                                return ua0[i][k];
                        }else{
                                return 0;
                        }
		}
		double getua1		(int i, int k)const{
			checkingrid(i,k);
                        if(ua0[i][k] == ua0[i][k] ){
                                return ua0[i][k];
                        }else{
                                return 0;
                        }
		}
		double getbx		(int i, int k)const{
			checkingrid(i,k);
                        if(bx[i][k] == bx[i][k] ){
                                return bx[i][k];
                        }else{
                                return 0;
                        }
		}
		double getby		(int i, int k)const{
			checkingrid(i,k);
                        if(by[i][k] == by[i][k] ){
                                return by[i][k];
                        }else{
                                return 0;
                        }
		}
		double getbz		(int i, int k)const{
			checkingrid(i,k);
                        if(bz[i][k] == bz[i][k] ){
                                return bz[i][k];
                        }else{
                                return 0;
                        }
		}

		double getna1		(int i, int k)const{
			checkingrid(i,k);
                        if(na1[i][k] == na1[i][k] ){
                                return na1[i][k];
                        }else{
                                return 0;
                        }
		}
		double getna1mi		(int i, int k)const{
			checkingrid(i,k);
                        if(na1[i][k]*mi == na1[i][k]*mi ){
                                return na1[i][k]*mi;
                        }else{
                                return 0;
                        }
		}
		double getmevap		(int i, int k)const{
			checkingrid(i,k);
                        if(mevap[i][k] == mevap[i][k] ){
                                return mevap[i][k];
                        }else{
                                return 0;
                        }
		}
		double getgridflag	(int i, int k)const{
			checkingrid(i,k);
                        if(gridflag[i][k] == gridflag[i][k] ){
                                return gridflag[i][k];
                        }else{
                                return 0;
                        }
		}	
		double getmi		()		const{return mi;}
		double getgamma		()		const{return gamma;}
		double getdl		()		const{return dl;}

		// Functions for printing
		void vtkcircle(double r, std::ofstream &fout); // Print the inside and the outside of the tokamak 
		void vtktoroid();
		void impurityprint(double totalmass);
		void datadump(double totalmass);
};

#endif
