#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include "Constants.h"

#define xmin 0.2
#define xmax 1.4
#define xstep 0.01
#define ymin -2.0
#define ymax 0.8
#define ystep 0.01

// Code to generate a .dat file of plasma background for the DTOKSU.cpp code.
// Specifies: Te: Electron Temperature, Ti: Ion Temperature, na0: Ion Density, na1: Electron Density
// po: Plasma Pressure, ua0: Ion Velocity, ua1: Electron Velocity.

double LinearRamp(double fx_init, double fx_fin, double fy_init, double fy_fin, const double &x, const double &y){
	double y_gradient = (fy_fin-fy_init)/(ymax-ymin);
	double x_gradient = (fx_fin-fx_init)/(xmax-xmin);

	double y_part = y_gradient*(y-ymin) + fy_init;
	double x_part = x_gradient*(x-xmin) + fx_init;

	return (y_part + x_part);
}

void MagneticField(std::string magfieldtype, double Bmax){
	std::string filename("BField.dat");
	std::ofstream fs;

	fs.open(filename);

	fs << "x\ty\tbbx\tbby\tbbz\n";

	double bbx(0.0), bby(0.0), bbz(0.0);
	double bbx_i(0.0), bby_i(0.0), bbz_i(0.0);
	double bbx_f(0.0), bby_f(0.0), bbz_f(0.0);

	if	( magfieldtype == "Constant"	){
		bbx = 0.0;
		bby = 0.0;
		bbz = Bmax;
	}else if( magfieldtype == "Linear"	){
		bbx_i  = 0.0; 
		bbx_f  = Bmax; 
		bby_i  = 0.0; 
		bby_f  = Bmax; 
		bbz_i  = 0.0; 
		bbz_f  = Bmax; 

		if( true ) {
			bbx_f = bbx_i;
			bby_f = bby_i;
//			bbz_f = bbz_i;
		}

	}else if( magfieldtype == "Blank"	){
	}else{
		std::cout << "\nInvalid Magnetic Field Type!\n";
 		return;
	}
	for(double y(ymin); y <= ymax; y += ystep){
		for(double x(xmin); x <= xmax; x += xstep){
			if( magfieldtype == "Linear" ){
				bbx = LinearRamp(bbx_i,bbx_f,bbx_i,bbx_f,x,y);
				bby = LinearRamp(bby_i,bby_f,bby_i,bby_f,x,y);
				bbz = LinearRamp(bbz_i,bbz_f,bbz_i,bbz_f,x,y);
			}
			fs << x << "\t" << y  << "\t" << bbx << "\t" << bby << "\t" << bbz << "\n";
		}
	}

	fs.close();
}

int main(){

	double echarge(1.60217662e-19);
	// Te: Electron Temperature, Ti: Ion Temperature, na0: Ion Density, na1: Electron Density
	// po: Plasma Pressure, ua0: Ion Velocity, ua1: Electron Velocity
	std::string PlasmaType = "Constant";	// Plasma Type, can be "Blank", "Constant", "Linear" or "TwoPointModel"

	double Te(0.0), Ti(0.0), na0(0.0), na1(0.0), po(0.0), ua0(0.0), ua1(0.0);

	double Te_i(0.0), Te_f(0.0), Ti_i(0.0), Ti_f(0.0), na0_i(0.0), na0_f(0.0);
	double na1_i(0.0), na1_f(0.0), po_i(0.0), po_f(0.0), ua0_i(0.0), ua0_f(0.0), ua1_i(0.0), ua1_f(0.0);

	double gamma = 7;		// 
	double qParallel = 10e6; 	// Wm^-2, Power flux
	double L = ymax-ymin;		// Length of System
	double k0e = 2000;		// Heat conductivity of electrons

	double Densityu = 1.0e17;	
	double TwoPointModelTempTt = 4*Mp*pow(qParallel,2)*pow(7*qParallel*L/(2*k0e),-(4.0/7.0))/(2*pow(echarge,3)*pow(gamma*Densityu,2));
	double TwoPointModelTempTu = pow(7*qParallel*L/(2*k0e),-(4.0/7.0));
	double Densityt = pow(Densityu*echarge,2)*gamma*pow(7*qParallel*L/(2*k0e),-(4.0/7.0))/(2*qParallel*Mp);


	if	( PlasmaType == "Constant"	){
		Te  = 50*echarge;//1.5*echarge;	// Convert ev to joules
		Ti  = 50*echarge;//1.188*echarge;	// Convert ev to joules
		na0 = 1.0e18;//3.8e18; 
		na1 = 1.0e18;
		po  = 1.0;		// Varies between -20 and 20
		ua0 = 10000.0; 		// Varies between -10,000 and 10,000 (roughly), generally smaller than ua1
		ua1 = 10000.0; 		// Varies between -10,000 and 10,000 (roughly)

	}else if( PlasmaType == "Linear"	){
		Te_i  = 1*echarge;
		Te_f  = 100*echarge;
		Ti_i  = 1*echarge;
		Ti_f  = 100*echarge;
		na0_i = 1.0e16;
		na0_f = 1.0e20;
		na1_i = 1.0e16;
		na1_f = 1.0e20;
		po_i = 0.0;
		po_f = 0.0;
		ua0_i = 0.0;
		ua0_f = 0.0;
		ua1_i = 0.0;
		ua1_f = 0.0;

		if( true ) { // Make Constant
			Te_f  = Te_i;
//			Ti_f  = Ti_i;
			na0_f = na0_i;
			na1_f = na1_i;
			po_f  = po_i;
			ua1_f = ua1_i;
			ua1_f = ua1_i;
		}

	}else if( PlasmaType == "TwoPointModel"	){ // Constant part of two point model
		na0 = Densityu;		// 3.8e18; 
		na1 = Densityu;
		Te  = Kb*TwoPointModelTempTt;	// Convert ev to joules
		Ti  = Kb*TwoPointModelTempTt;	// Convert ev to joules
		po  = 0.0;		// Varies between -20 and 20
		ua0 = 0.0; 		// Varies between -10,000 and 10,000 (roughly), generally smaller than ua1
		ua1 = 0.0; 		// Varies between -10,000 and 10,000 (roughly)
	}else if( PlasmaType == "Blank" 	){
	}else{
		std::cout << "\nInvalid Plasma Type!\n";
 		return 1;
	}


	std::string filename((PlasmaType + ".dat").c_str());
	std::ofstream fs;

	fs.open(filename);

	fs << "x\ty\tte\tti\tna0\tna1\tpo\tua0\tua1\n";

	for(double y(ymin); y < (ymax+ystep); y += ystep){
		for(double x(xmin); x < (xmax+xstep); x += xstep){
			if( PlasmaType == "Linear" ){	// Note, at the moment it's got the y component as constant

				Te  = LinearRamp(Te_i,Te_i,Te_i,Te_f,x,y);
				Ti  = LinearRamp(Ti_i,Ti_i,Ti_i,Ti_f,x,y);
				na0 = LinearRamp(na0_i,na0_i,na0_i,na0_f,x,y);
				na1 = LinearRamp(na1_i,na1_i,na1_i,na1_f,x,y);
				po  = LinearRamp(po_i,po_i,po_i,po_f,x,y);
				ua0 = LinearRamp(ua0_i,ua0_i,ua0_i,ua0_f,x,y);
				ua1 = LinearRamp(ua1_i,ua1_i,ua1_i,ua1_f,x,y);
			}else if( PlasmaType == "TwoPointModel" ){
				if( y < -1.6 ){
					// LINEAR RAMPS ARE WRONG!
					// THE DENOMINATOR IS THE ENTIRE Y RANGE!
					Te  = Kb*LinearRamp(0.0,0.0,TwoPointModelTempTt,TwoPointModelTempTu,x,y);
					Ti  = Kb*LinearRamp(0.0,0.0,TwoPointModelTempTt,TwoPointModelTempTu,x,y);
					na0 = LinearRamp(0.0,0.0,Densityt,Densityu,x,y);
					na1 = LinearRamp(0.0,0.0,Densityt,Densityu,x,y);
					po  = 0; 
					ua0 = sqrt(Ti/Mp);
					ua1 = 0;
				}else{
			                na0 = Densityu;         // 3.8e18; 
			                na1 = Densityu;
			                Te  = Kb*TwoPointModelTempTt;   // Convert ev to joules
			                Ti  = Kb*TwoPointModelTempTt;   // Convert ev to joules
				}
			}

			ua0 = sqrt(Ti/Mp);
			ua1 = sqrt(Ti/Mp);
			fs << x << "\t" << y  << "\t" << Te << "\t" << Ti << "\t" << na0 << "\t" 
				<< na1 << "\t" << po << "\t" << ua0 << "\t" << ua1 << "\n";
		}
	}

	fs.close();

	std::cout << "\nFinished Generating Plasma\n!";

	MagneticField("Constant",1.0);
	std::cout << "\nFinished Generating Magnetic Field!\n";
	return 0;
}
