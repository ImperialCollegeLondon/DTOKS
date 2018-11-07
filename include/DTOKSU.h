#ifndef __DTOKSU_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __DTOKSU_H_INCLUDED__

#include "HeatingModel.h"
#include "ForceModel.h"
#include "ChargingModel.h"

#include <algorithm>

static struct Boundary_Data BoundaryDefaults = {
	std::vector<std::pair<double,double>> ()
};

class DTOKSU{

	private:
		// Private member data
		double TotalTime;			// Seconds, total time taken to perform simulation

		Matter *Sample;				// Matter sample can be either Tungsten, Beryllium, Graphite or Iron

		HeatingModel HM;			// Heating Model 
		ForceModel FM;				// Force Model 
		ChargingModel CM;			// Charge Model

		Boundary_Data WallBound, CoreBound;

		std::ofstream MyFile;			// Output data file
	
		// Private Functions
		void print();				// Write to output data file
		void create_file(std::string filename);

		void SpecularReflection();
		bool Boundary_Check(bool InOrOut);

	public:
		static const unsigned int MN = 3;	// MODEL NUMBER, the number of physical models in DTOKS

//		DTOKSU();
		DTOKSU( std::array<float,MN> alvls, Matter *& sample, PlasmaData &pdata,
				std::array<bool,HMN> &heatmodels, std::array<bool,FMN> &forcemodels, 
				std::array<bool,CMN> &chargemodels);
		DTOKSU( std::array<float,MN> alvls, Matter *& sample, PlasmaGrid_Data &pgrid,
				std::array<bool,HMN> &heatmodels, std::array<bool,FMN> &forcemodels, 
				std::array<bool,CMN> &chargemodels);
		DTOKSU( std::array<float,MN> alvls, Matter *& sample, PlasmaGrid_Data &pgrid,
				PlasmaData &pdata,	std::array<bool,HMN> &heatmodels, 
				std::array<bool,FMN> &forcemodels, std::array<bool,CMN> &chargemodels);
		DTOKSU( std::array<float,MN> alvls, Matter *& sample, PlasmaGrid_Data &pgrid,
				PlasmaData &pdata, Boundary_Data &wbound, Boundary_Data &cbound,
				std::array<bool,HMN> &heatmodels, std::array<bool,FMN> &forcemodels, 
				std::array<bool,CMN> &chargemodels);

		~DTOKSU(){
		};

		int Run();
		void OpenFiles(std::string filename, unsigned int i);
		void CloseFiles();			// Close all model files
		void ResetModelTime(double HMTime, double FMTime, double CMTime);
	
		double 		get_HMTime()const	{ 	return HM.get_totaltime(); }
		double 		get_FMTime()const	{ 	return FM.get_totaltime(); }
		double 		get_CMTime()const	{ 	return CM.get_totaltime(); }
	
		threevector get_bfielddir()const{	return (FM.get_bfield());	}
};

#endif
