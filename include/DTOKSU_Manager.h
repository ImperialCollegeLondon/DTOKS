#ifndef __DTOKSU_Manager_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __DTOKSU_Manager_H_INCLUDED__

#include "DTOKSU.h"
#include <random>	// for std::normal_distribution<> etc.
#include <chrono>	// for chrono::high_resolution_clock::now().time_since_epoch().count();
#include <config4cpp/Configuration.h>	// For Configuration
#include <netcdfcpp.h>	// for reading netcdf files

class DTOKSU_Manager{

	private:
		// Private member data
		DTOKSU *Sim;
		Matter *Sample;
		PlasmaGrid_Data Pgrid;
		PlasmaData Pdata;

		int Config_Status;

		template<typename T> int input_function(int &argc, char* argv[], int &i, std::stringstream &ss0, T &Temp)const;
		void show_usage(std::string name)const;
		void config_message()const;

		int configure_plasmagrid(std::string plasma_dirname);
		int read_data(std::string plasma_dirname);
		int read_MPSIdata(std::string plasma_dirname);
	
	public:
		DTOKSU_Manager();
		DTOKSU_Manager( int argc, char* argv[] );
		DTOKSU_Manager( int argc, char* argv[], std::string = "Config_Files/DTOKSU_Config.cfg" );

		~DTOKSU_Manager(){
		};
		
		int Configure(std::string = "Config_Files/DTOKSU_Config.cfg");
		int Configure(int argc, char* argv[], std::string = "Config_Files/DTOKSU_Config.cfg");


		int Breakup();
		int Run();
};

#endif
