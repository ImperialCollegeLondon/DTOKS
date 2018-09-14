//#include "DTOKSU.h"
#include "DTOKSU_Manager.h"

int main(int argc, char* argv[]){

	std::cout << "\n * INITIALISING DTOKS * \n";
	DTOKSU_Manager MyManager;

//	int Config_Status = MyManager.Configure(argc,argv);
	int Config_Status = MyManager.Configure(argc,argv,"Config_Files/DTOKSU_Config.cfg");

	return MyManager.Breakup();
//	return 0;
}
