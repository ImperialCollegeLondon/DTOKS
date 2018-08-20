//#include "DTOKSU.h"
#include "DTOKSU_Manager.h"

int main(int argc, char* argv[]){

	std::cout << "\n * INITIALISING DTOKS * \n";
	DTOKSU_Manager MyManager;

	int Config_Status = MyManager.Configure(argc,argv);

	return MyManager.Breakup();
//	return 0;
}
