#include "DTOKSU_Manager.h"

int main(int argc, char* argv[]){

    std::cout << "\n * INITIALISING DTOKS * \n";
    DTOKSU_Manager MyManager;

    int Config_Status = MyManager.Configure(argc,argv);
    try{
//        int Config_Status = MyManager.Configure(argc,my_argv,"Config_Files/DTOKSU_Config_EAST.cfg");
        int run_status = MyManager.Run();
    }catch(std::exception &e){
        std::cout << "\nException Caught!\n";
        std::cout << e.what() << "\n";
    }
    return 0;
}
