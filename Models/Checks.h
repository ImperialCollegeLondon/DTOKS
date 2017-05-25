#ifndef __CHECKS_H_INCLUDED__   // if Matter.h hasn't been included yet...
#define __CHECKS_H_INCLUDED__

//#include <string>
//#include <iostream>

void CheckPos(double Value, std::string ErrorMssg){
	if( Value <= 0 ){
		std::cerr << "Error! " << ErrorMssg << " is Zero or negative. " << ErrorMssg << " = " << Value << "\n\n";
		assert(Value > 0);
	}
}

#endif 
