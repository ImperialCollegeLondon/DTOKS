#include <ctime>			// Time program
#include "TemperatureEvolutionEquation.h"
#include "NonConstHeatCapacity.h"

int main(){

	char Element[4]={'W','B','F','G'};// Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
	for(int i(0); i < 4; i ++){

		std::cout << "\nElement["<<i<<"] is = ";
		std::cout << Element[i];
//		int out = TemperatureEvolutionEquation(Element[i]);
		int out = NonConstHeatCapacity(Element[i]);
	
		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";
	
		std::cout << "\n\n*****\n";	
	}
	return 0;
}
