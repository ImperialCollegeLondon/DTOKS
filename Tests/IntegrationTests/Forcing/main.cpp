#include <ctime>			// Time program
#include "ConstantForceTest.h"

int main(){

	char Element[4]={'W','B','F','G'}; // Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
	int out(0);
	for(int i(0); i < 4; i ++){
		std::cout << "\nElement["<<i<<"] is = ";
		std::cout << Element[i];

		// Test bout 1,
		// Integration test one for all materials
		// Test if constant Lorentz force (Electric field with no magnetic field) and gravity
		// Produce expected results over 100 or so steps
		std::cout << "\nRunning Test 1...\n";
		out = ConstantForceTest(Element[i]);	
		std::cout << "\nFinished Running Test 1!\n";

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";

		std::cout << "\n\n*****\n\n\n";	
	}
	return 0;
}
