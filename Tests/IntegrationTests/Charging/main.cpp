//## Charging
//This contains tests which aim to establish the normalised potential of the dust grain at
//equilibrium. Since no analytic solution exists for the charge on a dust grain (ignoring
//the Lambert W function), the numerical solution is compared to a pre-calculated result

#include <ctime>			// Time program
#include "ChargeTest.h"

int main(){

	char Element[4]={'W','B','F','G'}; // Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
	int out(0);
	for(int i(0); i < 4; i ++){
		std::cout << "\nElement["<<i<<"] is = ";
		std::cout << Element[i];

// 		Charging test 1, ChargeTest :
//		This test simply compares the DTOKS OML result as calculated by the code with the same
//		expression coded identically to verify the results match
		std::cout << "\nRunning Test 1...\n";
		out = ChargeTest(Element[i]);	
		std::cout << "\nFinished Running Test 1!\n";

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";

		std::cout << "\n\n*****\n\n\n";	
	}
	return 0;
}
