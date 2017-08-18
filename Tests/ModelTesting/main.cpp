#include <ctime>			// Time program
#include "NeutralDragTest.h"
#include "NeutralHeatingTest.h"
#include "EvaporativeCoolingTest.h"
#include "VariableHeatCapacityTest.h"
#include "VariableEmissivityTest.h"

int main(){

	char Element[4]={'W','B','F','G'};// Element, (W) : Tungsten, (G) : Graphite, (B) : Beryllium or (F) : Iron.
	int out(0);
	for(int i(0); i < 1; i ++){

		std::cout << "\nElement["<<i<<"] is = ";
		std::cout << Element[i];

//		Error Estimate Test 1:
//		out = NeutralDragTest(Element[i],false);
//		out = NeutralDragTest(Element[i],true);

//		Error Estimate Test 2:
//		out = NeutralHeatingTest(Element[i],false);
//		out = NeutralHeatingTest(Element[i],true);

//		Error Estimate Test 3:
//		out = EvaporativeCoolingTest(Element[i],false);
//		out = EvaporativeCoolingTest(Element[i],true);

//		Error Estimate Test 4:
		out = VariableHeatCapacityTest(Element[i],false);
		out = VariableHeatCapacityTest(Element[i],true);

//		Error Estimate Test 5:
//		out = VariableEmissivityTest(Element[i],false);
//		out = VariableEmissivityTest(Element[i],true);

		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";
	
		std::cout << "\n\n*****\n";	
	}
	return 0;
}
