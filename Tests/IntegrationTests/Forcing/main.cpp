#include <ctime>			// Time program
#include "GravityTest.h"
#include "ConstantForceTest.h"
#include "ConstantMagneticForceTest.h"
#include "ConstantLorentzForceTest.h"
#include "ConstantLorentzGravityForceTest.h"
#include "NeutralDragTest.h"
#include "IonAndNeutralDragTest.h"

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
//		std::cout << "\nRunning Test 1...\n";
//		out = GravityTest(Element[i]);	
//		std::cout << "\nFinished Running Test 1!\n";

		// Test bout 2,
		// Integration test two for all materials
		// Test if constant Lorentz force (Electric field with no magnetic field) and gravity
		// Produce expected results over 100 or so steps
//		std::cout << "\nRunning Test 2...\n";
//		out = ConstantForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 2!\n";

		// Test bout 3,
		// Integration test three for all materials
		// Test if model of constant Magnetic field can be replicated analytically
		// Produce expected results over 1000000 or so steps
//		std::cout << "\nRunning Test 3...\n";
//		out = ConstantMagneticForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 3!\n";


		// Test bout 4,
		// Integration test four for all materials
		// Test if full Lorentz force (Electric field and magnetic field)
		// Produce expected results over 1000000 or so steps
//		std::cout << "\nRunning Test 4...\n";
//		out = ConstantLorentzForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 4!\n";

		// Test bout 5, THIS TEST DOESN'T WORK YET, Don't know why...
		// Integration test five for all materials
		// Test if full Lorentz force (Electric field and magnetic field) and gravity
		// Produce expected results over 1000000 or so steps
//		std::cout << "\nRunning Test 5...\n";
//		out = ConstantLorentzGravityForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 5!\n";

		// Test bout 6, 
		// Integration test six for all materials
		// Test if Neutral Drag force 
		// Produce expected results over 100 or so steps
//		std::cout << "\nRunning Test 6...\n";
//		out = NeutralDragTest(Element[i]);	
//		std::cout << "\nFinished Running Test 6!\n";

		// Test bout 7, 
		// Integration test seven for all materials
		// Test if Neutral and Ion Drag forces
		// Produce expected results over 100 or so steps
		std::cout << "\nRunning Test 7...\n";
		out = IonAndNeutralDragTest(Element[i]);	
		std::cout << "\nFinished Running Test 7!\n";


		if( out == 1 ) std::cout << "\n# PASSED!";
		if( out == 2 ) std::cout << "\n# WARNING! DEVIATION OF < 1%";
		if( out == -1 ) std::cout << "\n# FAILED!";

		std::cout << "\n\n*****\n\n\n";	
	}
	return 0;
}
