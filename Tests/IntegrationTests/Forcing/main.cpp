//## Forcing:
// This contains tests which aim to establish the position of the dust grain after a finite
// number of time steps. Comparisons have been made to the expected analytical result.
// Those without analytical result have not been considered.
// The dust is started with the following parameters:
// Size = 5e-8;                            // m
// Potential = 1;                          // Normalised Potential
// Temp = 280;                             // K
// InitialPosition = (0.0,0.0,0.0)         // m    
// The initial velocity is different depending on the test being conducted.
// A maximum number of time steps is set which is different for each term being tested.

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

// 		Forcing test 1, GravityTest :
		// Test if constant acceleration under gravity matches the analytical result
		// Produce expected results over 85 or so steps. The equation being solved is:
		// Velocity = InitialVelocity + gravity*Time
//		std::cout << "\nRunning Test 1...\n";
//		out = GravityTest(Element[i]);	
//		std::cout << "\nFinished Running Test 1!\n";


// 		Forcing test 2, ConstantForceTest :
		// Test if constant Lorentz force (Electric field with no magnetic field) and gravity
		// Produce expected results over 100 or so steps. The equation being solved is:
		// Velocity = InitialVelocity + (Electricfield*chargetomassratio + gravity)*Time
//		std::cout << "\nRunning Test 2...\n";
//		out = ConstantForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 2!\n";

// 		Forcing test 3, ConstantMagneticForceTest :
		// Test if model of constant Magnetic field can be replicated analytically
		// Produce expected results over 1000000 or so steps. The equation being solved is:
		// mass*dV/dt = (Velocity X MagneticField), Force equation!
//		std::cout << "\nRunning Test 3...\n";
//		out = ConstantMagneticForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 3!\n";


		// Test bout 4, DON'T THINK THIS TEST WORKS YET, don't know why...
// 		Forcing test 4, ConstantLorentzForceTest :
		// Test if full Lorentz force (Electric field and magnetic field)
		// Produce expected results over 1000000 or so steps. The Equation being solved is:
		// mass*dV/dt = (charge*ElectricField + velocity X MagneticField), Force equation!
//		std::cout << "\nRunning Test 4...\n";
//		out = ConstantLorentzForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 4!\n";

		// Test bout 5, THIS TEST DOESN'T WORK YET, Don't know why...
// 		Forcing test 5, ConstantLorentzGravityForceTest :
		// Test if full Lorentz force (Electric field and magnetic field) and gravity
		// Produce expected results over 1000000 or so steps. The equation being solved is:
		// mass*dV/dt = (charge*ElectricField + velocity X MagneticField) + mass*gravity, Force equation!
//		std::cout << "\nRunning Test 5...\n";
//		out = ConstantLorentzGravityForceTest(Element[i]);	
//		std::cout << "\nFinished Running Test 5!\n";

// 		Forcing test 6, NeutralDragTest :
		// Test if Neutral Drag force 
		// Produce expected results over 100 or so steps. The equation being solved is:
		// Velocity = InitialVelocity + NeutralDragAcceleration*Time;
//		std::cout << "\nRunning Test 6...\n";
//		out = NeutralDragTest(Element[i]);	
//		std::cout << "\nFinished Running Test 6!\n";

// 		Forcing test 7, IonAndNeutralDragTest :
		// Test if Neutral and Ion Drag forces
		// Produce expected results over 100 or so steps. The equation being solved is:
                // Velocity = InitialVelocity + (NeutralDragAcceleration + IonDragAcceleration)*Time;
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
