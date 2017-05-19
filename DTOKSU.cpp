//#define PAUSE
//#define DTOKSU_DEBUG
#include "DTOKSU.h"

// Default Constructor, no arguments
DTOKSU::DTOKSU(){
}


void DTOKSU::CreateFile(std::string filename){

	MyFile.open(filename);
	MyFile << "Yo";
	MyFile << "\n";
}

void DTOKSU::CheckTimeStep(){
	D_Debug( "\n\nIn DTOKSU::CheckTimeStep()" );
	assert(TimeStep > 0);
	D_Debug( "\nTimeStep = " << TimeStep );
	TotalTime += TimeStep;
}

void DTOKSU::Print(){
	D_Debug("\n\nIn DTOKSU::Print()");

	MyFile 	<< "\nSomeStuff";

	MyFile << "\n";
}

int DTOKSU::Run(){
	return 0;
}
// ************************************* \\
