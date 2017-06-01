//#define PAUSE
#define MODEL_DEBUG

#include "Model.h"

Model::Model():Sample(new Tungsten),Pdata(PlasmaDefaults),Accuracy(1.0){
	Mo_Debug("\n\nIn Model::Model():Sample(new Tungsten),Pdata(PlasmaDefaults),Accuracy(1.0)");
}

Model::Model( Matter *&sample, PlasmaData const& pdata, double accuracy )
		:Sample(sample),Pdata(pdata),Accuracy(accuracy){
	Mo_Debug("\n\nIn Model::Model( Matter *& sample, PlasmaData const& pdata ):Sample(sample),Pdata(pdata),Accuracy(accuracy)");
//	std::cout << "\nAccuracy = " << Accuracy;
//	std::cout << "\naccuracy = " << accuracy;
}

/*void Model::Reset_Data( std::shared_ptr <Matter> const& sample, PlasmaData const& pdata ){
	Sample = sample;
	Pdata = pdata;
}*/
