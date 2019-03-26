//#include "DTOKSU.h"


#include "DTOKSU_Manager.h"
//#include "CurrentTerms.h"
//#include "ForceTerms.h"
//#include "HeatTerms.h"

int main(int argc, char* argv[]){

/*
    std::vector<ForceTerm*> My_ForceTerms;
    std::vector<HeatingTerm*> My_HeatTerms;
    std::vector<ChagingTerm*> My_ChargeTerms;
    My_ForceTerms.push_back(new Term::Gravity());
    My_ForceTerms.push_back(new Term::Centrifugal());
    My_ForceTerms.push_back(new Term::LorentzForce());
    My_ForceTerms.push_back(new Term::SOMLIonDrag());
    My_ForceTerms.push_back(new Term::SMOMLIonDrag());
    My_ForceTerms.push_back(new Term::HybridIonDrag());
    My_ForceTerms.push_back(new Term::DTOKSIonDrag());
    My_ForceTerms.push_back(new Term::DUSTTIonDrag());
    My_ForceTerms.push_back(new Term::NeutralDrag());
    My_ForceTerms.push_back(new Term::RocketForce());

    My_HeatTerms.push_back(new Term::EmissivityModel());
    My_HeatTerms.push_back(new Term::EvaporationModel());
    My_HeatTerms.push_back(new Term::NewtonCooling());
    My_HeatTerms.push_back(new Term::NeutralHeatFlux());
    My_HeatTerms.push_back(new Term::SOMLIonHeatFlux());
    My_HeatTerms.push_back(new Term::SOMLNeutralRecombination());
    My_HeatTerms.push_back(new Term::SMOMLIonHeatFlux());
    My_HeatTerms.push_back(new Term::SMOMLNeutralRecombination());
    My_HeatTerms.push_back(new Term::SEE());
    My_HeatTerms.push_back(new Term::TEE());
    My_HeatTerms.push_back(new Term::PHLElectronHeatFlux());
    My_HeatTerms.push_back(new Term::OMLElectronHeatFlux());
    My_HeatTerms.push_back(new Term::DTOKSSEE());
    My_HeatTerms.push_back(new Term::DTOKSTEE());
    My_HeatTerms.push_back(new Term::DTOKSIonHeatFlux());
    My_HeatTerms.push_back(new Term::DTOKSNeutralRecombination());
    My_HeatTerms.push_back(new Term::DTOKSElectronHeatFlux());
    My_HeatTerms.push_back(new Term::DUSTTIonHeatFlux());

    My_ChargeTerms.push_back(new Term::solveOML());
    My_ChargeTerms.push_back(new Term::solvePosOML());
    My_ChargeTerms.push_back(new Term::solvePHL());
    My_ChargeTerms.push_back(new Term::solveCW());
    My_ChargeTerms.push_back(new Term::solveTHS());
    My_ChargeTerms.push_back(new Term::solveMOML());
    My_ChargeTerms.push_back(new Term::solveSOML());
    My_ChargeTerms.push_back(new Term::solveMOMLWEM());
    My_ChargeTerms.push_back(new Term::solveMOMLEM());

    threevector Acceleration(0.0,0.0,0.0);
    Matter* Sample = new Tungsten(1e-6,300);

    for( int i(0); i <My_ForceTerms.size(); i ++ ){
        Acceleration = Acceleration + My_ForceTerms[i]->Evaluate(*&Sample,&PlasmaDataDefaults);
        std::cout << (My_ForceTerms[i]->Evaluate(*&Sample,&PlasmaDataDefaults)) << "\t";
    }

    std::cout << Acceleration << "\n\n";
    
    double Power(0.0);
    for( int i(0); i <My_HeatTerms.size(); i ++ ){
        Power = Power + My_HeatTerms[i]->Evaluate(*&Sample,&PlasmaDataDefaults,Sample->get_temperature());
        std::cout << (My_HeatTerms[i]->Evaluate(*&Sample,&PlasmaDataDefaults,Sample->get_temperature())) << "\t";
    }

    std::cout << Power << "\n\n";
    
    double Potential(0.0);
    for( int i(0); i <My_ChargeTerms.size(); i ++ ){
        Potential = Potential + My_ChargeTerms[i]->Evaluate(*&Sample,&PlasmaDataDefaults,Potential);
        std::cout << (My_ChargeTerms[i]->Evaluate(*&Sample,&PlasmaDataDefaults,Potential)) << "\t";
    }

    std::cout << Power << "\n\n";*/

    std::cout << "\n * INITIALISING DTOKS * \n";
    DTOKSU_Manager MyManager;

//	int Config_Status = MyManager.Configure(argc,argv);
	try{
		int Config_Status = MyManager.Configure(argc,argv,"Config_Files/DTOKSU_Config.cfg");
		return MyManager.Run();
	}catch(std::exception &e){
		std::cout << "\nException Caught!\n";
		std::cout << e.what() << "\n";
	}
    return 0;
}
