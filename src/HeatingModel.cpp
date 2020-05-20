/** @file HeatingModel.cpp
 *  @brief Implementation of class for physics models relevant to dust heating
 *  
 *  Implement the member functions of the Heating class
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 */

#include "HeatingModel.h"
#include "Constants.h"
#include "Functions.h"

HeatingModel::HeatingModel():
Model(){
    H_Debug("\n\nIn HeatingModel::HeatingModel():Model()\n\n");
    Defaults();
    CreateFile("Default_Heating_filename.txt",false);
}

HeatingModel::HeatingModel(std::string filename, float accuracy, 
std::vector<HeatTerm*> heatterms, Matter *& sample, PlasmaData &pdata):
Model(filename,sample,pdata,accuracy){
    H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, "
        << "float accuracy, std::vector<HeatTerm*> heatterms, "
        << "Matter *& sample, PlasmaData const *&pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    Defaults();
    HeatTerms = heatterms;
    CreateFile(filename,false);
}

HeatingModel::HeatingModel(std::string filename, float accuracy, 
std::vector<HeatTerm*> heatterms, Matter *& sample, PlasmaData *pdata):
Model(filename,sample,*pdata,accuracy){
    H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, "
        << "float accuracy, std::vector<HeatTerm*> heatterms, "
        << "Matter *& sample, PlasmaData const *&pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    Defaults();
    HeatTerms = heatterms;
    CreateFile(filename,false);
}

HeatingModel::HeatingModel(std::string filename, float accuracy, 
std::vector<HeatTerm*> heatterms,Matter *& sample, PlasmaGrid_Data &pgrid):
Model(filename,sample,pgrid,accuracy){
    H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, "
        << "float accuracy, std::vector<HeatTerm*> heatterms, "
        << "Matter *& sample, PlasmaGrid const &pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    Defaults();
    HeatTerms = heatterms;
    CreateFile(filename,false);
}

HeatingModel::HeatingModel(std::string filename, float accuracy, 
std::vector<HeatTerm*> heatterms, Matter *& sample, PlasmaGrid_Data &pgrid,
PlasmaData &pdata):
Model(filename,sample,pgrid,pdata,accuracy){
    H_Debug("\n\nIn HeatingModel::HeatingModel(std::string filename, "
        << "float accuracy, std::vector<HeatTerm*> heatterms, "
        << "Matter *& sample, PlasmaGrid const &pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    Defaults();
    HeatTerms = heatterms;
    CreateFile(filename,false);
}

void HeatingModel::Defaults(){
    H_Debug("\tIn HeatingModel::Defaults()\n\n");
    PowerIncident = 0;                      //!< kW, Power Incident
    OldTemp = Sample->get_temperature();    //!< Set default OldTemp
    ThermalEquilibrium = false;
}

void HeatingModel::CreateFile(std::string filename){
    H_Debug("\tIn HeatingModel::CreateFile(std::string filename)\n\n");
    CreateFile(filename,false);
}

void HeatingModel::CreateFile(std::string filename, bool PrintPhaseData){
    H_Debug("\tIn HeatingModel::CreateFile(std::string filename, "
        << "bool PrintPhaseData)\n\n");
    FileName=filename;
    ModelDataFile.open(FileName);
    ModelDataFile << std::scientific << std::setprecision(16) << std::endl;
    ModelDataFile << "Time\tTemp\tMass\tDensity";

    if( PrintPhaseData )                
        ModelDataFile << "\tFusionE\tVapourE";
    if( Sample->get_c(1) == 'v' || Sample->get_c(1) == 'V' )    
        ModelDataFile << "\tLinearExpansion";
    if( Sample->get_c(2) == 'v' || Sample->get_c(2) == 'V' )    
        ModelDataFile << "\tCv";
    if( Sample->get_c(3) == 'v' || Sample->get_c(3) == 'V' )    
        ModelDataFile << "\tVapourP";
    
    //!< Loop over heat terms and print their names
    for(auto iter = HeatTerms.begin(); iter != HeatTerms.end(); ++iter) {
        ModelDataFile << "\t" << (*iter)->PrintName();
        if( (*iter)->PrintName() == "EmissivityModel" &&
            (Sample->get_c(0) == 'f' || Sample->get_c(0) == 'F') )   
            ModelDataFile << "\tEmissiv";
    }

    ModelDataFile << "\n";
    Print();
    ModelDataFile.close();
    ModelDataFile.clear();
}

double HeatingModel::ProbeTimeStep()const{
    H_Debug( "\tIn HeatingModel::ProbeTimeStep()\n\n" );

    //!< Take Eularian step to get initial time step
    H_Debug("\t"); 
    double TotalPower = CalculatePower(Sample->get_temperature());
    double timestep = fabs((Sample->get_mass()*Sample->get_heatcapacity()*
        Accuracy)/TotalPower);

    //!< Calculate timestep that produces mass change of less than 0.01% 
    //!< of current mass.
    //!< If this timestep is quicker than current step, change timestep
    if( TotalPower != 0 ){
        for(auto iter = HeatTerms.begin(); iter != HeatTerms.end(); ++iter) {
            if( (*iter)->PrintName() == "EvaporationModel" ){
                if( Sample->is_liquid() ){
                    double MassTimeStep = fabs((0.01*Sample->get_mass()*AvNo)
                        /fabs(Flux::EvaporationFlux(Sample,Pdata,
                        Sample->get_temperature())*Sample->get_atomicmass()));
                    if( MassTimeStep < timestep ){
                        H_Debug("\nMass is limiting time step\nMassTimeStep = " 
                            << MassTimeStep << "\ntimestep = " << timestep);
                        timestep = MassTimeStep;
                    }
                }
            }
        }
    }else{
        //!< Check thermal equilibrium hasn't been explicitly reached somehow.
        if( ContinuousPlasma ){ 
            static bool runOnce = true;
            WarnOnce(runOnce,"\nWarning! TotalPower = 0");
            std::cout << "\nThermalEquilibrium reached on condition (1): "
                << "TotalPower = 0.";
            timestep = 1;
        }else if( !ContinuousPlasma ){
            std::cout << "\nNo Net Power Region...";
            timestep = 10;
        }
    }
    //!< Check Thermal Equilibrium hasn't been reached for continuous plasma
    double DeltaTempTest = TotalPower*timestep/(Sample->get_mass()*
        Sample->get_heatcapacity());
    //!< If we're not boiling and in a continuous Plasma
    if( Sample->get_temperature() != Sample->get_superboilingtemp() && 
        Sample->get_temperature() != Sample->get_meltingtemp() 
        && ContinuousPlasma ){

        //!< If temperature changed sign this step
        if( ((Sample->get_temperature()-OldTemp > 0 && DeltaTempTest < 0) 
            || (Sample->get_temperature()-OldTemp < 0 && DeltaTempTest > 0)) ){
            std::cout << "\n\nThermal Equilibrium reached on condition (2):"
                << " Sign change of Temperature change!";
            timestep = 1;
        //!< If Temperature gradient is less than 1%
        }
        if( (fabs(DeltaTempTest/timestep) < 0.01) ){ 
            std::cout << "\n\nThermal Equilibrium reached on condition (3): "
                << " Temperature Gradient < 0.01!";
            timestep = 1;
        }
    }

    H1_Debug("\n\n\ttimestep = " << timestep << "\n\tSample->get_mass() = " 
        << Sample->get_mass() << "\n\tSample->get_heatcapacity() = " << 
        Sample->get_heatcapacity() << "\n\tTotalPower = " << TotalPower 
        << "\n\tAccuracy = " << Accuracy);
    assert(timestep > 0 && timestep != INFINITY && timestep == timestep);

    return timestep;
}

double HeatingModel::UpdateTimeStep(){
    H_Debug( "\tIn HeatingModel::UpdateTimeStep()\n\n" );

    TimeStep = ProbeTimeStep();
    OldTemp = Sample->get_temperature();
    if( TimeStep == 1 ) ThermalEquilibrium = true;

    return TimeStep;
}

void HeatingModel::Print(){
    H_Debug("\tIn HeatingModel::Print()\n\n");
    ModelDataFile.open(FileName,std::ofstream::app);
    ModelDataFile   << TotalTime << "\t" << Sample->get_temperature() << "\t" 
        << Sample->get_mass() << "\t" << Sample->get_density();

    //!< Print variable constants if they are varying
    if( Sample->get_c(1) == 'v' || Sample->get_c(1) == 'V' )    
        ModelDataFile << "\t" << Sample->get_linearexpansion();
    if( Sample->get_c(2) == 'v' || Sample->get_c(2) == 'V' )    
        ModelDataFile << "\t" << Sample->get_heatcapacity();
    if( Sample->get_c(3) == 'v' || Sample->get_c(3) == 'V' )    
        ModelDataFile << "\t" << Sample->get_vapourpressure();

    //!< Loop over heat terms and print their values
    for(auto iter = HeatTerms.begin(); iter != HeatTerms.end(); ++iter) {
        if( (*iter)->PrintName() == "EmissivityModel" ){
            ModelDataFile << "\t" << (*iter)->
                Evaluate(Sample, Pdata, Sample->get_temperature());
            if (Sample->get_c(0) == 'f' || Sample->get_c(0) == 'F'){
                ModelDataFile   << "\t" << Sample->get_emissivity();
            }
        }else if( (*iter)->PrintName() == "EvaporationModel" ){
            if( Sample->is_liquid() ){
                ModelDataFile << "\t" << (*iter)->
                    Evaluate(Sample, Pdata, Sample->get_temperature())*1000;
            }else{ //!< If evaporation is turned off
                ModelDataFile   << "\t" << 0;
            }
        }else{
            ModelDataFile << "\t" << (*iter)->
                Evaluate(Sample, Pdata, Sample->get_temperature());;
        }
    }

    ModelDataFile << "\n";
    ModelDataFile.close();
    ModelDataFile.clear();
}



const int HeatingModel::Vapourise(){
    H_Debug("\tIn HeatingModel::Vapourise()\n\n");
    //!< If the sample is gaseous or in TE 
    //!< (Given that the plasma is continuous), the model ends.
    while( !Sample->is_gas() ){
        Heat();
        UpdateTimeStep();
        if( ContinuousPlasma && ThermalEquilibrium )
            break;
    }

    int rValue(0);
    if( Sample->is_gas() 
        && Sample->get_superboilingtemp() <= Sample->get_temperature() ){
        std::cout << "\n\nSample has Boiled ";
        rValue = 1;
    }else if( Sample->is_gas() 
        && Sample->get_superboilingtemp() > Sample->get_temperature() ){
        std::cout << "\n\nSample has Evaporated ";
        rValue = 2;
    }else if( ThermalEquilibrium && ContinuousPlasma ){
        std::cout << "\n\nSample has reached Thermal Equilibrium in "
            << "Continuous Plasma.";
        rValue = 3;
    }
    std::cout << "\nat T = " << Sample->get_temperature() << "K in " 
        << TotalTime << "s!\n\n*********\n\n";
    ModelDataFile.close();
    //!<  0, running normally. 1; Sample boiled. 
    //!< 2; Sample Evaporated. 3; Thermal equilibrium.
    return rValue; 
}

void HeatingModel::UpdateRERN(){
    //!< If it's positive, Ions aren't backscattered
    double RE(0.0), RN(0.0);
    if( Sample->get_potential() >= 0.0 ){
        backscatter(Pdata->ElectronTemp,Pdata->IonTemp,Pdata->mi,
            Sample->get_potential(),Sample->get_elem(),RE,RN);
    }

    if( RN > 0.1 ){ //!< Uncomment when RN is calculated
        static bool runOnce = true;
        std::string Warning = "In HeatingModel::UpdateRERN()\nRN > 0.1. ";
        Warning += "Neutral Recombination affected by backscattering by more ";
        Warning += "than 10%!";
        WarnOnce(runOnce,Warning);
    }   
    if( RE > 0.1 ){ //!< Uncomment when RE is calculated
        static bool runOnce = true;
        std::string Warning = "In HeatingModel::UpdateRERN()\nRE > 0.1. ";
        Warning += "Ion Heat Flux affected by backscattering by more than 10%!";
        WarnOnce(runOnce,Warning);
    }
    Sample->set_rern(RE,RN);
}

void HeatingModel::Heat(double timestep){
    H_Debug("\tIn HeatingModel::Heat(double timestep)\n\n");
    
    //!< Make sure timestep input time is valid. Shouldn't exceed the timescale 
    //!< of the process.
    assert(timestep > 0 && timestep <= TimeStep );
    assert( Sample->get_mass() > 0 );   
    
    //!< Calculate total energy through RungeKutta4 method
    double TotalEnergy = RungeKutta4(timestep);
    H1_Debug( "\tTotalEnergy = " << TotalEnergy << "kJ\n");
    Sample->update_temperature(TotalEnergy);  //!< Update Temperature

    //!< Account for evaporative mass loss, if model is turned on, if it's a 
    //!< liquid and not boiling!
    for(auto iter = HeatTerms.begin(); iter != HeatTerms.end(); ++iter) {
        if( (*iter)->PrintName() == "EvaporationModel" )
            if( Sample->is_liquid()
                && (Sample->get_temperature() != Sample->get_boilingtemp()) )
                Sample->update_mass( 
                    (timestep*Flux::EvaporationFlux(Sample,Pdata,Sample->get_temperature())*
                    Sample->get_atomicmass())/AvNo );
    }

    H1_Debug("\n\t\tMass Loss = " 
        << (timestep*Flux::EvaporationFlux(Sample,Pdata,
        Sample->get_temperature())*Sample->get_atomicmass())/AvNo << "\n");

    if( !Sample->is_gas() )
        Sample->update();

    Print();  //!< Print data to file
    H_Debug("\t"); 

    TotalTime += timestep;
}

void HeatingModel::Heat(){
    H_Debug("\tIn HeatingModel::Heat()\n\n");
    Heat(TimeStep);
}

double HeatingModel::CalculatePower(double DustTemperature)const{
    H_Debug( "\tIn HeatingModel::CalculatePower(double DustTemperature = " 
        << DustTemperature << ")\n\n");
    //!< Rreduces the number of divisions
    double TotalPower = PowerIncident*1000;
    H1_Debug("\n\n\t\tPowerIncident = \t"    << PowerIncident*1000 << "W");
    
    //!< Loop over heat terms and print their names
    for(auto iter = HeatTerms.begin(); iter != HeatTerms.end(); ++iter) {
        if( (*iter)->PrintName() == "EvaporationModel" ){
            if( Sample->is_liquid() ){
                TotalPower += (*iter)
                    ->Evaluate(Sample, Pdata, DustTemperature)*1000;
            }
        }else{
            TotalPower 
                += (*iter)->Evaluate(Sample, Pdata, DustTemperature);
        }
        H1_Debug("\n\t\t" << (*iter)->PrintName() << " = "  
            << (*iter)->Evaluate(Sample, Pdata, DustTemperature)  << "W");
    }
    TotalPower = TotalPower/1000;

    H1_Debug("\n\t\tTotalPower = \t" << TotalPower*1000 << "W\n\n");
    return TotalPower;
}

double HeatingModel::RungeKutta4(double timestep){
    H_Debug( "\tIn HeatingModel::RungeKutta4(double timestep)\n\n");
    double k1 = CalculatePower(Sample->get_temperature()); 
    if(k1<0 && fabs(k1/2) > Sample->get_temperature()){
        std::cout << "\n\nThermal Equilibrium reached on condition (4):"
            << " k1 step negative and larger than Td!";
        ThermalEquilibrium = true;
        return 0;
    }
    double k2 = CalculatePower(Sample->get_temperature()+k1/2); 
    if( k2<0 && fabs(k2/2) > Sample->get_temperature() ){
        std::cout << "\n\nThermal Equilibrium reached on condition (4):"
            << " k2 step negative and larger than Td!";
        ThermalEquilibrium = true;
        return (timestep/6)*k1;
    }
    double k3 = CalculatePower(Sample->get_temperature()+k2/2);
    if( k3<0 && fabs(k3) > Sample->get_temperature() ){
        std::cout << "\n\nThermal Equilibrium reached on condition (4):"
            << " k3 step negative and larger than Td!";
        ThermalEquilibrium = true;
        return (timestep/6)*(k1+2*k2);
    }
    double k4 = CalculatePower(Sample->get_temperature()+k3);
    if( k4<0 && fabs(k4/2) > Sample->get_temperature() ){
        std::cout << "\n\nThermal Equilibrium reached on condition (4):"
            << " k4 step negative and larger than Td!";
        ThermalEquilibrium = true;
        return (timestep/6)*(k1+2*k2+2*k3);
    }
    H1_Debug( "\n\t\ttimestep = " << timestep << "\n\t\tk1 = " << k1 
        << "\n\t\tk2 =" << k2 << "\n\t\tk3 = " << k3 << "\n\t\tk4 = " << k4
        << "\n\n");
    return (timestep/6)*(k1+2*k2+2*k3+k4);
};
