/** @file ChargingModel.cpp
 *  @brief Implementation of class for physics models relevant to dust charging
 *  
 *  Implement the member functions of the ChargeModel class
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 */

#include "ChargingModel.h"

ChargingModel::ChargingModel():
ChargeModel(new Term::solveOML()),Model(){
    C_Debug("\n\nIn ChargingModel::ChargingModel():Model()\n\n");
    // Charging Models turned on of possible 3

    CreateFile("Data/default_cm_0.txt");
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, PlasmaData &pdata):
Model(filename,sample,pdata,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << " float accuracy, std::vector<CurrentTerm*> CurrentTerms, "
        << "Matter *& sample, PlasmaData *&pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    ChargeModel = chargemodel;
    CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, PlasmaData *pdata):
Model(filename,sample,*pdata,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << "float accuracy, std::vector<CurrentTerm*> CurrentTerms, "
        << "Matter *& sample, PlasmaData *&pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    ChargeModel = chargemodel;
    CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, PlasmaGrid_Data &pgrid):
Model(filename,sample,pgrid,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << "float accuracy, std::vector<CurrentTerm*> CurrentTerms, "
        << "Matter *& sample, PlasmaGrid_Data &pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    ChargeModel = chargemodel;
    CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> CurrentTerms, Matter *& sample, PlasmaGrid_Data &pgrid, 
PlasmaData &pdata):
Model(filename,sample,pgrid,pdata,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << "float accuracy, std::vector<CurrentTerm*> CurrentTerms, "
        << "Matter *& sample, PlasmaGrid_Data &pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    ChargeModel = chargemodel;
    CreateFile(filename);
}

void ChargingModel::CreateFile(std::string filename){
    C_Debug("\tIn ChargingModel::CreateFile(std::string filename)\n\n");
    FileName = filename;
    ModelDataFile.open(FileName);
    ModelDataFile << std::scientific << std::setprecision(16) << std::endl;
    ModelDataFile << "Time\tCharge\tSign\tDeltatot\tPotential\n";
    ModelDataFile.close();
    ModelDataFile.clear();
    Print();
}

void ChargingModel::Print(){
    C_Debug("\tIn ChargingModel::Print()\n\n");
    ModelDataFile.open(FileName,std::ofstream::app);
    ModelDataFile << TotalTime << "\t" 
        << -(4.0*PI*epsilon0*Sample->get_radius()*Sample->get_potential()*Kb*
        Pdata->ElectronTemp)/(echarge*echarge);
    if( Sample->is_positive() )  ModelDataFile << "\tPos";
    if( !Sample->is_positive() ) ModelDataFile << "\tNeg";
    ModelDataFile << "\t" << Sample->get_deltatot() << "\t" 
        << Sample->get_potential() << "\n";

    ModelDataFile.close();
    ModelDataFile.clear();
}

double ChargingModel::ProbeTimeStep()const{
    C_Debug( "\tIn ChargingModel::ProbeTimeStep()\n\n" );

    double timestep(1.0);

    //!< Tests have shown that Time step based on the electron plasma frequency. 
    if( Pdata->ElectronDensity != 0 )
        timestep = Accuracy*sqrt((epsilon0*Me)/
            (2*PI*Pdata->ElectronDensity*echarge*echarge));

    assert(timestep == timestep);
    assert(timestep > 0);

    return timestep;
}

double ChargingModel::UpdateTimeStep(){
    C_Debug( "\tIn ChargingModel::UpdateTimeStep()\n\n" );
    TimeStep = ProbeTimeStep();
    return TimeStep;
}

void ChargingModel::Charge(double timestep){
    C_Debug("\tIn ChargingModel::Charge(double timestep)\n\n");

    //!< Make sure timestep input time is valid. Shouldn't exceed the timescale
    //!< of the process.
    assert(timestep > 0);// && timestep <= TimeStep );
    
    
    double Potential = Sample->get_potential();

    // Evaluate the charging model
    double guess(0.0);
    do{
        double Current(0.0);
        for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); ++iter) {
            Current += (*iter)->Evaluate(Sample,Pdata,Potential);
            C_Debug( "\n\t\t" << (*iter)->PrintName << "=" 
                << (*iter)->Evaluate(Sample,Pdata) );
        }

        // IMPLEMENT ROOT FINDING METHOD, BISECTION?
        guess = Potential - (Current/CurrentDiff);
        Potential = guess;
        i ++;
    }while( fabs(Potential-guess) > Accuracy );


    double DSec = Flux::DeltaSec(Sample,Pdata);
    double DTherm = Flux::ThermFluxSchottky(Sample,Pdata,Potential)/
            Flux::OMLElectronFlux(Pdata,Potential);
    if( ChargeModel->PrintName() == "PHL" ){ 
        DTherm = Flux::ThermFluxSchottky(Sample,Pdata,Potential)/
            Flux::PHLElectronFlux(Sample,Pdata,Potential);
    }else if( ChargeModel->PrintName() == "DTOKSOML" ||
        ChargeModel->PrintName() == "DTOKSWell" ||
        ChargeModel->PrintName() == "MOMLWEM" ){
        DTherm = Flux::DeltaTherm();
    }

    //!< Have to calculate charge of grain here since it doesn't know about the 
    //!< Electron Temp and since potential is normalised.
    //!< This information has to be passed to the grain.
    double charge = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*
        Pdata->ElectronTemp)/echarge;
    Sample->update_charge(charge,Potential,DTherm,DSec);
    TotalTime += timestep;

    C_Debug("\t"); Print();
}

void ChargingModel::Charge(){
    C_Debug("\tIn ChargingModel::Charge()\n\n");
    Charge(TimeStep);
}
/*
void ChargingModel::Charge(double timestep){
     //!< Make sure timestep input time is valid. Shouldn't exceed the timescale
     //!< of the process.
     assert(timestep > 0);// && timestep <= TimeStep );
    double DSec = 0.0;
    double DTherm = 0.0;

    double Potential;
    if( UseModel[0] ){  //!< Dynamically choose charging model
        double DebyeRatio = Sample->get_radius()/
        sqrt((epsilon0*Kb*Pdata->ElectronTemp)/
            (Pdata->ElectronDensity*pow(echarge,2)));

        DSec = DeltaSec();
        DTherm = ThermFluxSchottky(Potential)/
            OMLElectronFlux(Sample->get_potential());
        if( (DSec+DTherm) <= 0.1 ){//!< Emission IS NOT important
            Potential = solveTHS();
        }else{ //!< Emission IS important
            //!< Small dust grains wrt the debye length
            if( DebyeRatio <= 1.0 ){ 
                //!< OMLWEM like Nikoleta's theory MOML-EM, DOESN'T EXIST YET
                //!< So instead, we do SOMLWEM
                Potential = solveSOML( Sample->get_potential());
            }else{ // large dust grains wrt the debye length
                //!< MOML-EM since emission is far more important than 
                //!< magnetic field effects HASN'T BEEN IMPLEMENTED YET, 
                //!< solve SMOMLWEM
                //std::cout << "\nSolving SMOMLWEM!";
                Potential = solveSMOML(Sample->get_potential());
            }
        }
        
    }else if( UseModel[1] ){ //!< OML Charging model
        Potential = solveOML(0.0,Potential);
        if( Potential < 0 ){ //!< In this case we have a positive grain!
            Potential = solvePosOML( 0.0, Sample->get_potential());
        }
    }else if( UseModel[2] ){ //!< MOML Charging model
        Potential = solveMOML();
        if( Potential < 0 ){ //!< In this case we have a positive grain!
            Potential = solvePosOML( 0.0, Sample->get_potential());
        }
    }else if( UseModel[3] ){ //!< SOML Charging model
    
    
    double Potential = Sample->get_potential();

        Potential = solveSOML(Sample->get_potential());
        DSec = DeltaSec();
        DTherm = ThermFluxSchottky(Potential)/
            OMLElectronFlux(Sample->get_potential());
    }else if( UseModel[4] ){ //!< SMOML Charging model
        Potential = solveSMOML(Sample->get_potential());
        DSec = DeltaSec();
        DTherm = ThermFluxSchottky(Potential)/
            OMLElectronFlux(Sample->get_potential());
    }else if(UseModel[5] ){  //!< CW Charging Model, Willis fit to Sceptic
        Potential = solveCW(Sample->get_potential());
        DSec = DeltaSec();
        DTherm = ThermFluxSchottky(Potential)/
            OMLElectronFlux(Sample->get_potential());
    }else if( UseModel[6] ){ //!< Use Patterchini, Hutchinson and Lapenta model
        Potential = solvePHL(Sample->get_potential());
        DSec = DeltaSec();
        DTherm = ThermFluxSchottky(Potential)/
            PHLElectronFlux(Sample->get_potential());
    }else if( UseModel[7] ){ //!< In this case, use THS Model
        Potential = solveTHS();
        DSec = DeltaSec();
        DTherm = ThermFluxSchottky(Potential)/
            PHLElectronFlux(Sample->get_potential());
    }else if( UseModel[8] ){ //!< Original DTOKS Charging scheme
        //!< WARNING! THIS SCHEME FOR THE CHARGING MODEL CREATES DISCONTINUITIES
        //!< WHEN FORMING A WELL!
        DSec = DeltaSec();
        DTherm = DeltaTherm();
        //!< If electron emission yield exceeds unity, we have potential well...
        if( (DSec + DTherm) >= 1.0 ){   
            //!< Take away factor of temperature ratios for depth of well
            //!< This is not explained further...
            Potential = solveOML( 0.0, Sample->get_potential()) 
                - (Kb*Sample->get_temperature())/(Pdata->ElectronTemp*Kb); 
        }else{ //!< If the grain is negative...
            //!< Calculate the potential with the normal current balance 
            //!< including electron emission
            Potential = solveOML( DSec + DTherm,Sample->get_potential());
            if( Potential < 0.0 ){
                //!< But! If it's now positive, our assumptions must be wrong!
                //!< So now we assume it's positive and calculate the potential 
                //!< with a well.
                Potential = solveOML(0.0,Sample->get_potential())-
                    (Kb*Sample->get_temperature())/(Pdata->ElectronTemp*Kb);
            }
        }
    //!< In this case, maintain a potential well for entire temperature range
    }else if( UseModel[9] ){ 
        Potential = solveOML( 0.0, Sample->get_potential()) 
            - Kb*Sample->get_temperature()/(Pdata->ElectronTemp*Kb);
        //!< In this case we have a positive grain! Assume well disappears
        if( Potential < 0 ){ 
            Potential = solvePosOML( 0.0, Sample->get_potential());
        }
    }else if( UseModel[10] ){ //!< MOMLWEM Charging model
        DSec = DeltaSec();
        DTherm = DeltaTherm();
        
        if( (DSec + DTherm) >= 1.0 ){//!< In this case we have a positive grain!
            Potential = solvePosOML(0.0, Sample->get_potential());
        }else{
            Potential = solveMOMLWEM(DSec + DTherm);
            if( Potential < 0.0 ){
                Potential = solvePosOML(0.0, Sample->get_potential());
            }
        }
//!< Following Semi-empirical fit to Sceptic results as detailed in Chris Willis
//!< Thesis,
//!< https://spiral.imperial.ac.uk/handle/10044/1/9329, pages 68-70
double ChargingModel::solveCW(double guess){
    C_Debug("\tIn ChargingModel::solveCW()\n\n");
    double A = Pdata->mi;
    double b = Pdata->IonTemp/Pdata->ElectronTemp;
    double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/
        (Pdata->ElectronDensity*pow(echarge,2)));
    double Rho = Sample->get_radius()/DebyeLength;
    double Rho_OML = 1.25*pow(b,0.4);
    double Rho_Upper = 50.0;
    double Potential(0.0);
    if( Rho <= Rho_OML ){ //!< This is the OML Limit
        if( b <=2 ){ //!< Ti <= 2.0*Te
            Potential = 0.405*log(Pdata->A)+(0.253+0.021*log(Pdata->A))*log(b)+
                2.454;
        }else{ //!< Ti > 2.0*Te
            Potential = 0.401*log(Pdata->A)+(-0.122+0.029*log(Pdata->A))*log(b)+
                2.698;
        }
        //!< This is the transition region
    }else if( Rho <= Rho_Upper && Rho > Rho_OML ){ 
        double Gradient = (log(Rho/Rho_Upper)/log(Rho_Upper/Rho_OML))+1.0;
        double DeltaPhi = 0.5*log(2.0*PI*(Me/Pdata->mi)*(1.0+5.0*b/3.0))*
            Gradient;
        int i(0);
        do{
            //!< If first, don't update guess. Avoids numerical instabilities.
            if( i > 0 ) 
                guess = Potential;
            Potential = guess - ( ( exp(-guess) - sqrt(b*(Me/Pdata->mi))*
                (1+guess/b-DeltaPhi/b))/(-exp(-guess)-sqrt((Me/Pdata->mi)/b)));
            i ++;
        }while( fabs(guess-Potential) > Accuracy );
        return guess;
    }else if( Rho > Rho_Upper ){ //!< This is the MOML limit
        if( b <=2 ){ //!< Ti <= 2.0*Te
            Potential = 0.456*log(Pdata->A)+3.179;
        }else{ //!< Ti > 2.0*Te
            Potential = 0.557*log(Pdata->A)-(0.386+0.024*log(Pdata->A))*log(b)+
                3.399;
        }
    }else{
        std::cerr << "\nError in ChargingModel::solveCW()!"
            << " Rho is poorly defined!\n\n";
    }
    return Potential;
}
*/