/** @file ChargingModel.cpp
 *  @brief Implementation of class for physics models relevant to dust charging
 *  
 *  Implement the member functions of the ChargeModel class
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 */

#include "ChargingModel.h"

ChargingModel::ChargingModel():Model(){
    C_Debug("\n\nIn ChargingModel::ChargingModel():Model()\n\n");
    // Charging Models turned on of possible 3

    CurrentTerms.push_back(new Term::OMLe());
    CurrentTerms.push_back(new Term::OMLi());
    CreateFile("Data/default_cm_0.txt");
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> currentterms, Matter *& sample, PlasmaData &pdata):
Model(filename,sample,pdata,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << " float accuracy, std::vector<CurrentTerm*> currentterms, "
        << "Matter *& sample, PlasmaData *&pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    CurrentTerms = currentterms;
    CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> currentterms, Matter *& sample, PlasmaData *pdata):
Model(filename,sample,*pdata,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << "float accuracy, std::vector<CurrentTerm*> currentterms, "
        << "Matter *& sample, PlasmaData *&pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    CurrentTerms = currentterms;
    CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> currentterms, Matter *& sample, PlasmaGrid_Data &pgrid):
Model(filename,sample,pgrid,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << "float accuracy, std::vector<CurrentTerm*> currentterms, "
        << "Matter *& sample, PlasmaGrid_Data &pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    CurrentTerms = currentterms;
    CreateFile(filename);
}

ChargingModel::ChargingModel(std::string filename, float accuracy, 
std::vector<CurrentTerm*> currentterms, Matter *& sample, PlasmaGrid_Data &pgrid, 
PlasmaData &pdata):
Model(filename,sample,pgrid,pdata,accuracy){
    C_Debug("\n\nIn ChargingModel::ChargingModel(std::string filename, "
        << "float accuracy, std::vector<CurrentTerm*> currentterms, "
        << "Matter *& sample, PlasmaGrid_Data &pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    CurrentTerms = currentterms;
    CreateFile(filename);
}

void ChargingModel::CreateFile(std::string filename){
    C_Debug("\tIn ChargingModel::CreateFile(std::string filename)\n\n");
    FileName = filename;
    ModelDataFile.open(FileName);
    ModelDataFile << std::scientific << std::setprecision(16) << std::endl;
//    ModelDataFile << "Time\tCharge\tSign\tDeltatot\tPotential\n";
    ModelDataFile << "Time\tCharge\tDeltatot\tPotential\n";
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
//    if( Sample->is_positive() )  ModelDataFile << "\t//Pos";
//    if( !Sample->is_positive() ) ModelDataFile << "\t//Neg";
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
    assert(timestep > 0);
    
    //!< Calculate Thermionic and secondary electron emission yields for use in
    //!< the heating models
    double DTherm(0.0), DSec(0.0);
    for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); ++iter) {
        if( (*iter)->PrintName() == "TEEcharge" ){
            DTherm = Flux::DeltaTherm(Sample,Pdata);
        }else if( (*iter)->PrintName() == "TEESchottky" ){
            DTherm = Flux::DeltaTherm(Sample,Pdata);
        }
        if( (*iter)->PrintName() == "SEEcharge" ){
            DSec = Flux::DeltaSec(Sample,Pdata);
        }
    }

    //!< Implement Bisection method to find root of current balance
    double Potential = Bisection();

    //!< Implement regular falsi method to find root of current balance
    //double Potential = RegulaFalsi();

    //!< Have to calculate charge of grain here since it doesn't know about the 
    //!< Electron Temp and since potential is normalised.
    if( Potential != 0 ){
        double charge = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*
            Pdata->ElectronTemp)/echarge;
        Sample->update_charge(charge,Potential,DTherm,DSec);
    }
    //!< Increment total time recorded by the model
    TotalTime += timestep;

    C_Debug("\t"); Print();
}

void ChargingModel::Charge(){
    C_Debug("\tIn ChargingModel::Charge()\n\n");
    Charge(TimeStep);
}

double ChargingModel::Bisection()const{
    C_Debug("\tIn ChargingModel::Bisection()\n\n");
    //!< Implement Bisection method to find root of current balance
    //!< Initial range of normalised potential based on mostly negative dust
    //!< typically with low magntiudes of potential. This may fail in unusual
    //!< cases.
    double a(-5.0), b(10.0);
    double Current1(0.0), Current2(0.0), Potential(0.0);
    int i(0), imax(1000);
    do{ //!< Do while difference in bounds is greater than accuracy
        //!< Take new x position as halfway between upper and lower bound
        Potential = (a+b)/2.0;

        //!< Sum all the terms in the current balance
        for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); 
            ++iter) {
            if( (*iter)->PrintName() == "SEEcharge" ){
                Current1 += (*iter)->Evaluate(Sample,Pdata,Potential)
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential);
                Current2 += (*iter)->Evaluate(Sample,Pdata,a)
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,a);
                C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                    << (*iter)->Evaluate(Sample,Pdata,Potential)
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential)  
                    << "\n" );
            }else{
                Current1 += (*iter)->Evaluate(Sample,Pdata,Potential);
                Current2 += (*iter)->Evaluate(Sample,Pdata,a);
                C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                    << (*iter)->Evaluate(Sample,Pdata,Potential) << "\n");
            }
        }
        //!< If the root is on the RHS of our midpoint
        if( Current1*Current2 > 0.0 )
            a = Potential;
        else //!< Else, it must be on LHS of our midpoint
            b = Potential;

        //!< Ensure we don't loop forever
        if( i > imax ){
            std::cerr << "ChargingModel::Bisection Root Finding failed to "
                << "converge in " << imax << " steps! Setting Potential = 0.0\n";
            Potential = 0.0;
            return Potential;
        }
        i ++; //!< Increment loop counter
    }while( fabs((b-a)/2.0) > Accuracy && Current1 != 0.0 );
    return Potential;
}

double ChargingModel::RegulaFalsi()const{
    C_Debug("\tIn ChargingModel::RegulaFalsi()\n\n");
    //!< Implement regula falsi method to find root of current balance
    //!< Initial range of normalised potential based on mostly negative dust
    //!< typically with low magntiudes of potential. This may fail in unusual
    //!< cases.

    double Current1(0.0), Current2(0.0);
    double a(-5.0), b(10.0);
    double amin(-5.0), bmax(10.0);
    //!< Sum all the terms in the current balance
    for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); 
        ++iter) {
        if( (*iter)->PrintName() == "SEEcharge" ){
            Current1 += (*iter)->Evaluate(Sample,Pdata,a)
                *CurrentTerms[0]->Evaluate(Sample,Pdata,a);
            Current2 += (*iter)->Evaluate(Sample,Pdata,b)
                *CurrentTerms[0]->Evaluate(Sample,Pdata,b);
            C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                << (*iter)->Evaluate(Sample,Pdata,a)
                *CurrentTerms[0]->Evaluate(Sample,Pdata,a)  
                << "\n" );
        }else{
            Current1 += (*iter)->Evaluate(Sample,Pdata,a);
            Current2 += (*iter)->Evaluate(Sample,Pdata,b);
            C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                << (*iter)->Evaluate(Sample,Pdata,a) << "\n");
        }
    }

    int side(0), imax(10000);
    double Potential(0.0);
    //!< Loop until we reach imax iterations
    for (int i = 0; i < imax; i++){
        //!< following line is regula falsi method
        Potential = (b*Current1-a*Current2)/(Current1-Current2);
        //!< If we're within accuracy, return result
        if (fabs(b-a) < Accuracy*fabs(b+a) && b < bmax && a > amin )
            return Potential;
        double Current3(0.0);

        //!< Sum all the terms in the current balance
        for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); 
            ++iter){
            if( (*iter)->PrintName() == "SEEcharge" ){
                Current1 += (*iter)->Evaluate(Sample,Pdata,Potential)
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential);
                C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                    << (*iter)->Evaluate(Sample,Pdata,Potential)
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential)  
                    << "\n" );
            }else{
                Current3 += (*iter)->Evaluate(Sample,Pdata,Potential);
                C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                    << (*iter)->Evaluate(Sample,Pdata,Potential) << "\n");
            }
        }

        if(Current3 * Current2 > 0){
            //!< Current3 and Current2 have same sign, copy Potential to b
            b = Potential; Current2 = Current3;
            if(side==-1){ 
                Current1 /= 2;
            }
            side = -1;
        }else if(Current3 * Current1 > 0){
            //!< Current3 and Current1 have same sign, copy Potential to a
            a = Potential;  Current1 = Current3;
            if(side==+1){
                Current2 /= 2;
            }
            side = +1;
        }else{
            //!< Current3 * Current is very small (looks like zero)
            if( b < bmax && a > amin  )
                return Potential;
        }
    }
    //!< If we reach this point we've looped more than imax times.
    std::cerr << "ChargingModel::RegulaFalsi Root Finding failed to "
                << "converge in " << imax << " steps! Setting Potential = 0.0\n";
    Potential = 0.0;
    return Potential;
}
