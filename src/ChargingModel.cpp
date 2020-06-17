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
    ModelDataFile << "Time\tCharge\tDeltatot\tPotential";

    for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); ++iter)
        ModelDataFile << "\t" << (*iter)->PrintName();
 
	ModelDataFile << "\n";
    ModelDataFile.close();
    ModelDataFile.clear();
    Print();
}

void ChargingModel::Print(){
    C_Debug("\tIn ChargingModel::Print()\n\n");
    ModelDataFile.open(FileName,std::ofstream::app);
    double Potential = Sample->get_potential();
    ModelDataFile << TotalTime << "\t" 
        << -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*
        Pdata->ElectronTemp)/(echarge*echarge);
//    if( Sample->is_positive() )  ModelDataFile << "\t//Pos";
//    if( !Sample->is_positive() ) ModelDataFile << "\t//Neg";

    ModelDataFile << "\t" << Sample->get_deltatot() << "\t" 
        << Potential;
    for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); 
        ++iter) {
        if( (*iter)->PrintName() == "SEEcharge" ){
            ModelDataFile << "\t" << (*iter)->Evaluate(Sample,Pdata,Potential)
                *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential);
        }else{
            ModelDataFile << "\t" << (*iter)->Evaluate(Sample,Pdata,Potential);
        }
    }
    ModelDataFile << "\n";
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
    //double Potential = Bisection();

    //!< Implement regular falsi method to find root of current balance
    double Potential = RegulaFalsi();

    //!< Have to calculate charge of grain here since it doesn't know about the 
    //!< Electron Temp and since potential is normalised.
    if( Potential != 0 ){
        double charge = -(4.0*PI*epsilon0*Sample->get_radius()*Potential*Kb*
            Pdata->ElectronTemp)/echarge;
        Sample->update_charge(charge,Potential,DTherm,DSec);
    }else{
	    double charge = -(4.0*PI*epsilon0*Sample->get_radius()
		    *Sample->get_potential()*Kb*Pdata->ElectronTemp)/echarge;
        Sample->update_charge(charge,Sample->get_potential(),DTherm,DSec);
	}
    //!< Increment total time recorded by the model
    TotalTime += timestep;

    C_Debug("\t"); 
	if( PrintSteps >= PrintInterval ){
	    Print();
		PrintSteps = 1;
	}else{
        PrintSteps ++;
	}
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
    double a(-500), b(10.0);
    double Current1(0.0), Current2(0.0), Potential(0.0);
//	double MaxCurrent;
    int i(0), imax(100000000);

	//!< in this case, we are near to zero within accuracy
    //if( Sample->get_deltatot() > 1.0 ){
	//    a = -1000.0;
	//	b = 0.0;
	//}else if( Sample->get_deltatot() <= 1.0 ){
	//    a = 0.0;
	//	b = 10.0;
	//}

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
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential) );
            }else{
                Current1 += (*iter)->Evaluate(Sample,Pdata,Potential);
                Current2 += (*iter)->Evaluate(Sample,Pdata,a);
                C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                    << (*iter)->Evaluate(Sample,Pdata,Potential) );
            }
//			if( fabs((*iter)->Evaluate(Sample,Pdata,Potential)) > fabs(MaxCurrent) )
//				MaxCurrent = (*iter)->Evaluate(Sample,Pdata,Potential);
        }
        C_Debug("\n\n");
        //!< If the root is on the RHS of our midpoint
        if( Current1*Current2 > 0.0 )
            a = Potential;
        else //!< Else, it must be on LHS of our midpoint
            b = Potential;

        //!< Ensure we don't loop forever
        if( i > imax ){
            std::cerr << "ChargingModel::Bisection Root Finding failed to "
                << "converge in " << imax << " steps! Returning previous potential.\n";
            return Potential;
        }
        i ++; //!< Increment loop counter
		std::cout << "\nCurrent1 = " << Current1 
//		    << "\nMaxCurrent = " << MaxCurrent
		    << "\nPotential = " << Potential
//		    << "\nfabs(Current1/MaxCurrent) = " << fabs(Current1/MaxCurrent)
			<< "\nAccuracy = " << Accuracy << "\n\n";
//    }while( (fabs((b-a)/2.0) > Accuracy || (fabs(Current1/MaxCurrent) > Accuracy)) && Current1 != 0.0 );
    }while( fabs((b-a)/2.0) > Accuracy && Current1 != 0.0 );
    //std::cin.get();
    return Potential;
}

double ChargingModel::RegulaFalsi()const{
    C_Debug("\tIn ChargingModel::RegulaFalsi()\n\n");
    //!< Implement regula falsi method to find root of current balance
    //!< Initial range of normalised potential based on mostly negative dust
    //!< typically with low magntiudes of potential. This may fail in unusual
    //!< cases.

    double Current1(0.0), Current2(0.0), Current3(0.0);
    double a(-20.0), b(5.0);
    double amin(-20.0), bmax(5.0);
//    std::cout << "\n\n*** BEGIN ROOT FINDING ***\n\n"; //std::cin.get();
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

        Current3 = 0.0;
        Potential = (b*Current1-a*Current2)/(Current1-Current2);

        //!< If we're within accuracy, return result
        if (fabs(b-a) < Accuracy*fabs(b+a) && b < bmax && a > amin ){
			//std::cout << "\nPotential difference small, returning potential = " << Potential;
            return Potential;
	    }


        //!< Sum all the terms in the current balance
		C_Debug("\n* Summing Currents *\n");
		C_Debug("* Potential = " << Potential << " *\n\n");
        for(auto iter = CurrentTerms.begin(); iter != CurrentTerms.end(); 
            ++iter){
            if( (*iter)->PrintName() == "SEEcharge" ){
                Current3 += (*iter)->Evaluate(Sample,Pdata,Potential)
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential);
                C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                    << (*iter)->Evaluate(Sample,Pdata,Potential)
                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential)  
                    << "\n" );
//			    if( (*iter)->Evaluate(Sample,Pdata,Potential)
//                    *CurrentTerms[0]->Evaluate(Sample,Pdata,Potential) == 0 ) {
//		            std::cerr << "\nCurrent : " << (*iter)->PrintName();
//			            << "Returned zero! Returning zero!\n";
					//std::cin.get();
//					return 0.0;
//				}
            }else{
                Current3 += (*iter)->Evaluate(Sample,Pdata,Potential);
                C_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
                    << (*iter)->Evaluate(Sample,Pdata,Potential) << "\n");
//				if( (*iter)->Evaluate(Sample,Pdata,Potential) == 0 ){
//		            std::cerr << "\nCurrent : " << (*iter)->PrintName();
//			            << "Returned zero! Returning zero!\n";
					//std::cin.get();
//				    return 0.0;
//			    }
            }
        }

//		std::cout << "\n\nPotential = " << Potential;
//		std::cout << "\na = " << a;
//		std::cout << "\nb = " << b;
//		std::cout << "\nCurrent1 = " << Current1;
//		std::cout << "\nCurrent2 = " << Current2;
//		std::cout << "\nCurrent3 = " << Current3;
//		std::cout << "\nCurrent1-Current2 = " << Current1-Current2;
//		std::cout << "\nb*Current1-a*Current2 = " << b*Current1-a*Current2;

        if(Current3 * Current2 > 0){
            //!< Current3 and Current2 have same sign, copy Potential to b
            b = Potential; Current2 = Current3;
            if(side==-1){ 
                Current1 /= 2.0;
            }
            side = -1;
        }else if(Current3 * Current1 > 0){
            //!< Current3 and Current1 have same sign, copy Potential to a
            a = Potential;  Current1 = Current3;
            if(side==+1){
                Current2 /= 2.0;
            }
            side = +1;
        }else{
            //!< Current3 * Current is very small (looks like zero)
            if( b < bmax && a > amin && Potential > a && Potential < b )
                return Potential;
        }
    }
    //!< If we reach this point we've looped more than imax times.
    std::cerr << "ChargingModel::RegulaFalsi Root Finding failed to "
                << "converge in " << imax << " steps! Setting Potential = 0.0\n";
    Potential = 0.0;
    return Potential;
}
