/** @file ForceModel.cpp
 *  @brief Implementation of class for physics models relevant to dust motion
 *  
 *  Implement the member functions of the ForceModel class
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 */

#include "ForceModel.h"

ForceModel::ForceModel():
Model(){
    F_Debug("\n\nIn ForceModel::ForceModel():Model()\n\n");
    OldTemp = Sample->get_temperature();
    CreateFile("Default_Force_filename.txt");
}

ForceModel::ForceModel(std::string filename, float accuracy, 
std::vector<ForceTerm*> forceterms, Matter *& sample, PlasmaData & pdata):
Model(filename,sample,pdata,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::vector<ForceTerm*> forceterms, "
        << "Matter *& sample, PlasmaData const *& pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    ForceTerms = forceterms;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, float accuracy, 
std::vector<ForceTerm*> forceterms, Matter *& sample, PlasmaData * pdata):
Model(filename,sample,*pdata,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::vector<ForceTerm*> forceterms, "
        << "Matter *& sample, PlasmaData const *& pdata) : "
        << "Model(sample,pdata,accuracy)\n\n");
    ForceTerms = forceterms;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, float accuracy,
std::vector<ForceTerm*> forceterms, Matter *& sample, PlasmaGrid_Data & pgrid):
Model(filename,sample,pgrid,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::vector<ForceTerm*> forceterms, "
        << "Matter *& sample, PlasmaGrid const& pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    ForceTerms = forceterms;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, float accuracy, 
std::vector<ForceTerm*> forceterms, Matter *& sample, PlasmaGrid_Data & pgrid, 
PlasmaData & pdata):
Model(filename,sample,pgrid,pdata,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::vector<ForceTerm*> forceterms, "
        << "Matter *& sample, PlasmaGrid const& pgrid) : "
        << "Model(sample,pgrid,accuracy)\n\n");
    ForceTerms = forceterms;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

void ForceModel::CreateFile(std::string filename){
    F_Debug("\tIn ForceModel::CreateFile(std::string filename)\n\n");
    FileName=filename;
    ModelDataFile.open(FileName);
    ModelDataFile << std::scientific << std::setprecision(16) << std::endl;
    ModelDataFile << "Time\tPotential\tPosition\tVelocity\tRotationFreq";

    //!< Loop over force terms and print their names
    for(auto iter = ForceTerms.begin(); iter != ForceTerms.end(); ++iter) {
        ModelDataFile << "\t" << (*iter)->PrintName();
    }
    
    ModelDataFile << "\n";
    ModelDataFile.close();
    ModelDataFile.clear();
    Print();
}

void ForceModel::Print(){
    F_Debug("\tIn ForceModel::Print()\n\n");
    ModelDataFile.open(FileName,std::ofstream::app);
    ModelDataFile << TotalTime << "\t" << Sample->get_potential() << "\t"
        << Sample->get_position() << "\t" << Sample->get_velocity() << "\t"
        << Sample->get_rotationalfreq();

    //!< Loop over force terms and print their evaluations
    for(auto iter = ForceTerms.begin(); iter != ForceTerms.end(); ++iter) {
        ModelDataFile << "\t" << (*iter)->Evaluate(Sample,Pdata);
    }

    ModelDataFile << "\n";
    ModelDataFile.close();
    ModelDataFile.clear();
}

double ForceModel::ProbeTimeStep()const{
    F_Debug( "\tIn ForceModel::ProbeTimeStep()const\n\n" );

    double timestep(0);
    threevector Acceleration = CalculateAcceleration();
    //!< For Accuracy = 1.0, requires change in velocity less than 10cm/s
    if( Acceleration.mag3() == 0 ){
        static bool runOnce = true;
        WarnOnce(runOnce,"Zero Acceleration!\ntimestep being set to unity");
        //!< Set arbitarily large time step
        timestep = 1;
    }else{
        timestep = (0.01*Accuracy)*(1.0/Acceleration.mag3());
    }

    //!< Check if the timestep should be shortened such that particles don't 
    //!< cross many grid cells in a single step
    //!< (This is often the case without this condition.)
    if( !ContinuousPlasma &&  Sample->get_velocity().mag3() != 0.0 
        && (get_dlx()*Accuracy/(2*Sample->get_velocity().mag3())) < timestep ){
        F_Debug("\ntimestep limited by grid size!");
        timestep = get_dlx()*Accuracy/(2*Sample->get_velocity().mag3());
    }
    
    //!< Check if the timestep is limited by the gyration of the particle in a
    //!< magnetic field.
    double GyromotionTimeStep = 
        Accuracy*Sample->get_velocity().mag3()*Sample->get_mass()
        *sqrt(1-(Pdata->MagneticField.getunit()*
        Sample->get_velocity().getunit()))
        /(echarge*Pdata->MagneticField.mag3());

    if( GyromotionTimeStep < timestep && GyromotionTimeStep > 0.0 ){
        std::cout << "\ntimestep limited by magnetic field (Gyromotion)\n";
        timestep = GyromotionTimeStep;
    }

    F1_Debug( "\t\tAcceleration = " << Acceleration << "\n\t\ttimestep = " 
        << timestep << "\n");
    assert(timestep == timestep);
    assert(timestep > 0);
    return timestep;
}

double ForceModel::UpdateTimeStep(){
    F_Debug( "\tIn ForceModel::UpdateTimeStep()\n\n" );
    TimeStep = ProbeTimeStep();
    
    return TimeStep;
}

void ForceModel::Force(double timestep){
    F_Debug("\tIn ForceModel::Force(double timestep)\n\n");

    //!< Make sure timestep input time is valid. Shouldn't exceed the timescale 
    //!< of the process.
    assert(timestep > 0 && timestep <= TimeStep );

    threevector Acceleration = CalculateAcceleration();
    OldTemp = Sample->get_temperature();    //!< Update OldTemp

    threevector ChangeInPosition(
        Sample->get_velocity().getx()*timestep,
        (Sample->get_velocity().gety()*timestep)/Sample->get_position().getx(),
        Sample->get_velocity().getz()*timestep);

    threevector ChangeInVelocity = Acceleration*timestep;

    //!< Assert change in absolute vel less than ten times accuracy
    assert( ChangeInVelocity.mag3() < 0.1*Accuracy );

    // Krasheninnikov, S. I. (2006). On dust spin up in uniform magnetized plasma. Physics of Plasmas, 13(11), 2004–2007.
//  double TimeOfSpinUp = Sample->get_radius()*sqrt(Pdata->mi/(Kb*Pdata->IonTemp))*Sample->get_density()/(Pdata->mi*Pdata->IonDensity);
//  TimeOfSpinUp = 1;
//  double RotationalSpeedUp = timestep*sqrt(Kb*Pdata->IonTemp/Pdata->mi)/(TimeOfSpinUp*Sample->get_radius());
//  if( Sample->get_radius() < sqrt(Kb*Pdata->IonTemp*Pdata->mi)/(echarge*Pdata->MagneticField.mag3()) ){
//          RotationalSpeedUp = (timestep*echarge*Pdata->MagneticField.mag3()/(TimeOfSpinUp*Pdata->mi));
//      F1_Debug( "\nREGIME TWO!" );
//  }
//  F1_Debug( "\nRho_{Ti} = " << sqrt(Kb*Pdata->IonTemp*Pdata->mi)/(echarge*Pdata->MagneticField.mag3()) <<
//       "\nV_{Ti} = " << sqrt(Kb*Pdata->IonTemp/Pdata->mi) << "\nOmega_{i} = " 
//      << echarge*Pdata->MagneticField.mag3()/Pdata->mi << "\ntimestep = " << timestep );
//  F1_Debug( "\nBfield.mag = " << Pdata->MagneticField.mag3() << "\nTimeOfSpinUp = " << TimeOfSpinUp
//       << "\nSpeedUp = " << RotationalSpeedUp << "\n" );


//  Introduced on 11/10/17, this is informed by over two months work on the theory and dynamics of dust rotation
//  due to ion collection. Even with the full theory, we find values of B that are too small for tokamak conditions
//  to lead to dust breakup. We need to find a way to increase the theoretical breakup speed.
//  double B = (5*sqrt(2*PI)*Pdata->IonDensity*Pdata->mi*sqrt((Kb*Pdata->IonTemp)/Pdata->mi)*pow(Sample->get_radius(),2))
//          /(2*Sample->get_mass()); 
//  double B = 1.0;
//  double RotationalSpeedUp = timestep*B*(2*(echarge*(Pdata->MagneticField.mag3())/Pdata->mi)-Sample->get_rotationalfreq());
    double RotationalSpeedUp = echarge*Pdata->MagneticField.mag3()*
        Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/Me)/
        (Sample->get_density()*Sample->get_radius())*timestep;


    Sample->update_motion(ChangeInPosition,ChangeInVelocity,RotationalSpeedUp);

    F1_Debug( "\nChangeInPosition : " << ChangeInPosition 
        << "\nChangeInVelocity : " << ChangeInVelocity << "\nAcceleration : " 
        << Acceleration << "\nTimeStep : " << TimeStep << "\n");
    F_Debug("\t"); Print();
    TotalTime += timestep;
}

void ForceModel::Force(){
    F_Debug("\tIn ForceModel::Force()\n\n");
    Force(TimeStep);
}

threevector ForceModel::CalculateAcceleration()const{
    F_Debug("\tIn ForceModel::CalculateAcceleration()const\n\n");
    threevector Accel(0.0,0.0,0.0);

    for(auto iter = ForceTerms.begin(); iter != ForceTerms.end(); ++iter) {
        Accel += (*iter)->Evaluate(Sample,Pdata);
        F1_Debug( "\n\t\t" << (*iter)->PrintName << "=" 
            << (*iter)->Evaluate(Sample,Pdata) );
    }

    F1_Debug( "\n\t\tAccel = " << Accel << "\n\n" );
    

    return Accel;
}