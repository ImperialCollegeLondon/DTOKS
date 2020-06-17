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
    CreateFile(filename);
}

void ForceModel::CreateFile(std::string filename){
    F_Debug("\tIn ForceModel::CreateFile(std::string filename)\n\n");
    FileName=filename;
    ModelDataFile.open(FileName);
    ModelDataFile << std::scientific << std::setprecision(16) << std::endl;
    ModelDataFile << "Time\tPosition\tVelocity\tRotationFreq";

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
    ModelDataFile << TotalTime << "\t" << Sample->get_position() << "\t" 
        << Sample->get_velocity() << "\t" << Sample->get_rotationalfreq();

    //!< Loop over force terms and print their evaluations
    for(auto iter = ForceTerms.begin(); iter != ForceTerms.end(); ++iter) {
        ModelDataFile << "\t" 
            << (*iter)->Evaluate(Sample,Pdata,Sample->get_velocity());
    }

    ModelDataFile << "\n";
    ModelDataFile.close();
    ModelDataFile.clear();
}

double ForceModel::ProbeTimeStep()const{
    F_Debug( "\tIn ForceModel::ProbeTimeStep()const\n\n" );

    double timestep(0);
    threevector Acceleration 
        = CalculateAcceleration(Sample->get_position(),Sample->get_velocity());

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
        Accuracy*Sample->get_mass()/(echarge*Pdata->MagneticField.mag3());
//    double GyromotionTimeStep = 
//        Accuracy*Sample->get_velocity().mag3()*Sample->get_mass()
//        *sqrt(1-(Pdata->MagneticField.getunit()*
//        Sample->get_velocity().getunit()))
//        /(echarge*Pdata->MagneticField.mag3());

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

    //!< Code for Euler time step, this was the original DTOKSU method
    threevector Acceleration 
        = CalculateAcceleration(Sample->get_position(),Sample->get_velocity());
    threevector ChangeInPosition(
        Sample->get_velocity().getx()*timestep,
        (Sample->get_velocity().gety()*timestep)/Sample->get_position().getx(),
        Sample->get_velocity().getz()*timestep);
    threevector ChangeInVelocity = Acceleration*timestep;

    //!< Code for 4th order Runge Kutta time step, higher accuracy and stability
    threevector xi = Sample->get_position();
    threevector vi = Sample->get_velocity();
    RungeKutta4(xi,vi,timestep);
    //std::cout << "\n\nEuler dx = " << ChangeInPosition;
    //std::cout << "\nEuler dv = " << ChangeInVelocity;
    //std::cout << "\nRK4 dx = " << xi-Sample->get_position();
    //std::cout << "\nRK4 dv = " << vi-Sample->get_velocity(); std::cin.get();
    
    //!< Assert change in absolute vel less than ten times accuracy
//    assert( ChangeInVelocity.mag3() < 0.1*Accuracy );
    assert( ChangeInVelocity.mag3() < Accuracy );

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
    F_Debug("\t"); 
	if( PrintSteps >= PrintInterval ){
	    Print();
		PrintSteps = 1;
	}else{
        PrintSteps ++;
	}
    TotalTime += timestep;
}

void ForceModel::Force(){
    F_Debug("\tIn ForceModel::Force()\n\n");
    Force(TimeStep);
}

threevector ForceModel::CalculateAcceleration(threevector position,
        threevector velocity)const{
    F_Debug("\tIn ForceModel::CalculateAcceleration(threevector position, "
        << "threevector velocity)const\n\n");

    //!< First term which is always present is the centrifugal acceleration due
    //!< to the cylindrical coordinate system.
    threevector Accel(
        velocity.gety()*velocity.gety()/position.getx(),
        -velocity.getx()*velocity.gety()/position.getx(),
        0.0);

    //!< Sum all other force terms for the velocity given.
    for(auto iter = ForceTerms.begin(); iter != ForceTerms.end(); ++iter) {
        Accel += (*iter)->Evaluate(Sample,Pdata,velocity);
        F1_Debug( "\n\t\t" << (*iter)->PrintName() << " = " 
            << (*iter)->Evaluate(Sample,Pdata,velocity) );
    }

    F1_Debug( "\n\t\tAccel = " << Accel << "\n\n" );

    return Accel;
}

void ForceModel::RungeKutta4(threevector &xf, threevector &vf, 
        double timestep)const{

    threevector xi = Sample->get_position();
    threevector vi = Sample->get_velocity();

    threevector k1x(
        vi.getx()*timestep,
        (vi.gety()*timestep)/xi.getx(), 
        vi.getz()*timestep);
    threevector k1v = timestep*CalculateAcceleration(xi,vi);

    threevector k2x(
        (vi.getx()+k1v.getx()/2.0)*timestep,
        ((vi.gety()+k1v.gety())*timestep)/(xi.getx()+k1x.getx()/2.0), 
        (vi.getz()+k1v.getz()/2.0)*timestep);
    threevector k2v 
        = timestep*CalculateAcceleration(xi+k1x*(1.0/2.0),vi+k1v*(1.0/2.0));

    threevector k3x(
        (vi.getx()+k2v.getx()/2.0)*timestep,
        ((vi.gety()+k2v.gety())*timestep)/(xi.getx()+k2x.getx()/2.0), 
        (vi.getz()+k2v.getz()/2.0)*timestep);
    threevector k3v 
        = timestep*CalculateAcceleration(xi+k2x*(1.0/2.0),vi+k2v*(1.0/2.0));

    threevector k4x(
        (vi.getx()+k3v.getx()/2.0)*timestep,
        ((vi.gety()+k3v.gety())*timestep)/(xi.getx()+k3x.getx()/2.0), 
        (vi.getz()+k3v.getz()/2.0)*timestep);
    threevector k4v = timestep*CalculateAcceleration(xi+k3x,vi+k3v);

    xf = xi + (k1x + 2.0*(k2x+k3x) + k4x)*(1.0/6.0);
    vf = vi + (k1v + 2.0*(k2v+k3v) + k4v)*(1.0/6.0);
}
