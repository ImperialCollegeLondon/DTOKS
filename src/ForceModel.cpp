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
    UseModel = {false,false,false,false,false,false};
    OldTemp = Sample->get_temperature();
    CreateFile("Default_Force_filename.txt");
}

ForceModel::ForceModel(std::string filename, float accuracy, 
std::array<bool,FMN> models, Matter *& sample, PlasmaData & pdata):
Model(sample,pdata,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::array<bool,3> models, Matter *& sample, "
        << " PlasmaData const *& pdata) : Model(sample,pdata,accuracy)\n\n");
    UseModel = models;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, float accuracy, 
std::array<bool,FMN> models, Matter *& sample, PlasmaData * pdata):
Model(sample,*pdata,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::array<bool,3> models, Matter *& sample, "
        << "PlasmaData const *& pdata) : Model(sample,pdata,accuracy)\n\n");
    UseModel = models;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, float accuracy,
std::array<bool,FMN> models, Matter *& sample, PlasmaGrid_Data & pgrid):
Model(sample,pgrid,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::array<bool,3> models, Matter *& sample, "
        << "PlasmaGrid const& pgrid) : Model(sample,pgrid,accuracy)\n\n");
    UseModel = models;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

ForceModel::ForceModel(std::string filename, float accuracy, 
std::array<bool,FMN> models, Matter *& sample, PlasmaGrid_Data & pgrid, 
PlasmaData & pdata):
Model(sample,pgrid,pdata,accuracy){
    F_Debug("\n\nIn ForceModel::ForceModel(std::string filename, "
        << "float accuracy, std::array<bool,3> models, Matter *& sample, "
        << "PlasmaGrid const& pgrid) : Model(sample,pgrid,accuracy)\n\n");
    UseModel = models;
    OldTemp = Sample->get_temperature();
    CreateFile(filename);
}

void ForceModel::CreateFile(std::string filename){
    F_Debug("\tIn ForceModel::CreateFile(std::string filename)\n\n");
    FileName=filename;
    ModelDataFile.open(FileName);
    ModelDataFile << std::scientific << std::setprecision(16) << std::endl;
    ModelDataFile << "Time\tPotential\tPosition\tVelocity\tRotationFreq";
    bool PrintGravity = true; // Lol
    if( UseModel[0] && PrintGravity )   ModelDataFile << "\tGravity";
    if( UseModel[1] )                   ModelDataFile << "\tCentrifugal";
    if( UseModel[2] )                   ModelDataFile << "\tLorentz";
    if( UseModel[3] )                   ModelDataFile << "\tSOMLIonDrag";
    if( UseModel[4] )                   ModelDataFile << "\tSMOMLIonDrag";
    if( UseModel[5] )                   ModelDataFile << "\tDTOKSIonDrag";
    if( UseModel[6] )                   ModelDataFile << "\tDUSTTIonDrag";
    if( UseModel[7] )                   ModelDataFile << "\tHybridIonDrag";
    if( UseModel[8] )                   ModelDataFile << "\tNeutralDrag";
    if( UseModel[9] )                   ModelDataFile << "\tRocketForce";
    
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
    bool PrintGravity = true; // Lol

    if( UseModel[0] && PrintGravity ) ModelDataFile << "\t" << Gravity(); 
    if( UseModel[1] )                 ModelDataFile << "\t" << Centrifugal();
    if( UseModel[2] )                 ModelDataFile << "\t" << LorentzForce();
    if( Pdata->IonTemp > 0.0 && Pdata->ElectronTemp > 0.0 
        && Pdata->mi > 0.0 && Pdata->IonDensity > 0.0  ){
        if( UseModel[3] )             ModelDataFile << "\t" << SOMLIonDrag();
        if( UseModel[4] )             ModelDataFile << "\t" << SMOMLIonDrag();
        if( UseModel[5] )             ModelDataFile << "\t" << DTOKSIonDrag();
        if( UseModel[6] )             ModelDataFile << "\t" << DUSTTIonDrag();
        if( UseModel[7] )             ModelDataFile << "\t" << HybridIonDrag();
    }else{
        threevector Zeros(0.0,0.0,0.0);
        if( UseModel[3] )             ModelDataFile << "\t" << Zeros;
        if( UseModel[4] )             ModelDataFile << "\t" << Zeros;
        if( UseModel[5] )             ModelDataFile << "\t" << Zeros;
        if( UseModel[6] )             ModelDataFile << "\t" << Zeros;
        if( UseModel[7] )             ModelDataFile << "\t" << Zeros;
    }
    if( Pdata->NeutralTemp > 0.0 && Pdata->NeutralDensity > 0.0 
        && Pdata->mi > 0.0 ){
        if( UseModel[8] )             ModelDataFile << "\t" << NeutralDrag();
        F1_Debug( "\n\t\tNeutralDrag = " << NeutralDrag() );
    }else if( UseModel[8] ){
        threevector Zeros(0.0,0.0,0.0);
        ModelDataFile << "\t" << Zeros;
    }
    if( UseModel[9] )                 ModelDataFile << "\t" << RocketForce();

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

    F1_Debug( "\n\t\tg = " << Gravity() );
    F1_Debug( "\n\t\tcentrifugal = " << Centrifugal() );
    F1_Debug( "\n\t\tlorentzforce = " << LorentzForce() );
    
    if( UseModel[0] ) Accel += Gravity();
    if( UseModel[1] ) Accel += Centrifugal();
    if( UseModel[2] ) Accel += LorentzForce();
    if( Pdata->IonTemp > 0.0 && Pdata->ElectronTemp > 0.0 && Pdata->mi > 0.0
        && Pdata->IonDensity > 0.0  ){
        if( UseModel[3] ) Accel += SOMLIonDrag();
        if( UseModel[4] ) Accel += SMOMLIonDrag();
        if( UseModel[5] ) Accel += DTOKSIonDrag();
        if( UseModel[6] ) Accel += DUSTTIonDrag();
        if( UseModel[7] ) Accel += HybridIonDrag();
        
        F1_Debug( "\n\t\tSOMLIonDrag = " << SOMLIonDrag() );
        F1_Debug( "\n\t\tSMOMLIonDrag = " << SMOMLIonDrag() );
        F1_Debug( "\n\t\tDTOKSIonDrag = " << DTOKSIonDrag() );
        F1_Debug( "\n\t\tDUSTTIonDrag = " << DUSTTIonDrag() );
        F1_Debug( "\n\t\tHybridIonDrag = " << HybridIonDrag() );
        
    }
    if( Pdata->NeutralTemp > 0.0 && Pdata->NeutralDensity > 0.0 
        && Pdata->mi > 0.0 ){
        if( UseModel[8] ) Accel += NeutralDrag();
        F1_Debug( "\n\t\tNeutralDrag = " << NeutralDrag() );
    }

    if( UseModel[9] ) Accel += RocketForce();
    F1_Debug( "\n\t\tRocketForce = " << RocketForce() );
    F1_Debug( "\n\t\tAccel = " << Accel << "\n\n" );
    

    return Accel;
}

threevector ForceModel::Gravity()const{
    F_Debug("\tIn ForceModel::Gravity()\n\n");
    return Pdata->Gravity;
}

threevector ForceModel::Centrifugal()const{
    F_Debug("\tIn ForceModel::Centrifugal()const\n\n");
    threevector returnval(
        Sample->get_velocity().gety()*Sample->get_velocity().gety()/
        Sample->get_position().getx(),
        -Sample->get_velocity().getx()*Sample->get_velocity().gety()/
        Sample->get_position().getx(),
        0.0);
    return returnval;
}

threevector ForceModel::LorentzForce()const{
    F_Debug("\tIn ForceModel::LorentzForce()const\n\n");
    //!< Dust grain charge to mass ratio
    double Charge = 4.0*PI*epsilon0*Sample->get_radius()*Kb*Pdata->ElectronTemp*
        Sample->get_potential()/echarge;
    double qtom = Charge/Sample->get_mass();

    // double ConvertKelvsToeV(8.621738e-5);
    // Estimate charge from potential difference
    // if( Sample->is_positive() ) 
    //      (qtom = 3.0*epsilon0*Kb*Sample->get_temperature())
    //          /(echarge*pow(Sample->get_radius(),2)*Sample->get_density());
    // else 
    //      qtom = -3.0*epsilon0*Pdata->ElectronTemp*ConvertKelvsToeV*
    //          Sample->get_potential()/
    //          (pow(Sample->get_radius(),2)*Sample->get_density());   
    //else qtom = -1000.0*3.0*epsilon0*V/(a*a*rho);//why it had Te???????
    //edo to eixa allaksei se ola ta runs sto After 28_Feb ara prepei na ksana 
    //ginoun
    // Google Translate: Here I had changed it to all the brides in After 28_Feb
    // so they have to be done again
    
    threevector returnvec = (Pdata->ElectricField+(Sample->get_velocity()^
        Pdata->MagneticField))*qtom;
    return returnvec;
}

threevector ForceModel::SOMLIonDrag()const{
    F_Debug("\tIn ForceModel::SOMLIonDrag()\n\n");
    return (Pdata->PlasmaVel-Sample->get_velocity())*
            Pdata->mi*(1.0/sqrt(Kb*Pdata->IonTemp/Pdata->mi))*
            SOMLIonFlux(Sample->get_potential())*(1.0/Sample->get_mass());
}

threevector ForceModel::SMOMLIonDrag()const{
    F_Debug("\tIn ForceModel::SMOMLIonDrag()\n\n");
    return (Pdata->PlasmaVel-Sample->get_velocity())*
            Pdata->mi*(1.0/sqrt(Kb*Pdata->IonTemp/Pdata->mi))*
            SMOMLIonFlux(Sample->get_potential())*(1.0/Sample->get_mass());
}

threevector ForceModel::HybridIonDrag()const{
    F_Debug("\tIn ForceModel::HybridIonDrag()\n\n");

    //!< Taken from :
    //!< S. A. Khrapak, A. V. Ivlev, S. K. Zhdanov, and G. E. Morfill, 
    //!< Phys. Plasmas 12, 1 (2005).
    double IonThermalVelocity = sqrt(Kb*Pdata->IonTemp/Pdata->mi);
    //!< Normalised ion flow velocity
    double u = Pdata->PlasmaVel.mag3()*(1.0/IonThermalVelocity);
    double Tau = Pdata->ElectronTemp/Pdata->IonTemp;

    if( u == 0.0 ){
        threevector Zero(0.0,0.0,0.0);
        return Zero;
    }

    double z = Sample->get_potential()*4.0*PI*epsilon0;
    double CoulombLogarithm = 17.0;     //!< Approximation of coulomb logarithm   


    double Coefficient = sqrt(2*PI)*pow(Sample->get_radius(),2.0)*
        Pdata->IonDensity*Pdata->mi*Pdata->PlasmaVel.square();
    double term1 = sqrt(PI/2.0)*erf(u/sqrt(2))*(1.0+u*u+(1.0-(1.0/(u*u)))*
        (1.0+2*Tau*z)+4*z*z*Tau*Tau*CoulombLogarithm/(u*u));
    double term2 = (1.0/u)*exp(-u*u/2.0)*(1.0+2.0*Tau*z+u*u-4*z*z*Tau*Tau*
        CoulombLogarithm);    
    
    threevector HybridDrag = Coefficient*(term1+term2)*(
        Pdata->PlasmaVel-Sample->get_velocity()).getunit();
    return HybridDrag*(1.0/Sample->get_mass());
}


threevector ForceModel::DTOKSIonDrag()const{
    F_Debug("\tIn ForceModel::DTOKSIonDrag()\n\n");
    threevector Fid(0.0,0.0,0.0);
    //!< Calculations for ion drag: Mach number, shielding length with fitting 
    //!< function and thermal scattering parameter
    //!< ION TEMPERATURE IN THIS FUNCTION IS IN ev.
    double ConvertKelvsToeV(8.621738e-5);
    threevector Mt(0.0,0.0,0.0);
    if( Pdata->IonTemp != 0 ) 
        Mt = (Pdata->PlasmaVel-Sample->get_velocity())*
            sqrt(Pdata->mi/(Kb*Pdata->IonTemp)); 
        F1_Debug("\nMt = " << Mt << "\tmi = " << Pdata->mi
            << "\nVp = " << Pdata->PlasmaVel
            << "\nVd = " << Sample->get_velocity());

        if( Pdata->IonDensity == 0 || Pdata->IonTemp == 0 
            || Pdata->ElectronTemp == 0  || Mt.mag3() == 0 ){ 

            Fid = threevector(0.0,0.0,0.0);
        }else{
            //!< Relative speed less than twice mach number, use Fortov et al 
            //!< theory with screening length 'Lambda'.
            if(Mt.mag3()<2.0){ 

                double lambda = sqrt(epsilon0/(Pdata->IonDensity*echarge*
                    exp(-Mt.mag3()*Mt.mag3()/2)*
                    (1.0/(Pdata->IonTemp*ConvertKelvsToeV))+
                    1.0/(Pdata->ElectronTemp*ConvertKelvsToeV)));
                double beta = Pdata->ElectronTemp*ConvertKelvsToeV*
                    Sample->get_radius()*fabs(Sample->get_potential())/
                    (lambda*Pdata->IonTemp*ConvertKelvsToeV);
                F1_Debug("\nlambda = " << lambda << "\nbeta = " << beta 
                    << "\nPot = " << Sample->get_potential());

                if(beta>13.0) 
                    std::cout << "nonlinear drag parameter" << std::endl;

                threevector FidS(0.0,0.0,0.0);
                if( Sample->get_potential() == 0 || beta == 0 ){
                    FidS = threevector(0.0,0.0,0.0);
                }else{
                    double Lambda = -exp(beta/2.0)*
                        Exponential_Integral_Ei(-beta/2.0); 
                    FidS = Mt*(sqrt(32*PI)/3.0*epsilon0*
                        pow(Pdata->IonTemp*ConvertKelvsToeV,2)*Lambda*
                        pow(beta,2));
                }

                threevector FidC =(Pdata->PlasmaVel-Sample->get_velocity())*4.0*
                    PI*pow(Sample->get_radius(),2)*Pdata->IonDensity*Pdata->mi*
                    sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*
                    exp(-Sample->get_potential()); 
                //for John's ion drag... I assume here and in other places in 
                //the calculation that the given potential is normalised to 
                //kTe/e Do The same but with the ion current instead of the 
                //electron current


                Fid=FidS+FidC;
                F1_Debug("\nFidS = " << FidS << "\nFidC = " << FidC);
        }else{ 
            //!< Relative speed greater than twice the mach number, 
            //!< use just plain collection area
            //double lambdadi = sqrt(epsilon0*Pdata->IonTemp*ConvertKelvsToeV)/
            //  sqrt(Pdata->IonDensity*echarge);
            //F_Debug("\nlambdadi = " << lambdadi);
            Fid = Mt*Mt.mag3()*PI*Pdata->IonTemp*ConvertKelvsToeV*
                pow(Sample->get_radius(),2)*Pdata->IonDensity*echarge;
        }
    }
    return Fid*(3/(4*PI*pow(Sample->get_radius(),3)*Sample->get_density()));
}

threevector ForceModel::DUSTTIonDrag()const{
    F_Debug("\tIn ForceModel::DUSTTIonDrag()const\n\n");
    //!< Pre-define commonly used quantities
    double Chi = Sample->get_potential(); 
    double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
    double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/
        sqrt(2.0*Kb*Pdata->IonTemp/Pdata->mi);
    double uzp = uz+sqrt(-Pdata->Z*Chi/Tau);
    double uzm = uz-sqrt(-Pdata->Z*Chi/Tau);
    double wzp = uz*uz+Pdata->Z*Chi/Tau;
    double wzm = uz*uz-Pdata->Z*Chi/Tau;

    double CircularArea = PI*Sample->get_radius()*Sample->get_radius();
    double G = (erf(uz)-2.0*uz*exp(-uz*uz))/(2.0*uz*uz*sqrt(PI));
    //!< Predefine scattering 
    double CoulombLogarithm = 17.0;     //!< Approximation of coulomb logarithm   
    double IonScatter = 2.0*CircularArea*Pdata->mi*Pdata->IonDensity*
        (Pdata->PlasmaVel-Sample->get_velocity()).mag3()*(Pdata->Z*Chi/Tau)*
        (Pdata->Z*Chi/Tau)*G*log(CoulombLogarithm);
    double IonCollect(0.0);
    if( Chi <= 0 ) 
        IonCollect = CircularArea*Pdata->mi*Pdata->IonDensity*
            sqrt(2.0*Kb*Pdata->IonTemp/Pdata->mi)*
            (Pdata->PlasmaVel-Sample->get_velocity()).mag3()*(1.0/(4.0*uz*uz))*
            ((1.0/sqrt(PI))*((1.0+2.0*uz*uz+(1.0-2.0*uz*uz)*
            sqrt(-Pdata->Z*Sample->get_potential()/Tau))*exp(-uzp*uzp)+
            (1.0+2.0*uz*uz-(1.0-2.0*uz*uz)*
            sqrt(-Pdata->Z*Sample->get_potential()/Tau))*exp(-uzm*uzm))+
            uz*(1.0+2*wzp-(1.0-2.0*wzm)/(2.0*uz*uz))*(erf(uzp)+erf(uzm)));
    else    
        IonCollect = CircularArea*Pdata->mi*Pdata->IonDensity*
            sqrt(2.0*Kb*Pdata->IonTemp/Pdata->mi)*
            (Pdata->PlasmaVel-Sample->get_velocity()).mag3()*(1.0/(2.0*uz*uz))*
            ((1.0/sqrt(PI))*(1.0+2.0*wzp)*exp(-uz*uz)+
            uz*(1.0+2*wzp-(1.0-2.0*wzm)/(2.0*uz*uz))*erf(uz));

    return (IonScatter+IonCollect)*
        ((Pdata->PlasmaVel-Sample->get_velocity()).getunit())*
        (1.0/Sample->get_mass());
}

threevector ForceModel::NeutralDrag()const{
    F_Debug("\tIn ForceModel::NeutralDrag()const\n\n");
    // Assuming OML flux of neutrals, with neutrals flowing with the
    // background plasma
    // return (Pdata->PlasmaVel-Sample->get_velocity())*Pdata->mi*sqrt(4*PI)*
    // NeutralFlux()*PI*pow(Sample->get_radius(),2)*(1.0/Sample->get_mass());

    // Assuming OML flux of neutrals, neutrals stationary with respect to dust 
    // grain and not flowing with plasma
    // return -1.0*Sample->get_velocity()*Pdata->mi*sqrt(4*PI)*NeutralFlux()*PI*
    // pow(Sample->get_radius(),2)*(1.0/Sample->get_mass());

    // (1) Pigarov A Yu, Krasheninnikov S I, Soboleva T K and Rognlien T D 2005 
    // Phys. Plasmas 12 122508
    // (2) Baines M J, Williams I P and Asebiomo A S 1965 Mon. Not. R. Astron. 
    // Soc. 130 63
    // Assuming DUSTT flux of neutrals flowing 
    //double ua = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/
    // sqrt(2.0*Kb*Pdata->NeutralTemp/Pdata->mi);

    //!< Assuming DUSTT flux of neutrals, neutrals stationary
    double ua = -Sample->get_velocity().mag3()/sqrt(2.0*Kb*Pdata->NeutralTemp
        /Pdata->mi);

    if( ua == 0.0 ){
        threevector Zeros(0.0,0.0,0.0);
        return Zeros;
    }else{
        return PI*Sample->get_radius()*Sample->get_radius()*Pdata->mi*
            Pdata->NeutralDensity*sqrt(2.0*Kb*Pdata->NeutralTemp/Pdata->mi)*
            (1.0/ua)*((1.0/sqrt(PI))*(ua+1/(2.0*ua))*exp(-ua*ua)+
            (1.0+ua*ua-1.0/(4.0*ua*ua))*erf(ua))*-1.0*Sample->get_velocity()*
            (1.0/Sample->get_mass());
    }
}

threevector ForceModel::RocketForce()const{
    F_Debug("\tIn ForceModel::RocketForce()const\n\n");
    threevector returnvec(0.0,0.0,0.0);

    if( Sample->is_liquid() ){
        double Pv_plus = 
            Sample->probe_vapourpressure(Sample->get_temperature());
        double Pv_minus = Sample->probe_vapourpressure(OldTemp);
        returnvec = Sample->get_surfacearea()*((Pv_plus-Pv_minus)/
            Sample->get_radius())*Pdata->MagneticField.getunit();
    }
    return returnvec*(1.0/Sample->get_mass());;
}
