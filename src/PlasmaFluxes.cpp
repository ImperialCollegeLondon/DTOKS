/** @file PlasmaFluxes.cpp
 *  @brief Contains functions defining the default plasma fluxes
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */
//#define MODEL_DEBUG
#include "PlasmaFluxes.h"

namespace Flux{
//!< The flux of ions onto a charged sphere following classic OML theory.
//!< The formula can be found in:
//!< P. K. Shukla and A. A. Mamun, 
//!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
//!< Equation (2.2.6) and (2.2.7).
double OMLIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    PF_Debug("\n\t\tIn OMLIonFlux:Term()");

    double IonFlux(0.0);

    if( Potential >= 0 )
        //!< P. K. Shukla and A. A. Mamun, 
        //!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
        //!< Equation (2.2.6).
        IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))*
            (1+Pdata->Z*Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
    else
        //!< P. K. Shukla and A. A. Mamun, 
        //!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
        //!< Equation (2.2.7).
        IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
                *exp(Potential*(Pdata->ElectronTemp/Pdata->IonTemp));

    //!< Sanity check sensible return value
    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;

    std::string Warning = "\nError in OMLIonFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nOMLIonFlux() Return value: IonFlux = " << IonFlux);
    PF_Debug("\nPdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))= " 
        << Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi)));
    PF_Debug("\nPotential*(Pdata->ElectronTemp/Pdata->IonTemp)= " 
        << Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of ions onto a negatively charged sphere following MOML theory.
//!< The formula can be found in:
//!< C. T. N. Willis, M. Coppins, M. Bacharis, and J. E. Allen, 
//!< Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 85, (2012).
//!< Equation (3).
//!< For positively charged dust, we resort to OML
double MOMLIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    PF_Debug("\n\t\tIn MOMLIonFlux:Term()");

    double IonFlux(0.0);
    double HeatCapacityRatio(5.0/3.0);

    if( Potential >= 0 ){
        //!< C. T. N. Willis, M. Coppins, M. Bacharis, and J. E. Allen, 
        //!< Phys. Rev. E - Stat. Nonlinear, Soft Matter Phys. 85, (2012).
        //!< Equation (3).
        IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
            *(1.0-(Pdata->ElectronTemp/Pdata->IonTemp)*
            (-Pdata->Z*Potential-0.5*log((2.0*PI*Me/Pdata->mi)*
            (1.0+HeatCapacityRatio*Pdata->IonTemp/Pdata->ElectronTemp))));
    }else
        //!< For positively charged dust, we resort to OML
        IonFlux = OMLIonFlux(Sample,Pdata,Potential);
    //!< Sanity check sensible return value
    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in MOMLIonFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nMOMLIonFlux() Return value: IonFlux = " << IonFlux);
    PF_Debug("\nPdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))= " 
        << Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi)));
    PF_Debug("\n(Pdata->ElectronTemp/Pdata->IonTemp)= " 
        << (Pdata->ElectronTemp/Pdata->IonTemp));
    PF_Debug("\nPdata->Z*Potential = " << Pdata->Z*Potential);
    PF_Debug("\nlog( arg ) = " << log((2.0*PI*Me/Pdata->mi)*
        (1.0+HeatCapacityRatio*Pdata->IonTemp/Pdata->ElectronTemp)));
    PF_Debug("\narg = " << (2.0*PI*Me/Pdata->mi)*
        (1.0+HeatCapacityRatio*Pdata->IonTemp/Pdata->ElectronTemp));
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of ions onto a negatively charged sphere following SOML theory.
//!< For negativel charged dust, the formula can be found in:
//!< D. Thomas, Theory and Simulation of the Charging of Dust in Plasmas, 2016.
//!< Equation (2.132), (2.133) and (2.134)
//!< P. K. Shukla and A. A. Mamun, 
//!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
//!< Equation (2.2.9).
//!< R. D. Smirnov, A. Y. Pigarov, M. Rosenberg, S. I. Krasheninnikov, and 
//!< D. a Mendis, Plasma Phys. Control. Fusion 49, 347 (2007).
//1< Equation (1).

//!< For the case of positively charged dust, the formula can be found in:
//!< R. D. Smirnov, A. Y. Pigarov, M. Rosenberg, S. I. Krasheninnikov, and 
//!< D. a Mendis, Plasma Phys. Control. Fusion 49, 347 (2007).
//!< Equation (2).

//!< For the case of no flow, to avoid dividing by zero we return OMLIonFlux.
double SOMLIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    PF_Debug( "\n\t\tIn SOMLIonFlux:Term()\n\n");
    
    //!< Calculate the Ion Thermal Velocity
    double IonThermalVelocity = sqrt((2.0*Kb*Pdata->IonTemp)/Pdata->mi);
    //!< uz is the relative velocity normalised to ion thermal speed
    double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()
        /IonThermalVelocity;
    //!< Tau is the ion to electron temperature ratio
    double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
    double IonFlux(0.0);
    if( uz == 0.0 ){
        //!< For no flow case, avoid dividing by zero return OMLIonFlux.
        return Flux::OMLIonFlux(Sample,Pdata,Potential);
    }else{

        if( Potential >= 0.0 ){
            //!< For negativel charged dust, the formula can be found in:
            //!< D. Thomas, Theory and Simulation of the Charging of Dust in 
            //!< Plasmas, 2016.
            //!< Equation (2.132), (2.133) and (2.134)
            //!< P. K. Shukla and A. A. Mamun, 
            //!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
            //!< Equation (2.2.9).
            double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+
                exp(-uz*uz)/2.0;
            double s2 = sqrt(PI)*erf(uz)/(2.0*uz);
            IonFlux = Pdata->IonDensity*
                (IonThermalVelocity/sqrt(4.0*PI))*
                (s1+s2*Pdata->Z*Potential/Tau);
        }else{
            //!< For the case of positively charged dust:
            //!< R. D. Smirnov, A. Y. Pigarov, M. Rosenberg, 
            //!< S. I. Krasheninnikov, and D. a Mendis, Plasma Phys. Control. 
            //!< Fusion 49, 347 (2007).
            //!< Equation (2).
            double uzp = uz+sqrt(-Pdata->Z*Potential/Tau);
            double uzm = uz-sqrt(-Pdata->Z*Potential/Tau);
            IonFlux = Pdata->IonDensity*IonThermalVelocity*
                (1.0/(4.0*uz))*
                ((1.0+2.0*(uz*uz+Pdata->Z*Potential/Tau))*
                (erf(uzp)+erf(uzm))+(2.0/sqrt(PI))*
                (uzp*exp(-uzm*uzm)+uzm*exp(-uzp*uzp)));
        }
    }
    PF_Debug("\nPdata->IonDensity = " << Pdata->IonDensity);
    PF_Debug("\nIonThermalVelocity = " << IonThermalVelocity);
    PF_Debug("\nuz = " << uz);
    PF_Debug("\nTau = " << Tau);
    PF_Debug("\nPdata->Z*Potential = " << Pdata->Z*Potential);
    //!< Sanity check sensible return value
    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in SOMLIonFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nSOMLIonFlux() Return value: IonFlux = " << IonFlux);
    PF_Debug("\nPdata->IonDensity = " << Pdata->IonDensity);
    PF_Debug("\nIonThermalVelocity = " << IonThermalVelocity);
    PF_Debug("\nuz = " << uz);
    PF_Debug("\nTau = " << Tau);
    PF_Debug("\nPdata->Z*Potential = " << Pdata->Z*Potential);
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of ions onto a negatively charged sphere following SMOML theory.
//!< For negatively charged dust, the formula can be found in:
//!< D. Thomas, Theory and Simulation of the Charging of Dust in Plasmas, 2016.
//!< Equation (2.133), (2.134) and (2.139)
//!< C. M. Bray, (2014).
//!< Equation (2.90)

//!< For no flow case, avoid dividing by zero return MOMLIonFlux.
//!< For Positive dust case, do SOML
double SMOMLIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    PF_Debug( "\n\t\tIn SMOMLIonFlux:Term()\n\n");
    double MassRatio = Pdata->mi/Me;
    //!< Calculate the Ion Thermal Velocity
    double IonThermalVelocity = sqrt((2.0*Kb*Pdata->IonTemp)/Pdata->mi);
    //!< uz is the relative velocity normalised to ion thermal speed
    double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()
        /IonThermalVelocity;
    //!< Tau is the ion to electron temperature ratio
    double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
    double IonFlux(0.0);
    if( uz == 0.0 ){ 
        //!< For no flow case, avoid dividing by zero return MOMLIonFlux.
        IonFlux = Flux::MOMLIonFlux(Sample,Pdata,Potential);
    }else{
        if( Potential >= 0.0 ){
            //!< For negatively charged dust, the formula can be found in:
            //!< D. Thomas, Theory and Simulation of the Charging of Dust in 
            //!< Plasmas, 2016.
            //!< Equation (2.133), (2.133) and (2.139)
            double HeatCapacityRatio = 5.0/3.0;
            double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+
                exp(-uz*uz)/2.0;
            double s2 = sqrt(PI)*erf(uz)/(2.0*uz);
            IonFlux = Pdata->IonDensity*(IonThermalVelocity/sqrt(4.0*PI))*
                (s1-(s2/Tau)*(-Potential*Pdata->Z-
                0.5*log(2.0*PI*(1.0+HeatCapacityRatio*Tau)/MassRatio)));
        }else{ 
            //!< For Positive dust, resort to SOML
            IonFlux = Flux::SOMLIonFlux(Sample,Pdata,Potential);
        }
    }
    //!< Sanity check sensible return value
    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in SMOMLIonFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nSMOMLIonFlux() Return value: IonFlux = " << IonFlux);
    PF_Debug("\nPdata->IonDensity = " << Pdata->IonDensity);
    PF_Debug("\nIonThermalVelocity = " << IonThermalVelocity);
    PF_Debug("\nuz = " << uz);
    PF_Debug("\nMassRatio = " << MassRatio);
    PF_Debug("\nTau = " << Tau);
    PF_Debug("\nPdata->Z*Potential = " << Pdata->Z*Potential);
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of electrons on a negatively charged sphere following PHL.
//!< For negatively charged dust, the formula can be found in:
//!< L. Patacchini, I. H. Hutchinson, and G. Lapenta, Phys. Plasmas 14, (2007).
//!< Equations (8) to equation (15) inclusive

//!< For Positive dust case, do OML
double PHLElectronFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    PF_Debug( "\n\t\tIn PHLElectronFlux:Term()\n\n");

    //!< Tau is the ion to electron temperature ratio
    double Tau = Pdata->ElectronTemp/Pdata->IonTemp;
    //!< Beta is the dust radius to ion gyro-radius ratio
    double Beta = Sample->get_radius()
            /(sqrt(PI*Pdata->ElectronTemp*Me)/(2.0*echarge*echarge*
            Pdata->MagneticField*Pdata->MagneticField));
    double MassRatio = Pdata->mi/Me;

    if( Beta/MassRatio > 0.01 ){
        //!< PHL give a limited range for their model
        static bool runOnce = true;
        std::string Warning = "Beta/MassRatio > 0.01 in solvePHL! Model may ";
        Warning += "not be valid in this range! see Fig 11. of  L. Patacchini,";
        Warning += " I. H. Hutchinson, and G. Lapenta, ";
        Warning += "Phys. Plasmas 14, (2007).";
        WarnOnce(runOnce,Warning);
    }
    double AtomicNumber = Pdata->Z; 
    //!< Calculate the electron debye length of the plasma
    double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)
        /(Pdata->ElectronDensity*pow(echarge,2)));

    //!< Calculate the result of equation (8)
    double z = Beta/(1.0+Beta);
    double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+
        1.150*z*z*z*z*z;
    //!< Calculate the result of equation (15)
    double eta = (Potential/Beta)*(1.0+(Beta/4.0)*
        (1-exp(-4.0/(DebyeLength*Beta))));

    double w(1.0);
    if( Beta == 0.0 || eta == -1.0 ){
        w = 1.0;
    }else if( std::isnan(eta) ){
        PF_Debug("\nWarning! w being set to 1.0 "
            << "(Assuming high B field limit) but Phi/Beta is nan while "
            << "Beta != 0.");
        w = 1.0;
    }else{
        //!< Calculate the result of equation (11)
        w = eta/(1+eta);
    } 
    
    //!< Calculate the result of equation (12)
    double A = 0.678*w+1.543*w*w-1.212*w*w*w;
    double ElecFlux(0.0);
    if( Potential >= 0.0 ){ 
        //!< For negatively charged dust, the formula can be found in:
        //!< Solve equation (13)
        ElecFlux = Pdata->ElectronDensity*
            sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(A+(1.0-A)*i_star)*
            exp(-Potential);
    }else{ //!< For positive dust, do OML
        ElecFlux = Flux::OMLElectronFlux(Pdata,Potential);
    }
    //!< Sanity check sensible return value
    if(ElecFlux >= Underflows::Flux && ElecFlux != INFINITY 
        && ElecFlux == ElecFlux && ElecFlux < Overflows::Flux ){
        return ElecFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in PHLElectronFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

//!< The flux of Ions following the original model of DTOKS
//!< The equations defining the flux can be found in:
//!< M. Bacharis, M. Coppins, and J. E. Allen, Phys. Plasmas 17, (2010).
//!< Equations (8), (9) and (10)
double DTOKSIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    PF_Debug("\n\t\tIn DTOKSIonFlux:Term()\n\n");
    double IonFlux=0;

    if( Sample->is_positive() ) 
        IonFlux = Flux::DTOKSElectronFlux(Pdata, Potential); 
    else    
        IonFlux = Flux::DTOKSElectronFlux(Pdata, Potential)*
            (1-Sample->get_deltatot());

    //!< Sanity check sensible return value
    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in DTOKSIonFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nDTOKSIonFlux() Return value: IonFlux = " << IonFlux);
    PF_Debug("\nSample->is_positive() = " << Sample->is_positive());
    PF_Debug("\nSample->get_deltatot() = " << Sample->get_deltatot());
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of Electrons following the original model of DTOKS
//!< The equations defining the flux can be found in:
//!< M. Bacharis, M. Coppins, and J. E. Allen, Phys. Plasmas 17, (2010).
//!< Equations (8), (9) and (10)
double DTOKSElectronFlux(const std::shared_ptr<PlasmaData> Pdata, 
        const double Potential){
    PF_Debug("\n\t\tIn DTOKSElectronFlux:Term()\n\n");
    double ElecFlux(0.0);

    ElecFlux = Pdata->ElectronDensity*exp(-Potential)*
        sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));

    if(ElecFlux >= Underflows::Flux && ElecFlux != INFINITY 
        && ElecFlux == ElecFlux && ElecFlux < Overflows::Flux ){
        return ElecFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in DTOKSElectronFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nDTOKSElectronFlux() Return value: ElecFlux = " << ElecFlux);
    PF_Debug("\nElectronThermalVel = " << sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me)));
    PF_Debug("\nPotential = " << Potential);
    PF_Debug("\nDensity = " << Pdata->ElectronDensity);
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of electrons onto a charged sphere following classic OHML theory.
//!< The formula can be found in:
//!< P. K. Shukla and A. A. Mamun, 
//!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
//!< Equation (2.2.6) and (2.2.7).
double OMLElectronFlux(const std::shared_ptr<PlasmaData> Pdata, 
        const double Potential){
    PF_Debug("\n\t\tIn OMLElectronFlux:Term()\n\n");

    double ElecFlux(0.0);
    if( Potential < 0.0 )
        //!< P. K. Shukla and A. A. Mamun, 
        //!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
        //!< Equation (2.2.6).
        ElecFlux = Pdata->ElectronDensity*
            sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me))*(1-Potential);
    else
        //!< P. K. Shukla and A. A. Mamun, 
        //!< Introduction to Dusty Plasma Physics (CRC Press, 2015).
        //!< Equation (2.2.7).
        ElecFlux = Pdata->ElectronDensity*exp(-Potential)*
            sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
    
    //!< Sanity check sensible return value
    if(ElecFlux >= Underflows::Flux && ElecFlux != INFINITY 
        && ElecFlux == ElecFlux && ElecFlux < Overflows::Flux ){
        return ElecFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in OMLElectronFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nOMLElectronFlux() Return value: ElecFlux = " << ElecFlux);
    PF_Debug("\nElectronThermalVel = " << sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me)));
    PF_Debug("\nElectronTemp = " << Pdata->ElectronTemp);
    PF_Debug("\nPotential = " << Potential);
    PF_Debug("\nDensity = " << Pdata->ElectronDensity);
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of neutrals onto a sphere for a stationary maxwellian
//!< For an equation, see the following:
//!< A. Y. Pigarov, S. I. Krasheninnikov, T. K. Soboleva, and T. D. Rognlien, 
//!< Phys. Plasmas 12, 1 (2005).
//!< Pg. 7, top right hand side of page
double NeutralFlux(const std::shared_ptr<PlasmaData> Pdata){
    PF_Debug("\n\t\tIn NeutralFlux:Term()\n\n");

    double NeutFlux = Pdata->NeutralDensity*
            sqrt(Kb*Pdata->NeutralTemp/(2*PI*Pdata->mi));

    //!< Sanity check sensible return value
    if(NeutFlux >= Underflows::Flux && NeutFlux != INFINITY 
        && NeutFlux == NeutFlux && NeutFlux < Overflows::Flux ){
        return NeutFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in NeutralFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nNeutralFlux() Return value: NeutFlux = " << NeutFlux);
    PF_Debug("\nNeutralThermalVel = " << sqrt(Kb*Pdata->NeutralTemp/(2*PI*Pdata->mi)));
    PF_Debug("\nDensity = " << Pdata->NeutralDensity);
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The flux of particles from the sphere for a stationary maxwellian, see:
//!< S. I. Krasheninnikov, R. D. Smirnov, and D. L. Rudakov, 
//!< Plasma Phys. Control. Fusion 53, 083001 (2011).
//!< Equation (35)
//!< Note here that the ambient vapour pressure has also been accounted for
double EvaporationFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    PF_Debug("\n\t\tIn EvaporationFlux:Term()\n\n");
    double AmbientPressure = Pdata->NeutralDensity*Kb*Pdata->NeutralTemp;

    double StickCoeff = 1.0;
    double EvapFlux = (StickCoeff*Sample->get_surfacearea()*AvNo*
        (Sample->get_vapourpressure()-AmbientPressure))/
        sqrt(2*PI*Sample->get_atomicmass()*R*DustTemperature);

    //!< Sanity check sensible return value
    //!< Missing "EvapFlux >= 0.0 &&" permits dust gaining mass from 
    //!< condensation of nearby cold nuclei.
    if( EvapFlux != INFINITY && EvapFlux == EvapFlux ){
        return EvapFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in EvaporationFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nEvaporationFlux() Return value: EvapFlux = " << EvapFlux);
    PF_Debug("\nAmbientPressure = " << AmbientPressure);
    PF_Debug("\nStickCoeff = " << StickCoeff);
    PF_Debug("\nSample->get_vapourpressure() = " << Sample->get_vapourpressure());
    PF_Debug("\nSample->get_surfacearea() = " << Sample->get_surfacearea());
    PF_Debug("\nAmbientPressure = " << AmbientPressure;   ) 
    PF_Debug("\nDustTemperature = " << DustTemperature;) 
    PF_Debug("\nDenominator = "
        << sqrt(2*PI*Sample->get_atomicmass()*R*DustTemperature)); 
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The thermionic electron emission yield coefficient of DTOKS
//!< The equations defining the yield can be found in:
//!< M. Bacharis, M. Coppins, and J. E. Allen, Phys. Plasmas 17, (2010).
//!< Equations (5) and (7)
double DeltaTherm(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata){
    C_Debug("\n\t\tIn DeltaTherm:Term()\n\n");

    double dtherm = (Richardson*Sample->get_temperature()*
        Sample->get_temperature()*exp(-(Sample->get_workfunction()*echarge)/
        (Kb*Sample->get_temperature())))/
        (echarge*OMLElectronFlux(Pdata, Sample->get_potential()));

    //!< Sanity check sensible return value
    if(dtherm >= 0.0 && dtherm == dtherm && dtherm != INFINITY ){
        return dtherm;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in DeltaTherm()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nDeltaTherm() Return value: dtherm = " << dtherm);
    PF_Debug("\nRichardson = " << Richardson);
    PF_Debug("\nSample->get_temperature() = " << Sample->get_temperature());
    PF_Debug("\nSample->get_workfunction() = " << Sample->get_workfunction());
    PF_Debug("\nSample->get_temperature() = " << Sample->get_temperature());
    PF_Debug("\nexp( arg ) = " << exp((Sample->get_workfunction()*echarge))/
        (Kb*Sample->get_temperature()));
    PF_Debug("\narg = " << (Sample->get_workfunction()*echarge)/
        (Kb*Sample->get_temperature()));
    PF_Debug("\necharge*OMLElectronFlux() = " 
        << (echarge*OMLElectronFlux(Pdata, Sample->get_potential()))); 
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The thermionic electron emission flux of DTOKS
//!< The equations defining the flux can be found in:
//!< R. D. Smirnov, A. Y. Pigarov, M. Rosenberg, S. I. Krasheninnikov, and 
//!< D. a Mendis, Plasma Phys. Control. Fusion 49, 347 (2007).
//!< Equation (3) with zero potential
double ThermFlux(const Matter* Sample){
    C_Debug("\n\t\tIn DeltaTherm:Term()\n\n");

    double ThermFlux = (Richardson*Sample->get_temperature()*
        Sample->get_temperature()*exp(-(Sample->get_workfunction()*echarge)/
        (Kb*Sample->get_temperature())))/echarge;
    //!< Sanity check sensible return value
    if(ThermFlux >= Underflows::Flux && ThermFlux == ThermFlux 
        && ThermFlux != INFINITY && ThermFlux < Overflows::Flux ){
        return ThermFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in ThermFlux()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nThermFlux() Return value: ThermFlux = " << ThermFlux);
    PF_Debug("\nRichardson = " << Richardson);
    PF_Debug("\nSample->get_temperature() = " << Sample->get_temperature());
    PF_Debug("\nSample->get_workfunction() = " << Sample->get_workfunction());
    PF_Debug("\nSample->get_temperature() = " << Sample->get_temperature());
    PF_Debug("\nexp( arg ) = " << exp(-(Sample->get_workfunction()*echarge)/
        (Kb*Sample->get_temperature())));
    PF_Debug("\narg = " << -(Sample->get_workfunction()*echarge)/
        (Kb*Sample->get_temperature()));
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The thermionic electron emission flux with Schottky correction
//!< The equations defining the flux can be found in:
//!< R. D. Smirnov, A. Y. Pigarov, M. Rosenberg, S. I. Krasheninnikov, and 
//!< D. a Mendis, Plasma Phys. Control. Fusion 49, 347 (2007).
//!< Equation (3) with zero potential
double ThermFluxSchottky(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential ){
    C_Debug("\n\t\tIn DeltaTherm:Term()\n\n");
    //!< Returns the flux of electrons due to Thermionic emission
    //!< Following the Richard-Dushmann formula with Schottky Correction.
    //!< Negative Dust grains have the normal Richard-Dushmann form, dependant 
    //!< on work funciton
    //!< Positive Dust grains are the same flux multiplied by (1-phi)
    //!< e^(-e phi/(kb Td))
    double ThermFlux(0.0);
    if( Potential >= 0.0 ){
        ThermFlux = Richardson*Sample->get_temperature()*
            Sample->get_temperature()*exp(-echarge*Sample->get_workfunction()/
            (Kb*Sample->get_temperature()))/echarge;
    }else{

        ThermFlux = Richardson*Sample->get_temperature()*
            Sample->get_temperature()*(1.0-Potential
            *(Pdata->ElectronTemp/Sample->get_temperature()))
            *exp((-echarge*Sample->get_workfunction()-
            Potential*Kb*Pdata->ElectronTemp)/(Kb*Sample->get_temperature()))/
            echarge;
    }
    //!< Sanity check sensible return value
    if(ThermFlux >= Underflows::Flux && ThermFlux == ThermFlux 
        && ThermFlux != INFINITY && ThermFlux < Overflows::Flux ){
        return ThermFlux;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in ThermFluxSchottky()!";
    Warning += " Return value badly specified\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nThermFluxSchottky() Return value: ThermFlux = " << ThermFlux);
    PF_Debug("\nRichardson = " << Richardson);
    PF_Debug("\nSample->get_temperature() = " << Sample->get_temperature());
    PF_Debug("\nSample->get_workfunction() = " << Sample->get_workfunction());
    PF_Debug("\nSample->get_temperature() = " << Sample->get_temperature());
    PF_Debug("\nPotential = " << Potential);
    PF_Debug("\nTe/Td = " << (Pdata->ElectronTemp/Sample->get_temperature()));
    PF_Debug("\nexp( arg ) = " << exp((-echarge*Sample->get_workfunction()-
            Potential*Kb*Pdata->ElectronTemp)/(Kb*Sample->get_temperature())));
    PF_Debug("\narg = " << (-echarge*Sample->get_workfunction()-
            Potential*Kb*Pdata->ElectronTemp)/(Kb*Sample->get_temperature()));
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

//!< The secondary electron emission yield coefficient of DTOKS
//!< The equations defining the yield can be found in:
//!< M. Bacharis, M. Coppins, and J. E. Allen, Phys. Plasmas 17, (2010).
//!< Equation (1-4)
double DeltaSec(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata){
    C_Debug("\n\t\tIn DeltaSec:Term()\n\n");
    double ConvertKtoev(8.6173303e-5);
    double DeltaSec = sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem());
    //!< Sanity check sensible return value
    if(DeltaSec >= 0.0 && DeltaSec == DeltaSec && DeltaSec != INFINITY ){
        return DeltaSec;
    }
    //!< If return value is not well defined, print error and return 0.
    static bool runOnce = true;
    std::string Warning = "\nError in DeltaSec()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    PF_Debug("\n\nDeltaSec() Return value: DeltaSec = " << DeltaSec);
    PF_Debug("\nSample->get_elem() = " << Sample->get_elem());
    PF_Debug("\nPdata->ElectronTemp*ConvertKtoev = " 
        << Pdata->ElectronTemp*ConvertKtoev);
    PF_Debug("\nsec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem()) = " 
        << sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem()));
    PF_Debug("\nReturning zero!\n");
    return 0.0;
}

}
