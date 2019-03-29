/** @file PlasmaFluxes.cpp
 *  @brief Contains functions defining the default plasma fluxes
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */
//#define MODEL_DEBUG
#include "PlasmaFluxes.h"

namespace Flux{

double OMLIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug("\n\t\tIn OMLIonFlux:Term()");

    double IonFlux(0.0);

    if( Potential >= 0 ) 
        IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
                *(1+Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
    else    
        IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
                *exp(Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
    std::cout << "\nPdata->IonDensity = " << Pdata->IonDensity;
    std::cout << "\nPdata->IonTemp = " << Pdata->IonTemp;
    std::cout << "\nPdata->mi = " << Pdata->mi;
    std::cout << "\nPdata->ElectronTemp = " << Pdata->ElectronTemp;
    std::cout << "\nPotential = " << Potential;

    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in OMLIonFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double SOMLIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug( "\n\t\tIn SOMLIonFlux:Term()\n\n");
    
    double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
    double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()
        /IonThermalVelocity;
    double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
    double IonFlux(0.0);
    if( uz == 0.0 ){
        return Flux::OMLIonFlux(Sample,Pdata,Potential);
    }else{

        if( Potential >= 0.0 ){
            double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+
                exp(-uz*uz)/2.0;
            double s2 = sqrt(PI)*erf(uz)/(2.0*uz);

            IonFlux = Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*
                    (s1+s2*Pdata->Z*Potential/Tau);
        }else{
            double uzp = uz+sqrt(-Pdata->Z*Potential/Tau);
            double uzm = uz-sqrt(-Pdata->Z*Potential/Tau);

            IonFlux = Pdata->IonDensity*IonThermalVelocity*
                (1.0/(8.0*sqrt(2.0)*uz))*
                ((1.0+2.0*(uz*uz+Pdata->Z*Potential/Tau))*
                (erf(uzp)+erf(uzm))+(2.0/sqrt(PI))*
                (uzp*exp(-uzm*uzm)+uzm*exp(-uzp*uzp)));
        }
    }
    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in SOMLIonFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double SMOMLIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug( "\n\t\tIn SMOMLIonFlux:Term()\n\n");
    double HeatCapacityRatio = 5.0/3.0;
    double MassRatio = Pdata->mi/Me;

    double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
    double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
    double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()
        /IonThermalVelocity;
    double IonFlux(0.0);
    if( uz == 0.0 ){ //!< Handle zero-relative velocity case
        IonFlux = Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0))*
            (1.0+(1.0/Tau)*(Potential*Pdata->Z-
            0.5*log(2.0*PI*(1.0+HeatCapacityRatio*Tau)/MassRatio)));
    }else{
        if( Potential >= 0.0 ){ //!< For negative dust, do SMOML
            double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+
                exp(-uz*uz)/2.0;
            double s2 = sqrt(PI)*erf(uz)/(2.0*uz);
            IonFlux = Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*
                (s1+(s2/Tau)*(Potential*Pdata->Z-
                0.5*log(2.0*PI*(1.0+HeatCapacityRatio*Tau)/MassRatio)));
        }else{  //!< For Positive dust, do SOML
            double uzp = uz+sqrt(-Pdata->Z*Potential/Tau);
            double uzm = uz-sqrt(-Pdata->Z*Potential/Tau);
            IonFlux = Pdata->IonDensity*IonThermalVelocity*
                (1.0/(8.0*sqrt(2.0)*uz))*
                ((1.0+2.0*(uz*uz+Pdata->Z*Potential/Tau))*(erf(uzp)+erf(uzm))
                +(2.0/sqrt(PI))*(uzp*exp(-uzm*uzm)+uzm*exp(-uzp*uzp)));
        }
    }
    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in SMOMLIonFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double PHLElectronFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug( "\n\t\tIn PHLElectronFlux:Term()\n\n");

    double Tau = Pdata->ElectronTemp/Pdata->IonTemp;
    double Beta = Sample->get_radius()
            /(sqrt(PI*Pdata->ElectronTemp*Me)/(2.0*echarge*echarge*
            Pdata->MagneticField*Pdata->MagneticField));
    double MassRatio = Pdata->mi/Me;

    if( Beta/MassRatio > 0.01 ){
        static bool runOnce = true;
        std::string Warning = "Beta/MassRatio > 0.01 in solvePHL! Model may ";
        Warning += "not be valid in this range! see Fig 11. of  L. Patacchini,";
        Warning += " I. H. Hutchinson, and G. Lapenta, ";
        Warning += "Phys. Plasmas 14, (2007).";
        WarnOnce(runOnce,Warning);
    }
    double AtomicNumber = Pdata->Z; 
    double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)
        /(Pdata->ElectronDensity*pow(echarge,2)));

    double z = Beta/(1.0+Beta);
    double i_star = 1.0-0.0946*z-0.305*z*z+0.950*z*z*z-2.2*z*z*z*z+
        1.150*z*z*z*z*z;
    double eta = (Potential/Beta)*(1.0+(Beta/4.0)*
        (1-exp(-4.0/(DebyeLength*Beta))));

    double w(1.0);
    if( Beta == 0.0 || eta == -1.0 ){
        w = 1.0;
    }else if( std::isnan(eta) ){
        std::cout << "\nWarning! w being set to 1.0 "
            << "(Assuming high B field limit) but Phi/Beta is nan while "
            << "Beta != 0.";
        w = 1.0;
    }else{
        w = eta/(1+eta);
    } 
    
    double A = 0.678*w+1.543*w*w-1.212*w*w*w;
    double ElecFlux(0.0);
    if( Potential >= 0.0 ){ //!< For negative dust, do PHL
        ElecFlux = Pdata->ElectronDensity*
            sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(A+(1.0-A)*i_star)*
            exp(-Potential);
    }else{ //!< For positive dust, do OML
        ElecFlux = Pdata->ElectronDensity*
            sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*(1.0-Potential);
    }
    if(ElecFlux >= Underflows::Flux && ElecFlux != INFINITY 
        && ElecFlux == ElecFlux && ElecFlux < Overflows::Flux ){
        return ElecFlux;
    }

    static bool runOnce = true;
    std::string Warning = "\nError in PHLElectronFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double DTOKSIonFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug("\n\t\tIn DTOKSIonFlux:Term()\n\n");
    double IonFlux=0;

    //!< Positive grain, DeltaTot() > 1
    if( Sample->is_positive() ) 
        IonFlux = Flux::DTOKSElectronFlux(Pdata, Potential); 
    else    
        IonFlux = Flux::DTOKSElectronFlux(Pdata, Potential)*
            (1-Sample->get_deltatot());

    if(IonFlux >= Underflows::Flux && IonFlux != INFINITY && IonFlux == IonFlux 
        && IonFlux < Overflows::Flux ){
        return IonFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in DTOKSIonFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double DTOKSElectronFlux(const std::shared_ptr<PlasmaData> Pdata, 
        const double Potential){
    Mo_Debug("\n\t\tIn DTOKSElectronFlux:Term()\n\n");
    double ElecFlux(0.0);

    ElecFlux = Pdata->ElectronDensity*exp(-Potential)*
        sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));

    if(ElecFlux >= Underflows::Flux && ElecFlux != INFINITY 
        && ElecFlux == ElecFlux && ElecFlux < Overflows::Flux ){
        return ElecFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in DTOKSElectronFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double OMLElectronFlux(const std::shared_ptr<PlasmaData> Pdata, 
        const double Potential){
    Mo_Debug("\n\t\tIn OMLElectronFlux:Term()\n\n");

    double ElecFlux(0.0);
    if( Potential < 0.0 ){
        ElecFlux = Pdata->ElectronDensity*
            sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me))*(1-Potential);
    }else{  
        ElecFlux = Pdata->ElectronDensity*exp(-Potential)*
            sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
    }
    if(ElecFlux >= Underflows::Flux && ElecFlux != INFINITY 
        && ElecFlux == ElecFlux && ElecFlux < Overflows::Flux ){
        return ElecFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in OMLElectronFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double NeutralFlux(const std::shared_ptr<PlasmaData> Pdata){
    Mo_Debug("\n\t\tIn NeutralFlux:Term()\n\n");

    double NeutFlux = Pdata->NeutralDensity*
        sqrt(Kb*Pdata->NeutralTemp/(2*PI*Pdata->mi));

    if(NeutFlux >= Underflows::Flux && NeutFlux != INFINITY 
        && NeutFlux == NeutFlux && NeutFlux < Overflows::Flux ){
        return NeutFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in NeutralFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double EvaporationFlux(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    Mo_Debug("\n\t\tIn EvaporationFlux:Term()\n\n");
    double AmbientPressure = Pdata->NeutralDensity*Kb*Pdata->NeutralTemp;

    double StickCoeff = 1.0;
    double EvapFlux = (StickCoeff*Sample->get_surfacearea()*AvNo*
        (Sample->get_vapourpressure()-AmbientPressure))/
        sqrt(2*PI*Sample->get_atomicmass()*R*DustTemperature);

    //!< Missing "EvapFlux >= 0.0 &&" permits dust gaining mass from 
    //!< condensation of nearby cold nuclei.
    if( EvapFlux >= 0.0 && EvapFlux != INFINITY && EvapFlux == EvapFlux ){
        return EvapFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in EvaporationFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return 0.0;
}

double DeltaTherm(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata){
    C_Debug("\n\t\tIn DeltaTherm:Term()\n\n");

    double dtherm = (Richardson*Sample->get_temperature()*
        Sample->get_temperature()*exp(-(Sample->get_workfunction()*echarge)/
        (Kb*Sample->get_temperature())))/
        (echarge*OMLElectronFlux(Pdata, Sample->get_potential()));

    if(dtherm >= 0.0 && dtherm == dtherm && dtherm != INFINITY ){
        return dtherm;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in DeltaTherm()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return dtherm;
}

double ThermFlux(const Matter* Sample){
    C_Debug("\n\t\tIn DeltaTherm:Term()\n\n");

    double ThermFlux = (Richardson*Sample->get_temperature()*
        Sample->get_temperature()*exp(-(Sample->get_workfunction()*echarge)/
        (Kb*Sample->get_temperature())))/echarge;
    if(ThermFlux >= 0.0 && ThermFlux == ThermFlux && ThermFlux != INFINITY ){
        return ThermFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in ThermFlux()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return ThermFlux;
}

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
    if(ThermFlux >= 0.0 && ThermFlux == ThermFlux && ThermFlux != INFINITY ){
        return ThermFlux;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in ThermFluxSchottky()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return ThermFlux;
}

double DeltaSec(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata){
    C_Debug("\n\t\tIn DeltaSec:Term()\n\n");
    double ConvertKtoev(8.6173303e-5);
    double DeltaSec = sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem());
    if(DeltaSec >= 0.0 && DeltaSec == DeltaSec && DeltaSec != INFINITY ){
        return DeltaSec;
    }
    static bool runOnce = true;
    std::string Warning = "\nError in DeltaSec()!";
    Warning += " Return value badly specified\nReturning zero!\n";
    WarnOnce(runOnce,Warning);
    return DeltaSec;
}

}