/** @file PlasmaFluxes.cpp
 *  @brief Contains functions defining the default plasma fluxes
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#include "PlasmaFluxes.h"

namespace Flux{

double OMLIonFlux(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    H_Debug("\n\tIn IonFlux():");

    double IonFlux=0;

    if( Potential >= 0 ) 
        IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
                *(1+Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
    else    
        IonFlux = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi))
                *exp(Potential*(Pdata->ElectronTemp/Pdata->IonTemp));

    assert(IonFlux >= 0);
    return IonFlux;
}

double SOMLIonFlux(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug( "\tSOMLIonFlux:Term\n\n");
    double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
    double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()
        /IonThermalVelocity;
    double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
    if( uz == 0.0 ){
        return Flux::OMLIonFlux(Sample,Pdata,Potential);
    }
    if( Potential >= 0.0 ){
        double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+exp(-uz*uz)/2.0;
        double s2 = sqrt(PI)*erf(uz)/(2.0*uz);
        return Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*
                (s1+s2*Pdata->Z*Potential/Tau);
    }else{
        double uzp = uz+sqrt(-Pdata->Z*Potential/Tau);
        double uzm = uz-sqrt(-Pdata->Z*Potential/Tau);
        return Pdata->IonDensity*IonThermalVelocity*(1.0/(8.0*sqrt(2.0)*uz))*
                ((1.0+2.0*(uz*uz+Pdata->Z*Potential/Tau))*(erf(uzp)+erf(uzm))
                +(2.0/sqrt(PI))*(uzp*exp(-uzm*uzm)+uzm*exp(-uzp*uzp)));
    }
}

double SMOMLIonFlux(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug( "\tSMOMLIonFlux:Term\n\n");
    double HeatCapacityRatio = 5.0/3.0;
    double MassRatio = Pdata->mi/Me;
    double Tau = Pdata->IonTemp/Pdata->ElectronTemp;
    double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
    double uz = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()
        /IonThermalVelocity;
    if( uz == 0.0 ){ //!< Handle zero-relative velocity case
        return Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0))*(1.0+(1.0/Tau)*
            (Potential*Pdata->Z-0.5*log(2.0*PI*(1.0+HeatCapacityRatio*Tau)
            /MassRatio)));
    }
    if( Potential >= 0.0 ){ //!< For negative dust, do SMOML
        double s1 = sqrt(PI)*(1.0+2.0*uz*uz)*erf(uz)/(4.0*uz)+exp(-uz*uz)/2.0;
        double s2 = sqrt(PI)*erf(uz)/(2.0*uz);
        return Pdata->IonDensity*(IonThermalVelocity/sqrt(2.0*PI))*(s1+(s2/Tau)*
            (Potential*Pdata->Z-0.5*log(2.0*PI*(1.0+HeatCapacityRatio*Tau)
            /MassRatio)));
    }else{  //!< For Positive dust, do SOML
        double uzp = uz+sqrt(-Pdata->Z*Potential/Tau);
        double uzm = uz-sqrt(-Pdata->Z*Potential/Tau);
        return Pdata->IonDensity*IonThermalVelocity*(1.0/(8.0*sqrt(2.0)*uz))*
                ((1.0+2.0*(uz*uz+Pdata->Z*Potential/Tau))*(erf(uzp)+erf(uzm))
                +(2.0/sqrt(PI))*(uzp*exp(-uzm*uzm)+uzm*exp(-uzp*uzp)));
    }
}

double PHLElectronFlux(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    Mo_Debug( "\tPHLElectronFlux:Term\n\n");
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

    if( Potential >= 0.0 ){ //!< For negative dust, do PHL
        return Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*
            (A+(1.0-A)*i_star)*exp(-Potential);
    }else{ //!< For positive dust, do OML
        return Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2.0*PI*Me))*
            (1.0-Potential);
    }
}

double DTOKSIonFlux(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
H_Debug("\tIn DTOKSIonFlux:Term\n\n");
    double IonFlux=0;

    //!< Positive grain, DeltaTot() > 1
    if( Sample->is_positive() ) IonFlux = Flux::DTOKSElectronFlux(Pdata, Potential); 
    else    IonFlux = Flux::DTOKSElectronFlux(Pdata, Potential)*(1-Sample->get_deltatot());
    assert(IonFlux >= 0);
    return IonFlux;
}

double DTOKSElectronFlux(const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    H_Debug("\tIn DTOKSElectronFlux:Term\n\n");
    return Pdata->ElectronDensity*exp(-Potential)*
        sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
}

double OMLElectronFlux(const std::shared_ptr<PlasmaData> Pdata, const double Potential){
    H_Debug("\tIn OMLElectronFlux:Term(const std::shared_ptr<PlasmaData> Pdata, const double Potential)\n\n");
    if( Potential < 0.0 ){
        return Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me))*
            (1-Potential);
    }else{  
        return Pdata->ElectronDensity*exp(-Potential)*
            sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
    }
}

double NeutralFlux(const std::shared_ptr<PlasmaData> Pdata){
    H_Debug("\tIn NeutralFlux(const std::shared_ptr<PlasmaData> Pdata)\n\n");

    return Pdata->NeutralDensity*sqrt(Kb*Pdata->NeutralTemp/(2*PI*Pdata->mi));
}

double EvaporationFlux(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn EvaporationFlux():\n\n");
    double AmbientPressure = Pdata->NeutralDensity*Kb*Pdata->NeutralTemp;

    double StickCoeff = 1.0;
    return (StickCoeff*Sample->get_surfacearea()*AvNo*
        (Sample->get_vapourpressure()-AmbientPressure))/
        sqrt(2*PI*Sample->get_atomicmass()*R*DustTemperature);
}

double DeltaTherm(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata){
    C_Debug("\tIn DeltaTherm(const Matter* Sample)\n\n");

    double dtherm = (Richardson*Sample->get_temperature()*
        Sample->get_temperature()*
        exp(-(Sample->get_workfunction()*echarge)/
            (Kb*Sample->get_temperature())))/
            (echarge*OMLElectronFlux(Pdata, Sample->get_potential()));
    assert(dtherm >= 0.0 && dtherm == dtherm && dtherm != INFINITY );
    return dtherm;
}

double ThermFlux(const Matter* Sample){
    C_Debug("\tIn DeltaTherm()\n\n");

    double ThermFlux = (Richardson*Sample->get_temperature()*
        Sample->get_temperature()*
        exp(-(Sample->get_workfunction()*echarge)/
            (Kb*Sample->get_temperature())))/
            echarge;
    assert(ThermFlux >= 0.0 && ThermFlux == ThermFlux && ThermFlux != INFINITY );
    return ThermFlux;
}

double ThermFluxSchottky(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential ){
    C_Debug("\tIn DeltaTherm( const double Potential )\n\n");
    //!< Returns the flux of electrons due to Thermionic emission
    //!< Following the Richard-Dushmann formula with Schottky Correction.
    //!< Negative Dust grains have the normal Richard-Dushmann form, dependant 
    //!< on work funciton
    //!< Positive Dust grains are the same flux multiplied by (1-phi)
    //!< e^(-e phi/(kb Td))
    if( Potential >= 0.0 ){
        return Richardson*Sample->get_temperature()*
            Sample->get_temperature()*
            exp(-echarge*Sample->get_workfunction()/
            (Kb*Sample->get_temperature()))/echarge;
    }else{

        return Richardson*Sample->get_temperature()*
            Sample->get_temperature()*(1.0-Potential
            *(Pdata->ElectronTemp/Sample->get_temperature()))
            *exp((-echarge*Sample->get_workfunction()-
            Potential*Kb*Pdata->ElectronTemp)/(Kb*Sample->get_temperature()))/
            echarge;
    }
}

double DeltaSec(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata){
    C_Debug("\tIn DeltaSec()\n\n");
    double ConvertKtoev(8.6173303e-5);
    return sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem());
}

}