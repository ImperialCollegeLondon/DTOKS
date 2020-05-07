/** @file HeatTerms.cpp
 *  @brief Contains function definitions for heating terms
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#include "HeatTerms.h"

// ***************************** HEATING MODELS ***************************** //
namespace Term{
    
//!< Using Stefan-Boltzmann Law, returns Energy lost per second in Kila Joules
double EmissivityModel::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::EmissivityModel(const double DustTemperature):"
        << "\n\n");
    //!< Energy emitted from a sample converted to kJ
    return -Sample->get_emissivity()*Sample->get_surfacearea()*Sigma
            *(pow(DustTemperature,4)-pow(Pdata->AmbientTemp,4));
}

//!< http://users.wfu.edu/ucerkb/Nan242/L06-Vacuum_Evaporation.pdf, 
//!< https://en.wikipedia.org/wiki/Hertz%E2%80%93Knudsen_equation
//!< Using Hertz–Knudsen equation, returns Energy lost per second in Kila Joules
//!< MASS LOSS EQUATION 
//!< http://www.leb.eei.uni-erlangen.de/
//!< winterakademie/2006/result/content/course01/pdf/0102.pdf
double EvaporationModel::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::EvaporationModel():\n\n");

    //!< Approximate emitted energy as mean maxwell
    //!< This used to be multiplied by 1000 which is why it was significant.
    //!< Converted to kJ.
    double MaxwellEnergy = (3.0*Kb*DustTemperature/(2.0*1000.0)); 
    //double MaxwellMeanVelocity = sqrt((8*Kb*DustTemperature)/
    //(AtomicMass*1000*AMU));

    double EvapFlux = Flux::EvaporationFlux(Sample, Pdata, DustTemperature);


    //!< See GrainStructs.h for more info on 'bondenergy'. 
    //!< Added to account for energy lost by breaking bonds.
    return -EvapFlux*(MaxwellEnergy+Sample->get_bondenergy()/AvNo); 
}

//!< VERY APPROXIMATE MODEL: Atmosphere assumed to be 300 degrees always,
//!< rough heat transfer coefficient is use
//!< https://en.wikipedia.org/wiki/Newton%27s_law_of_cooling
 double NewtonCooling::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::NewtonCooling():\n\n");
    static bool runOnce = true;
    std::string Warning = "In HeatingModel::NewtonCooling():\nHeatTransair ";
    Warning += "Coefficient wrong for Tungsten, Beryllium, Graphite, Helium,";
    Warning += " Lithium & Molybdenum.";
    WarnOnce(runOnce,Warning);

    return (Sample->get_heattransair()*Sample->get_surfacearea()*
        (DustTemperature-Pdata->AmbientTemp)); 
}   

double NeutralHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::NeutralHeatFlux(const double DustTemperature)"
        << ":\n\n");
    return Sample->get_surfacearea()*Flux::NeutralFlux(Pdata)*
        (Pdata->NeutralTemp-DustTemperature)*Kb;
}

// ************************** SOML/OML/PHL MODELS ************************** //

double SEE::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::SEE():\n\n");
    double ConvertKtoev(8.6173303e-5);
    return -Sample->get_surfacearea()*Flux::PHLElectronFlux(Sample,Pdata,Sample->get_potential())*
        sec(Pdata->ElectronTemp*ConvertKtoev,Sample->get_elem())*echarge*
        (3.0+Sample->get_workfunction()); 
}

double TEE::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::TEE():\n\n");
    if( Sample->get_potential() >= 0.0){
        return -Sample->get_surfacearea()*Richardson*pow(DustTemperature,2)*
            exp(-Sample->get_workfunction()*echarge/(Kb*DustTemperature))*
            (2.0*Kb*DustTemperature+echarge*Sample->get_workfunction())/echarge;
    }else{
        //<! See [1] V. E. Fortov, A. V. Ivlev, S. A. Khrapak, A. G. Khrapak, 
        //<! and G. E. Morfill, Phys. Rep. 421, 1 (2005).
        //<! Page 23, section 3.3
        return -Sample->get_surfacearea()*Richardson*pow(DustTemperature,2)*
            (1.0-Sample->get_potential()*
            (Pdata->ElectronTemp/Sample->get_temperature()))*
            exp((Sample->get_potential()*Kb*Pdata->ElectronTemp-
            Sample->get_workfunction()*echarge)/(Kb*DustTemperature))*
            (2.0*Kb*DustTemperature+echarge*Sample->get_workfunction())/echarge;
    }
}

double PHLElectronHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){ 

    H_Debug("\n\tIn HeatingModel::PHLElectronHeatFlux():\n\n");
    //!< Only for a negative grain
    return Sample->get_surfacearea()*Flux::PHLElectronFlux(Sample,Pdata,Sample->get_potential())*
        Kb*(Pdata->ElectronTemp-DustTemperature);
}

double SOMLIonHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){ 
    H_Debug("\n\tIn HeatingModel::SOMLIonHeatFlux(const double DustTemperature):"
        <<"\n\n");
    // Assuming Re = 0
    return Sample->get_surfacearea()*(1.0-Sample->get_re())*
        Flux::SOMLIonFlux(Sample,Pdata,Sample->get_potential())*Pdata->IonTemp*Kb; 
}

double SOMLNeutralRecombination::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature)
{
    H_Debug("\n\tIn HeatingModel::SOMLNeutralRecombination():\n\n");
    //!< Neutral Recombination assuming Rn=0; fraction of backscattered 
    //!< ions/neutrals is zero.
    return Sample->get_surfacearea()*(1.0-Sample->get_rn())*
        (14.7*echarge - 2.0*Kb*DustTemperature)*
        Flux::SOMLIonFlux(Sample,Pdata,Sample->get_potential()); 
}

double SMOMLIonHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){ 
    H_Debug("\n\tIn HeatingModel::SMOMLIonHeatFlux(const double DustTemperature):"
        << "\n\n");
    // Assuming Re = 0
    return Sample->get_surfacearea()*(1.0-Sample->get_re())*
        Flux::SMOMLIonFlux(Sample,Pdata,Sample->get_potential())*Pdata->IonTemp*Kb;
}

double SMOMLNeutralRecombination::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature)
{
    H_Debug("\n\tIn HeatingModel::SMOMLNeutralRecombination():\n\n");
    // Neutral Recombination assuming Rn=0; fraction of backscattered 
    // ions/neutrals is zero.
    return Sample->get_surfacearea()*(1.0-Sample->get_rn())*
        (14.7*echarge - 2.0*Kb*DustTemperature)*
        Flux::SMOMLIonFlux(Sample,Pdata,Sample->get_potential()); 
}


double OMLElectronHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){ 
    H_Debug("\n\tIn HeatingModel::OMLElectronHeatFlux():\n\n");
    // Only for a negative grain
    return Sample->get_surfacearea()*Flux::OMLElectronFlux(Pdata,Sample->get_potential())*
        Kb*(Pdata->ElectronTemp-DustTemperature);
}

// ****************************** DTOKS MODELS ****************************** //

 double DTOKSSEE::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::DTOKSSEE():\n\n");
    //!< Electrons released all re-captured by positive dust grain
    double SEE=0; 
    if( !Sample->is_positive() ) // Dust grain is negative
        SEE = Sample->get_surfacearea()*
            Flux::OMLElectronFlux(Pdata,Sample->get_potential())*Sample->get_deltasec()*
            echarge*(3.0+Sample->get_workfunction()); 
    //!< else, electrons captured by positive grain, return zero
    return -SEE; 
}

double DTOKSTEE::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::DTOKSTEE():\n\n");
    //!< Electrons released all re-captured by positive dust grain
    double TEE=0;
    if( !Sample->is_positive() ) // Dust grain is negative
        TEE = Sample->get_surfacearea()*Sample->get_deltatherm()*
            Flux::OMLElectronFlux(Pdata,Sample->get_potential())*
            (2.0*Kb*DustTemperature+echarge*Sample->get_workfunction()); 
    //!< else, electrons captured by positive grain, return zero

    return -TEE;
}

double DTOKSIonHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){ 
    H_Debug("\n\tIn HeatingModel::DTOKSIonHeatFlux(const double DustTemperature):"
        << "\n\n");
    // Assuming Re = 0
    return (Sample->get_surfacearea()*(1.0-Sample->get_re())*
        Flux::DTOKSIonFlux(Sample,Pdata,Sample->get_potential())*Pdata->IonTemp*Kb)*
        (2.0+2.0*Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp)+
        pow(Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp),2.0))/
        (1.0+Sample->get_potential()*(Pdata->ElectronTemp/Pdata->IonTemp)); 
    // Convert from Joules to KJ
}


double DTOKSNeutralRecombination::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::NeutralRecombination():\n\n");
    // Neutral Recombination assuming Rn=0; fraction of backscattered 
    // ions/neutrals is zero.
    return Sample->get_surfacearea()*(1.0-Sample->get_rn())*
        (14.7*echarge - Kb*DustTemperature)*
        Flux::DTOKSIonFlux(Sample,Pdata,Sample->get_potential()); 
}

double DTOKSElectronHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){ 
    H_Debug("\n\tIn HeatingModel::ElectronHeatFlux():\n\n");
    // Only for a negative grain
    return 2.0*sqrt(2.0*PI)*Sample->get_radius()*Sample->get_radius()
        *Flux::DTOKSElectronFlux(Pdata,Sample->get_potential())*Kb*
        (Pdata->ElectronTemp-DustTemperature);
}

// ****************************** DUSTT MODELS ****************************** //


//!< This is the heat flux of ions incident on the dust grain surface ONLY
//!< Calculated from equation (31) & (32), page 29, from the following reference
//!< S.I. Krasheninnikov, R.D. Smirnov, and D.L. Rudakov,
//!< Plasma Phys. Control. Fusion 53, 083001 (2011).
double DUSTTIonHeatFlux::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature){
    H_Debug("\n\tIn HeatingModel::DUSTTIonHeatFlux(const double DustTemperature):"
        << "\n\n");
    // Assuming Re = 0
//  H1_Debug( "\nSample->get_re() = " << Sample->get_re() );
    if( Sample->get_re() > 0.1 ){ // Uncomment when Sample->get_re() is calculated
        static bool runOnce = true;
        std::string Warning = "In HeatingModel::DUSTTIonHeatFlux(double ";
        Warning += "DustTemperature)\nSample->get_re() > 0.1. Ion Heat Flux affected by ";
        Warning += "backscattering by more than 10%!";
        WarnOnce(runOnce,Warning);
    }
    double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
    double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
    double RelativeVelocity = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/
        IonThermalVelocity;
    double RenormalisedPotential = (Pdata->Z*Sample->get_potential())/
        TemperatureRatio;
    double Term1(0.0);
    double Term2(0.0), Term2Coeff(0.0);
    if( RelativeVelocity != 0.0 ){
        Term2Coeff = (1.0/(4.0*RelativeVelocity))*(3.0+12.0*RelativeVelocity*
            RelativeVelocity+4.0*pow(RelativeVelocity,4.0)+
            2.0*RenormalisedPotential*
            (1+2.0*RelativeVelocity*RelativeVelocity));
    }else if(RelativeVelocity == 0.0){
        return  Sample->get_surfacearea()*(1.0-Sample->get_re())*Pdata->IonDensity*
            Pdata->IonTemp*Kb*IonThermalVelocity;
    }
    if( RenormalisedPotential >= 0.0 ){ // Negative dust grain
        Term1 =(1.0/(2.0*sqrt(PI)))*(5.0+2.0*RelativeVelocity*RelativeVelocity
            -2.0*RenormalisedPotential)*exp(-RelativeVelocity*RelativeVelocity); 
        Term2 = erf(RelativeVelocity)*Term2Coeff; 
    }else{ // Positive dust grain
        assert( Sample->get_potential() <= 0.0 );
        double Coeff1 = 5.0; 
        double Coeff2 = 5.0; 
        if( RelativeVelocity != 0.0 ){
            Coeff1 = 5.0+2.0*RelativeVelocity*RelativeVelocity
                -((3.0+2.0*RelativeVelocity*RelativeVelocity)/RelativeVelocity)*
                sqrt(-RenormalisedPotential);
            Coeff2 = 5.0+2.0*RelativeVelocity*RelativeVelocity
                +((3.0+2.0*RelativeVelocity*RelativeVelocity)/RelativeVelocity)*
                sqrt(-RenormalisedPotential);
        }

        double uip = RelativeVelocity+sqrt(-RenormalisedPotential);
        double uim = RelativeVelocity-sqrt(-RenormalisedPotential);
        Term1 = (1.0/(4.0*sqrt(PI)))*(Coeff1*exp(-uip*uip) +
            Coeff2*exp(-uim*uim)); 
        Term2 = (Term2Coeff/2.0)*(erf(uip)+erf(uim));
    }
    
    assert((Term1 + Term2) > 0 && (Term1 + Term2) != INFINITY 
        && (Term1 + Term2) == (Term1 + Term2));

    return Sample->get_surfacearea()*(1.0-Sample->get_re())*Pdata->IonDensity
        *Pdata->IonTemp*Kb*IonThermalVelocity*(Term1 + Term2);
    // Convert from Joules to KJ
}

}
