/** @file ForceTerms.cpp
 *  @brief Contains function definitions for force terms
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#include "ForceTerms.h"

namespace Term{

threevector Gravity::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct Centrifugal::Evaluate()\n\n");
    threevector return_Vec(0.0,0.0,-9.81);
    return return_Vec;
};

threevector Centrifugal::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct Centrifugal::Evaluate()\n\n");
    threevector returnval(
        Sample->get_velocity().gety()*Sample->get_velocity().gety()/
        Sample->get_position().getx(),
        -Sample->get_velocity().getx()*Sample->get_velocity().gety()/
        Sample->get_position().getx(),
        0.0);
    return returnval;
}

threevector LorentzForce::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct LorentzForce::Evaluate()\n\n");
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

threevector SOMLIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct SOMLIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata)\n\n");
    return (Pdata->PlasmaVel-Sample->get_velocity())*
            Pdata->mi*(1.0/sqrt(Kb*Pdata->IonTemp/Pdata->mi))*
            Flux::SOMLIonFlux(Sample,Pdata,Sample->get_potential())*(1.0/Sample->get_mass());
}

threevector SMOMLIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn SMOMLIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata)\n\n");
    return (Pdata->PlasmaVel-Sample->get_velocity())*
            Pdata->mi*(1.0/sqrt(Kb*Pdata->IonTemp/Pdata->mi))*
            Flux::SMOMLIonFlux(Sample,Pdata,Sample->get_potential())*(1.0/Sample->get_mass());
}

threevector DTOKSIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct DTOKSIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata)\n\n");
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

threevector DUSTTIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct DUSTTIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata)\n\n");
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

threevector HybridIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct HybridIonDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata)\n\n");

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

threevector NeutralDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct NeutralDrag::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata)\n\n");
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

threevector RocketForce::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata){
    F_Debug("\tIn struct RocketForce::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata)\n\n");
        threevector returnvec(0.0,0.0,0.0);

    if( Sample->is_liquid() ){
        double Pv_plus = 
            Sample->probe_vapourpressure(Sample->get_temperature());
        double Pv_minus = Sample->probe_vapourpressure(OldTemp);
        OldTemp = Sample->get_temperature();
        returnvec = Sample->get_surfacearea()*((Pv_plus-Pv_minus)/
            Sample->get_radius())*Pdata->MagneticField.getunit();
    }
    return returnvec*(1.0/Sample->get_mass());
}

}