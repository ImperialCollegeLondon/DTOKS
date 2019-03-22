/** @file ChargingTerms.cpp
 *  @brief Contains function definitions for charging terms
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#include "CurrentTerms.h"

namespace Term{
    double OMLe::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::OMLElectronFlux(Pdata,Potential);
    }
    double PHLe::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::PHLElectronFlux(Sample,Pdata,Potential);
    }
    double OMLi::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::OMLIonFlux(Sample,Pdata,Potential);
    }
    double SOMLi::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::SOMLIonFlux(Sample,Pdata,Potential);
    }
    double SMOMLi::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::SMOMLIonFlux(Sample,Pdata,Potential);
    }
    double TEE::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::ThermFlux(Sample);
    }
    double TEESchottky::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::ThermFluxSchottky(Sample,Pdata,Potential);
    }
    double THSe::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){

        double TiTe = Pdata->IonTemp/Pdata->ElectronTemp;
        double MassRatio = Pdata->mi/Me;
        double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/
            (Pdata->ElectronDensity*pow(echarge,2)));
        double DebyeLength_Tilde = DebyeLength/Sample->get_radius();
        double Betai=Sample->get_radius()/(sqrt(2.0*Kb*Pdata->IonTemp/
            (PI*Pdata->mi))/(echarge*Pdata->MagneticField.mag3()/Pdata->mi));
        double Betae=Betai*sqrt(TiTe*MassRatio);

        double a = 1.522;
        double b = -0.9321;
        double c = 0.3333;
        double d = 1.641;

        double TDR=erf(a*pow(DebyeLength,b)*pow(TiTe,c));

        double e = 2.82;
        double f = 0.4772;
        double g = 0.1315;
        double h = 0.7364;
        
        double Repelled_Species_Current=(exp(-g*Betae)+e*exp(-DebyeLength_Tilde)*
            pow((Betai/(Betai+1)),h));
        return Repelled_Species_Current;
        
        double Coeff = Pdata->Z*sqrt(TiTe/MassRatio);
        double Attacted_Species_Current_T1=Coeff*(exp(-f*Betai)*(TDR+(1-TDR)*d*
            (1.0/sqrt(TiTe)+1.0/sqrt(DebyeLength_Tilde)))+
            e*pow((Betai/(Betai+1)),h));
        double Attacted_Species_Current_T2=(Coeff/TiTe)*exp(-f*Betai)*TDR;

        double a1= Attacted_Species_Current_T1/Repelled_Species_Current;
        double b1= Attacted_Species_Current_T2/Repelled_Species_Current;

        return fabs(LambertW(exp(a1/b1)/b1)-a1/b1);
    }
    double THSi::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){

        double TiTe = Pdata->IonTemp/Pdata->ElectronTemp;
        double MassRatio = Pdata->mi/Me;
        double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/
            (Pdata->ElectronDensity*pow(echarge,2)));
        double DebyeLength_Tilde = DebyeLength/Sample->get_radius();
        double Betai=Sample->get_radius()/(sqrt(2.0*Kb*Pdata->IonTemp/
            (PI*Pdata->mi))/(echarge*Pdata->MagneticField.mag3()/Pdata->mi));
        double Betae=Betai*sqrt(TiTe*MassRatio);

        double a = 1.522;
        double b = -0.9321;
        double c = 0.3333;
        double d = 1.641;

        double TDR=erf(a*pow(DebyeLength,b)*pow(TiTe,c));

        double e = 2.82;
        double f = 0.4772;
        double g = 0.1315;
        double h = 0.7364;

        double Coeff = Pdata->Z*sqrt(TiTe/MassRatio);
        double Attacted_Species_Current_T1=Coeff*(exp(-f*Betai)*(TDR+(1-TDR)*d*
            (1.0/sqrt(TiTe)+1.0/sqrt(DebyeLength_Tilde)))+
            e*pow((Betai/(Betai+1)),h));
        double Attacted_Species_Current_T2=(Coeff/TiTe)*exp(-f*Betai)*TDR;

        return Attacted_Species_Current_T1+Attacted_Species_Current_T2*Potential;
    }
    double DTOKSi::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::DTOKSIonFlux(Sample,Pdata,Potential);
    }
    double DTOKSe::Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential){
        return Flux::DTOKSElectronFlux(Pdata,Potential);
    }
}