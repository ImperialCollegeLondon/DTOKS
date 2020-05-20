/** @file ChargingTerms.cpp
 *  @brief Contains function definitions for charging terms
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#include "CurrentTerms.h"

namespace Term{
    double OMLe::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return -Flux::OMLElectronFlux(Pdata,Potential);
    }
    double PHLe::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return -Flux::PHLElectronFlux(Sample,Pdata,Potential);
    }
    double OMLi::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return Pdata->Z*Flux::OMLIonFlux(Sample,Pdata,Potential);
    }
    double MOMLi::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return Pdata->Z*Flux::MOMLIonFlux(Sample,Pdata,Potential);
    }
    double SOMLi::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return Pdata->Z*Flux::SOMLIonFlux(Sample,Pdata,Potential);
    }
    double SMOMLi::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return Pdata->Z*Flux::SMOMLIonFlux(Sample,Pdata,Potential);
    }
    double TEEcharge::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return Flux::ThermFlux(Sample);
    }
    double TEESchottky::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return Flux::ThermFluxSchottky(Sample,Pdata,Potential);
    }
    double SEEcharge::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return -Flux::DeltaSec(Sample,Pdata);
    }
    double THSe::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        double TiTe = Pdata->IonTemp/Pdata->ElectronTemp;
        double MassRatio = Pdata->mi/Me;
        double Betae=Sample->get_radius()/(sqrt(Kb*Pdata->ElectronTemp/
            Me)/(echarge*Pdata->MagneticField.mag3()/Me));

        double a = 0.2293;
        
        double Repelled_Species_Flux=
		    Pdata->ElectronDensity*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me))*(1.0+exp(-a*Betae));
        return -Repelled_Species_Flux;
    }
    double THSi::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        double TiTe = Pdata->IonTemp/Pdata->ElectronTemp;
        double MassRatio = Pdata->mi/Me;
        double Betai=Sample->get_radius()/(sqrt(Kb*Pdata->IonTemp/
            (Pdata->mi))/(Pdata->Z*echarge*Pdata->MagneticField.mag3()/Pdata->mi));

        double b = 1.559;
        double c = 0.5636;

        double Coeff = Pdata->IonDensity*sqrt(Kb*Pdata->IonTemp/(2*PI*Pdata->mi));
        double Attacted_Species_Flux_T1=Coeff*(1+exp(-b*pow(Betai,c)));
        double Attacted_Species_Flux_T2=(Coeff*Pdata->Z/TiTe)*2.0*exp(-b*pow(Betai,c));

        return Attacted_Species_Flux_T1+Attacted_Species_Flux_T2*Potential;

    }
    double DTOKSi::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return Pdata->Z*Flux::DTOKSIonFlux(Sample,Pdata,Potential);
    }
    double DTOKSe::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        return -Flux::DTOKSElectronFlux(Pdata,Potential);
    }
    double CW::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        //!< Following Semi-empirical fit to Sceptic results as detailed in Chris Willis
        //!< Thesis,
        //!< https://spiral.imperial.ac.uk/handle/10044/1/9329, pages 68-70
        double A = Pdata->mi;
        double b = Pdata->IonTemp/Pdata->ElectronTemp;
        double DebyeLength=sqrt((epsilon0*Kb*Pdata->ElectronTemp)/
            (Pdata->ElectronDensity*pow(echarge,2)));
        double Rho = Sample->get_radius()/DebyeLength;
        double Rho_OML = 1.25*pow(b,0.4);
        double Rho_Upper = 50.0;
        double Current(0.0);
        if( Rho <= Rho_OML ){ //!< This is the OML Limit
            if( b <=2 ){ //!< Ti <= 2.0*Te
                Current = 0.405*log(Pdata->A)+(0.253+0.021*log(Pdata->A))*log(b)+
                    2.454 - Potential;
            }else{ //!< Ti > 2.0*Te
                Current = 0.401*log(Pdata->A)+(-0.122+0.029*log(Pdata->A))*log(b)+
                    2.698 - Potential;
            }
            //!< This is the transition region
        }else if( Rho <= Rho_Upper && Rho > Rho_OML ){ 
            double Gradient = (log(Rho/Rho_Upper)/log(Rho_Upper/Rho_OML))+1.0;
            double DeltaPhi = 0.5*log(2.0*PI*(Me/Pdata->mi)*(1.0+5.0*b/3.0))*
                Gradient;
            return exp(-Potential)-sqrt(b*(Me/Pdata->mi))*(1+Potential/b-DeltaPhi/b);
        }else if( Rho > Rho_Upper ){ //!< This is the MOML limit
            if( b <=2 ){ //!< Ti <= 2.0*Te
                Current = 0.456*log(Pdata->A)+3.179-Potential;
            }else{ //!< Ti > 2.0*Te
                Current = 0.557*log(Pdata->A)-(0.386+0.024*log(Pdata->A))*log(b)+
                    3.399 - Potential;
            }
        }else{
            std::cerr << "\nError in ChargingModel::solveCW()!"
                << " Rho is poorly defined!\n\n";
        }
        return Current;
    }

    double MOMLWEM::Evaluate(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata, const double Potential){
        // Solve the Modified orbital motion limited potential for large emitting dust grains.
        // See the paper by Minas and Nikoleta, equation (1) and (2)
        // N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
        double DeltaTot = Sample->get_deltatot();
        if( DeltaTot >= 1.0 ){//!< In this case we have a positive grain!
            return Pdata->Z*Flux::OMLIonFlux(Sample,Pdata,Potential)+Flux::OMLElectronFlux(Pdata,Potential);
        }else{
            double HeatCapacityRatio = 1.0;
            double TemperatureRatio = Pdata->IonTemp/Pdata->ElectronTemp;
            double MassRatio = Pdata->mi/Me;
            double Ionization = Pdata->Z;       // Ionization
            double IonThermalVelocity = sqrt((Kb*Pdata->IonTemp)/Pdata->mi);
            double PlasmaFlowSpeed = (Pdata->PlasmaVel-Sample->get_velocity()).mag3()/IonThermalVelocity;
            
            double Delta_Phi_em = 0.5*log((2.0*PI/MassRatio)*
                (1+HeatCapacityRatio*TemperatureRatio)/pow(1.0-DeltaTot,2.0));
            // Uncomment following line to compare MOMLWEM results with figure (1) of paper:
            // N. Rizopoulou and M. Bacharis, Phys. Plasmas 25, (2018).
            //std::cout << DeltaTot << "\t" << Delta_Phi_em << "\n";
            double Arg = sqrt(2*PI*TemperatureRatio*(1+HeatCapacityRatio*TemperatureRatio))
                *exp(TemperatureRatio);
            double Potential = -1.0*(TemperatureRatio/Ionization+Delta_Phi_em/Ionization-LambertW(Arg));
            if( Potential < 0.0 ){
                return Pdata->Z*Flux::OMLIonFlux(Sample,Pdata,Potential)+Flux::OMLElectronFlux(Pdata,Potential);
            }
            double IonCurrent = sqrt(TemperatureRatio/MassRatio)*
                (1.0-(-Potential-Delta_Phi_em)/TemperatureRatio);
            double ElectronCurrent = (1.0-DeltaTot)*
                exp(-Potential);
            double TotalCurr = Pdata->Z*IonCurrent - ElectronCurrent;
            //!< Sanity check sensible return value
            if(TotalCurr >= Underflows::Flux && TotalCurr != INFINITY 
                && TotalCurr == TotalCurr && TotalCurr < Overflows::Flux ){
                return TotalCurr;
            }
            //!< If return value is not well defined, print error and return 0.
            static bool runOnce = true;
            std::string Warning = "\nError in MOMLWEM:Evaluate()!";
            Warning += " Return value badly specified\nReturning zero!\n";
            WarnOnce(runOnce,Warning);
            return 0.0;
        }
    }
}
