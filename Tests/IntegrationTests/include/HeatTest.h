#include "HeatingModel.h"

int HeatTest(char Element, std::string HeatType, double accuracy){
    clock_t begin = clock();
    // ********************************************************** //
    // FIRST, define program default behaviour

    // Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
    char EmissivityModel = 'c';     // Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
    char ExpansionModel = 'c';  // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
                                                    // and (z)ero expansion
    char HeatCapacityModel = 'c';   // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
    char BoilingModel = 'y';    // Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
                                                    // and (t)homson
    char BreakupModel = 'n';    // Possible values 'r', 'e', 'b'  and 'n': Corresponding to (r)otational, (e)lectrostatic, (b)oth and (n)o

    std::string Name="constant";    // Describes heating model

    // Parameters describing the heating model
    double Power=1e-11;     // Kilo-Watts power in addition to heating model powers
    double Radius=5e-8;     // m
    double Temp=300;        // K
    double TimeStep=1e-10;  // s
    Matter *Sample;         // Define the sample matter type

    // Set to true all heating models that are wanted
    // Set to true all heating models that are wanted
    bool RadiativeCooling = false;
    bool EvaporativeCooling = false;
    bool NewtonCooling = false;
    bool NeutralHeatFlux = false;
    bool SOMLIonHeatFlux = false;
    bool SOMLNeutralRecombination = false;
    bool SMOMLIonHeatFlux = false;
    bool SMOMLNeutralRecombination = false;
    bool TEE = false;
    bool SEE = false;
    bool PHLElectronHeatFlux = false;
    bool OMLElectronHeatFlux = false;
    bool DTOKSTEE = false;
    bool DTOKSSEE = false;
    bool DTOKSIonHeatFlux = false;
    bool DTOKSNeutralRecomb = false;
    bool DTOKSElectronHeatFlux = false;
    bool DUSTTIonHeatFlux = false;

    PlasmaData *Pdata = new PlasmaData;
    Pdata->NeutralDensity = 3e19;    // m^-3, Neutral density
    Pdata->ElectronDensity = 8e17;   // m^-3, Electron density
    Pdata->IonDensity = 8e17;        // m^-3, Electron density
    double Potential = 1;            // arb, assumed negative, potential normalised to dust temperature, (-e*phi)/(Kb*Td)
    Pdata->IonTemp = 10*1.16e4;      // K, Ion Temperature
    Pdata->ElectronTemp = 10*1.16e4; // K, Electron Temperature, convert from eV
    Pdata->NeutralTemp = 10*1.16e4;  // K, Neutral Temperature, convert from eV
    Pdata->mi           = Mp;        // kg, Ion Mass
    Pdata->A            = 1.0;       // arb, Atomic Number
    Pdata->AmbientTemp = 0;

    if( HeatType == "ConstantHeating" ){
    }else if( HeatType == "ThermalRadiation" ){
        RadiativeCooling = true;
    }else if( HeatType == "ElectronHeatingRadiation" ){
        RadiativeCooling = true;
        DTOKSElectronHeatFlux = true;
        NeutralHeatFlux = true;
    }else if( HeatType == "PlasmaHeatingRadiation" ){
        RadiativeCooling = true;
        NeutralHeatFlux = true;
        DTOKSElectronHeatFlux = true;
        DTOKSIonHeatFlux = true;
    }else if( HeatType == "PlasmaHeatingNeutralRecomb" ){
        NeutralHeatFlux = true;
        DTOKSElectronHeatFlux = true;
        DTOKSIonHeatFlux = true;
        DTOKSNeutralRecomb = true;
    }else{
        std::cout << "\n\nInput not recognised! Exiting program\n.";
        return -1;
    }

    std::vector<HeatTerm*> HeatTerms;
    // Models and ConstModels are placed in an array in this order:
    std::array<bool, 18> HeatModels = 
        {RadiativeCooling, EvaporativeCooling, NewtonCooling, NeutralHeatFlux,
            SOMLIonHeatFlux, SOMLNeutralRecombination, SMOMLIonHeatFlux, 
            SMOMLNeutralRecombination, TEE, SEE, PHLElectronHeatFlux,
            OMLElectronHeatFlux, DTOKSTEE, DTOKSSEE, DTOKSIonHeatFlux,
            DTOKSNeutralRecomb, DTOKSElectronHeatFlux, DUSTTIonHeatFlux };

    if( HeatModels[0] )  HeatTerms.push_back(new Term::EmissivityModel());
    if( HeatModels[1] )  HeatTerms.push_back(new Term::EvaporationModel());
    if( HeatModels[2] )  HeatTerms.push_back(new Term::NewtonCooling());
    if( HeatModels[3] )  HeatTerms.push_back(new Term::NeutralHeatFlux());
    if( HeatModels[11] ) HeatTerms.push_back(new Term::OMLElectronHeatFlux());
    else if( HeatModels[10] ) HeatTerms.push_back(new Term::PHLElectronHeatFlux());
    else if( HeatModels[16] ) HeatTerms.push_back(new Term::DTOKSElectronHeatFlux());
    if( HeatModels[4] )  HeatTerms.push_back(new Term::SOMLIonHeatFlux());
    else if( HeatModels[6] )  HeatTerms.push_back(new Term::SMOMLIonHeatFlux());
    else if( HeatModels[14] ) HeatTerms.push_back(new Term::DTOKSIonHeatFlux());
    else if( HeatModels[17] ) HeatTerms.push_back(new Term::DUSTTIonHeatFlux());
    if( HeatModels[5] )  HeatTerms.push_back(new Term::SOMLNeutralRecombination());
    else if( HeatModels[7] )  HeatTerms.push_back(new Term::SMOMLNeutralRecombination());
    else if( HeatModels[15] ) HeatTerms.push_back(new Term::DTOKSNeutralRecombination());
    if( HeatModels[8] )  HeatTerms.push_back(new Term::SEE());
    else if( HeatModels[12] ) HeatTerms.push_back(new Term::DTOKSSEE());
    if( HeatModels[9] )  HeatTerms.push_back(new Term::TEE());
    else if( HeatModels[13] ) HeatTerms.push_back(new Term::DTOKSTEE());

    std::array<char, CM> ConstModels =
        { EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel,BreakupModel};

    if  (Element == 'W'){ 
        Sample = new Tungsten(Radius,Temp,ConstModels);
    }else if (Element == 'B'){ 
        Sample = new Beryllium(Radius,Temp,ConstModels);
    }else if (Element == 'F'){
        Sample = new Iron(Radius,Temp,ConstModels);
    }else if (Element == 'G'){
        Sample = new Graphite(Radius,Temp,ConstModels);
    }else if (Element == 'D'){
        Sample = new Deuterium(Radius,Temp,ConstModels);
    }else if (Element == 'M'){
        Sample = new Molybdenum(Radius,Temp,ConstModels);
    }else if (Element == 'L'){
        Sample = new Lithium(Radius,Temp,ConstModels);
    }else{ 
        std::cerr << "\nInvalid Option entered";
        return -1;
    }

    double mass = Sample->get_mass();
    Sample->set_potential(Potential);
    std::string filename = "Tests/IntegrationTests/Data/out_" + HeatType + "_Test.txt";
    double SA = Sample->get_surfacearea();
    HeatingModel MyModel(filename,accuracy,HeatTerms,Sample,Pdata);

    MyModel.set_PowerIncident(Power);
    MyModel.UpdateTimeStep();
    MyModel.Vapourise();
    double ModelTime = MyModel.get_totaltime();
    double AnalyticTime = (0), AnalyticFinalTemp = (0);
    if( HeatType == "ConstantHeating" ){
        double p1 = (mass/Power)*((Sample->get_superboilingtemp()-Temp)*Sample->get_heatcapacity());
        double p2 = (mass*Sample->get_latentfusion())/Power;
        double p3 = (mass*Sample->get_latentvapour())/Power;
        std::cout << "\np1 = " << p1 << "\tp2 = " << p2 << "\tp3 = " << p3;
        AnalyticTime = p1 + p2 + p3;
        AnalyticFinalTemp=Sample->get_superboilingtemp();
    }else if( HeatType == "ThermalRadiation" ){
        
        double a = Power;
        double b = Sample->get_emissivity()*SA*Sigma/1000;
    
        double FinalTemp = pow(a/b,1.0/4.0)-0.000000001;
    
        double ti(0), tf(0), t1(0), t2(0), t3(0), t4(0);
        if( FinalTemp >= Sample->get_boilingtemp() ){
            FinalTemp = Sample->get_boilingtemp();
            t4 = Sample->get_latentvapour()*mass/(a-b*pow(Sample->get_boilingtemp(),4));
        }
        if( Element != 'G' ){
            if( FinalTemp < Sample->get_meltingtemp() ){
                tf=(atan(FinalTemp*pow(b/a,1.0/4.0))+atanh(FinalTemp*pow(b/a,1.0/4.0)))
                            /(2*pow(pow(a,3)*b,1.0/4.0));
                ti=(atan(Temp*pow(b/a,1.0/4.0))+atanh(Temp*pow(b/a,1.0/4.0)))
                            /(2*pow(pow(a,3)*b,1.0/4.0));
                t1 = mass*Sample->get_heatcapacity()*(tf-ti);
            }else{
                tf=(atan(Sample->get_meltingtemp()*pow(b/a,1.0/4.0))+atanh(Sample->get_meltingtemp()*pow(b/a,1.0/4.0)))
                            /(2*pow(pow(a,3)*b,1.0/4.0));
                ti=(atan(Temp*pow(b/a,1.0/4.0))+atanh(Temp*pow(b/a,1.0/4.0)))
                            /(2*pow(pow(a,3)*b,1.0/4.0));
                t1 = mass*Sample->get_heatcapacity()*(tf-ti);
        
                t2 = Sample->get_latentfusion()*mass/(a-b*Sample->get_meltingtemp());
    
                tf=(atan(FinalTemp*pow(b/a,1.0/4.0))+atanh(FinalTemp*pow(b/a,1.0/4.0)))
                            /(2*pow(pow(a,3)*b,1.0/4.0));
                ti=(atan(Sample->get_meltingtemp()*pow(b/a,1.0/4.0))+atanh(Sample->get_meltingtemp()*pow(b/a,1.0/4.0)))
                            /(2*pow(pow(a,3)*b,1.0/4.0));   
        
                t3 = mass*Sample->get_heatcapacity()*(tf-ti);
            }
        }else{
            tf=(atan(FinalTemp*pow(b/a,1.0/4.0))+atanh(FinalTemp*pow(b/a,1.0/4.0)))
                        /(2*pow(pow(a,3)*b,1.0/4.0));
            ti=(atan(Temp*pow(b/a,1.0/4.0))+atanh(Temp*pow(b/a,1.0/4.0)))
                        /(2*pow(pow(a,3)*b,1.0/4.0));
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
        }
    
        if( t1 != t1 || t3 != t3 ){
            std::cout << "\nThermal Equilibrium reached before end of test!";
            return -1;
        }
        AnalyticFinalTemp = pow(a/b,1.0/4.0);
        //std::cout << "\n\nt1 = " << t1 << "\nt2 = " << t2 << "\nt3 = " << t3 << "\nt4 = " << t4;
        
        AnalyticTime = (t1 + t2 + t3 + t4)*(1.0/3.0);
        std::cout << "\nWARNING! UNEXPLAINED FACTOR OF 1/3ADDED\n\n";

    }else if( HeatType == "ElectronHeatingRadiation" ){
        double ElectronFlux = Pdata->ElectronDensity*exp(-Potential)*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
        double NeutralFlux = Pdata->NeutralDensity*sqrt(Kb*Pdata->NeutralTemp/(2*PI*Mp));
        
        double ElectronFluxPower = Sample->get_surfacearea()*2*ElectronFlux*Pdata->ElectronTemp*Kb/1000; // Convert from Joules to KJ
        double NeutralFluxPower = Sample->get_surfacearea()*2*NeutralFlux*Pdata->NeutralTemp*Kb/1000; // Convert from Joules to KJ
    
        double a = Power+ElectronFluxPower+NeutralFluxPower;
    
        double b = Sample->get_emissivity()*Sample->get_surfacearea()*Sigma/1000;
        double ti(0), tf(0), t1(0), t2(0), t3(0), t4(0);
        double FinalTemp = pow(a/b,0.25)-0.000000001;
    //  std::cout << "\nFinalTemp = " << FinalTemp; std::cin.get();
        if( FinalTemp < Sample->get_meltingtemp() ){ 
            tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
        }else{
            if( FinalTemp > Sample->get_boilingtemp() )
                FinalTemp = Sample->get_boilingtemp();
            tf=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
            t2 = Sample->get_latentfusion()*mass/(a-b*Sample->get_meltingtemp()); 
            tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            ti=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
    //      std::cout << "\ntf is " << tf;
            t3 = mass*Sample->get_heatcapacity()*(tf-ti);
            if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
                t4 = Sample->get_latentvapour()*mass/(a-b*Sample->get_boilingtemp());
            }
    
        }
    
    
        if( Element == 'G' ){
            if( FinalTemp > Sample->get_boilingtemp() )
            FinalTemp = Sample->get_boilingtemp();
        
            // Only one phase transition
                    tf = (atan(FinalTemp*pow(b/a,0.25)) + atanh(FinalTemp*pow(b/a,0.25)))/(2*pow(pow(a,3)*b,0.25));
            ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))/(2*pow(pow(a,3)*b,0.25));
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
    
            if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
                            t4 = Sample->get_latentvapour()*mass/(a-b*Sample->get_boilingtemp());
                    }
            t2 = 0;
            t3 = 0;
        }
    
//      std::cout << "\nt1 = " << t1 << "\nt2 = " << t2 << "\nt3 = " << t3 << "\nt4 = " << t4; 
    
        if( t1 != t1 && t3 != t3 ){
            std::cout << "\nt1 AND t3 are nan!";
            return -1;
        }else if( t1 != t1 ){
            std::cout << "\nt1 is a nan!";
            return -1;
        }else if( t3 != t3 ){
            std::cout << "\nt3 is a nan!";
            return -1;
        }
        AnalyticFinalTemp = pow(a/b,1.0/4.0);
        AnalyticTime = t1 + t2 + t3 + t4;

    }else if( HeatType == "PlasmaHeatingRadiation" ){
        double ElectronFlux = Pdata->ElectronDensity*exp(-Potential)*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
        double NeutralFlux = Pdata->NeutralDensity*sqrt(Kb*Pdata->NeutralTemp/(2*PI*Mp));
        double IonFlux = ElectronFlux;
    
        double ElectronFluxPower = SA*2*ElectronFlux*Pdata->ElectronTemp*Kb/1000; // Convert from Joules to KJ
        double NeutralFluxPower = SA*2*NeutralFlux*Pdata->NeutralTemp*Kb/1000; // Convert from Joules to KJ
        double IonFluxPower = (SA*IonFlux*Pdata->IonTemp*Kb/1000) // Convert from Joules to KJ
        *(2+2*Potential*(Pdata->ElectronTemp/Pdata->IonTemp)+pow(Potential*(Pdata->ElectronTemp/Pdata->IonTemp),2))/(1+Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
    
        double a = Power+ElectronFluxPower+NeutralFluxPower+IonFluxPower;
    
        double b = Sample->get_emissivity()*SA*Sigma/1000;
        double ti(0), tf(0), t1(0), t2(0), t3(0), t4(0);
        double FinalTemp = pow(a/b,0.25)-0.000000001;
    //  std::cout << "\nFinalTemp = " << FinalTemp; std::cin.get();
        if( FinalTemp < Sample->get_meltingtemp() ){ 
            tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
        }else{
            if( FinalTemp > Sample->get_boilingtemp() )
                FinalTemp = Sample->get_boilingtemp();
            tf=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            ti=(atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
            t2 = Sample->get_latentfusion()*mass/(a-b*Sample->get_meltingtemp()); 
            tf=(atan(FinalTemp*pow(b/a,0.25))+atanh(FinalTemp*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
            ti=(atan(Sample->get_meltingtemp()*pow(b/a,0.25))+atanh(Sample->get_meltingtemp()*pow(b/a,0.25)))
                        /(2*pow(pow(a,3)*b,0.25));
    //      std::cout << "\ntf is " << tf;
            t3 = mass*Sample->get_heatcapacity()*(tf-ti);
            if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
                t4 = Sample->get_latentvapour()*mass/(a-b*Sample->get_boilingtemp());
            }
    
            std::cout << "\nt1 = " << t1 << "\nt2 = " << t2 << "\nt3 = " << t3 << "\nt4 = " << t4; 
        }
    
        if( Element == 'G' ){
            if( FinalTemp > Sample->get_boilingtemp() )
                FinalTemp = Sample->get_boilingtemp();
        
            // Only one phase transition
            //double p1 = atan(FinalTemp*pow(b/a,0.25));
            //double p2 = atanh(FinalTemp*pow(b/a,0.25));
            //double p3 = 2*pow(pow(a,3)*b,0.25);
                    tf = (atan(FinalTemp*pow(b/a,0.25)) + atanh(FinalTemp*pow(b/a,0.25)))/(2*pow(pow(a,3)*b,0.25));
            ti = (atan(Temp*pow(b/a,0.25))+atanh(Temp*pow(b/a,0.25)))/(2*pow(pow(a,3)*b,0.25));
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
    
            if( FinalTemp >= Sample->get_boilingtemp() ){
                            t4 = Sample->get_latentvapour()*mass/(a-b*Sample->get_boilingtemp());
            }
            t2 = 0;
            t3 = 0;
    //      std::cout << "\nt1 = " << t1 << "\nt2 = " << t2 << "\nt3 = " << t3 << "\nt4 = " << t4; 
    //      std::cout << "\ntf = " << tf << "\nti = " << ti << "\nFinalTemp = " << FinalTemp; 
        }
    
        if( t1 != t1 && t3 != t3 ){
            std::cout << "\nt1 AND t3 are nan!";
            return -1;
        }else if( t1 != t1 ){
            std::cout << "\nt1 is a nan!";
            return -1;
        }else if( t3 != t3 ){
            std::cout << "\nt3 is a nan!";
            return -1;
        }
        AnalyticFinalTemp = pow(a/b,1.0/4.0);
        AnalyticTime = t1 + t2 + t3 + t4;

    }else if( HeatType == "PlasmaHeatingNeutralRecomb" ){
        double ElectronFlux = Pdata->ElectronDensity*exp(-Potential)*sqrt(Kb*Pdata->ElectronTemp/(2*PI*Me));
        double NeutralFlux = Pdata->NeutralDensity*sqrt(Kb*Pdata->NeutralTemp/(2*PI*Mp));
        double IonFlux = ElectronFlux;
    
        double ElectronFluxPower = Sample->get_surfacearea()*2*ElectronFlux*Pdata->ElectronTemp*Kb/1000; // Convert from Joules to KJ
        double NeutralFluxPower = Sample->get_surfacearea()*2*NeutralFlux*Pdata->NeutralTemp*Kb/1000; // Convert from Joules to KJ
        double IonFluxPower = (Sample->get_surfacearea()*IonFlux*Pdata->IonTemp*Kb/1000) // Convert from Joules to KJ
        *(2+2*Potential*(Pdata->ElectronTemp/Pdata->IonTemp)+pow(Potential*(Pdata->ElectronTemp/Pdata->IonTemp),2))
        /(1+Potential*(Pdata->ElectronTemp/Pdata->IonTemp));
        double NeutralRecombPower = Sample->get_surfacearea()*14.7*echarge*IonFlux/1000; // Convert from J to kJ
    
    //  std::cout << "\nElectronFluxPower = " << ElectronFluxPower << "\nNeutralFluxPower = " << NeutralFluxPower << "\nIonFluxPower = " << IonFluxPower << "\nNeutralRecombPower = " << NeutralRecombPower;
    
        double a = Power+ElectronFluxPower+NeutralFluxPower+IonFluxPower+NeutralRecombPower;
        double b = Sample->get_surfacearea()*2*Kb*NeutralFlux/1000; // Convert from J to kJ
    
        double ti(0), tf(0), t1(0), t2(0), t3(0), t4(0);
        double FinalTemp = a/b;
    //  std::cout << "\nFinalTemp = " << FinalTemp; std::cin.get();
        if( FinalTemp < Sample->get_meltingtemp() ){ 
            tf=-log(a-b*FinalTemp)/b;
            ti=-log(a-b*Temp)/b;
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
        }else{
            if( FinalTemp > Sample->get_boilingtemp() )
                FinalTemp = Sample->get_boilingtemp()-0.000000001;
            tf=-log(a-b*Sample->get_meltingtemp())/b;
            ti=-log(a-b*Temp)/b;
    
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
            t2 = Sample->get_latentfusion()*mass/(a-b*Sample->get_meltingtemp()); 
            tf=-log(a-b*FinalTemp)/b;
            ti=-log(a-b*Sample->get_meltingtemp())/b;
    
            t3 = mass*Sample->get_heatcapacity()*(tf-ti);
            if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
                t4 = Sample->get_latentvapour()*mass/(a-b*Sample->get_boilingtemp());
            }
    
        }
    
        if( Element == 'G' ){
            if( FinalTemp > Sample->get_boilingtemp() )
            FinalTemp = Sample->get_boilingtemp();
        
            // Only one phase transition
            tf=-log(a-b*FinalTemp)/b;
            ti=-log(a-b*Temp)/b;
            t1 = mass*Sample->get_heatcapacity()*(tf-ti);
    
            if( Sample->get_temperature() >= Sample->get_boilingtemp() ){
                            t4 = Sample->get_latentvapour()*mass/(a-b*Sample->get_boilingtemp());
                    }
            t2 = 0;
            t3 = 0;
    //      std::cout << "\nt1 = " << t1 << "\nt2 = " << t2 << "\nt3 = " << t3 << "\nt4 = " << t4; 
    
        }
    
        if( t1 != t1 && t3 != t3 ){
            std::cout << "\nt1 AND t3 are nan!";
            return -1;
        }else if( t1 != t1 ){
            std::cout << "\nt1 is a nan!";
            return -1;
        }else if( t3 != t3 ){
            std::cout << "\nt3 is a nan!";
            return -1;
        }
        AnalyticFinalTemp = pow(a/b,1.0/4.0);
        AnalyticTime = t1 + t2 + t3 + t4;

    }   

    double ReturnVal = 0;
    if( ModelTime == AnalyticTime )             ReturnVal = 1;
    else if( fabs(1-ModelTime/AnalyticTime) < 0.01 )    ReturnVal = 2;
    else                            ReturnVal = -1;

    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
    std::cout << "\n\n*****\nModelFinalTemperature = " << Sample->get_temperature() << "K : AnalyticFinalTemperature = " << AnalyticFinalTemp << "K";
    std::cout << "\n\n*****\nModelTime = " << ModelTime << "s : AnalyticTime = " << AnalyticTime << "s";
    std::cout << "\nPercentage Deviation = " << fabs(100-100*ModelTime/AnalyticTime) <<"%\n*****\n\n";

    delete Sample;

    return ReturnVal;
}
