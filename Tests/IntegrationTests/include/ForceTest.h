#include "ForceModel.h"

int ForceTest(char Element, std::string ForceType, double accuracy){
    clock_t begin = clock();
    // ********************************************************** //
    // FIRST, define program default behaviour

    // Define the behaviour of the models for the temperature dependant constants, the time step and the 'Name' variable.
    char EmissivityModel = 'c';     // Possible values 'c', 'v' and 'f': Corresponding to (c)onstant, (v)ariable and from (f)ile
    char ExpansionModel = 'c';     // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable, (s)et 
                                                    // and (z)ero expansion
    char HeatCapacityModel = 'c';     // Possible values 'c', 'v' and 's': Corresponding to (c)onstant, (v)ariable and (s)et
    char BoilingModel = 'y';     // Possible values 'y', 'n', 's' and 't': Corresponding to (y)es, (n)o, (s)uper 
                                                    // and (t)homson
    char BreakupModel = 'n';    // Possible values 'r', 'e', 'b'  and 'n': Corresponding to (r)otational, (e)lectrostatic, (b)oth and (n)o

     // Parameters describing the heating model
    double Radius=5e-8;         // m
    double Temp=280;        // K
    double Potential = 2.5;        // Normalised Potential
    Matter *Sample;            // Define the sample matter type

    // Set to true all Force models that are wanted
    bool Gravity = false;        // IS OFF!
    bool Lorentz = false;        // IS OFF!
    bool SOMLIonDrag = false;    // IS OFF!
    bool SMOMLIonDrag = false;    // IS OFF!
    bool HybridIonDrag = false;    // IS OFF!
    bool DTOKSIonDrag = false;    // IS OFF!
    bool DUSTTIonDrag = false;    // IS OFF!
    bool NeutralDrag = false;    // Is OFF!
    bool RocketForce = false;    // Is OFF!

    PlasmaData Pdata;

//    threevector PlasmaVelocity(10, 5, 3); // Taken from initial for DTOKS
//    Pdata.PlasmaVel = PlasmaVelocity;
    Pdata.ElectronTemp = 10*1.16e4;    // K, Electron Temperature, convert from eV
    threevector GravityForce(0, 0, -9.81);
    Pdata.Gravity = GravityForce;
//    threevector Efield(1.0, -2.0, 3.0);
    Pdata.mi           = Mp;        // kg, Ion Mass
    Pdata.Z           = 1.0;        // arb, Ion Charge
    Pdata.A            = 1.0;       // arb, Atomic Number
    threevector Efield(1.0, 0.0, 0.0);
    Pdata.ElectricField = Efield;
    threevector Bfield(0.0, 0.0, 10.0);
    Pdata.MagneticField = Bfield;

    if( ForceType == "Gravity" ){
        Gravity = true;
    }else if( ForceType == "EFieldGravity" ){
        Gravity = true;
        Lorentz = true;
        threevector zeros(0.0, 0.0, 0.0);
        Pdata.MagneticField = zeros;
    }else if( ForceType == "BField" ){
        Lorentz = true;
        threevector zeros(0.0, 0.0, 0.0);
        Pdata.ElectricField = zeros;
    }else if( ForceType == "Lorentz" ){
        Lorentz = true;
    }else if( ForceType == "LorentzGravity" ){
        Gravity = true;
        Lorentz = true;
        threevector Bfield(0.0, 1.0, 0.0);
    }else{
        std::cout << "\n\nInput not recognised! Exiting program\n.";
        return -1;
    }
    std::array<bool,9> ForceModels  = {Gravity,Lorentz,SOMLIonDrag,SMOMLIonDrag,
        HybridIonDrag,DTOKSIonDrag,DUSTTIonDrag,NeutralDrag,RocketForce};
    std::vector<ForceTerm*> ForceTerms;
    if( ForceModels[0] ) ForceTerms.push_back(new Term::Gravity());
    if( ForceModels[1] ) ForceTerms.push_back(new Term::LorentzForce());
    if( ForceModels[2] ) ForceTerms.push_back(new Term::SOMLIonDrag());
    else if( ForceModels[3] ) ForceTerms.push_back(new Term::SMOMLIonDrag());
    else if( ForceModels[4] ) ForceTerms.push_back(new Term::HybridIonDrag());
    else if( ForceModels[5] ) ForceTerms.push_back(new Term::DTOKSIonDrag());
    else if( ForceModels[6] ) ForceTerms.push_back(new Term::DUSTTIonDrag());
    if( ForceModels[7] ) ForceTerms.push_back(new Term::NeutralDrag());
    if( ForceModels[8] ) ForceTerms.push_back(new Term::RocketForce());

    // Models and ConstModels are placed in an array in this order:
    std::array<char, CM> ConstModels =
        { EmissivityModel,ExpansionModel,HeatCapacityModel,BoilingModel,BreakupModel};

    if    (Element == 'W'){ 
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
    Sample->set_potential(Potential);
    // Determine initial conditions depending on the force type
    threevector vinit(0.0,0.0,1.0);
    threevector rinit(10.0,0.0,1.0);

    if( ForceType == "BField" || ForceType == "Lorentz"
	    || ForceType == "LorentzGravity" ){ // Start at large radius and with positive radial and z velocity
        vinit.sety(1.0);
//        vinit.setz(1.0);
//        rinit.setx(10.0);
    }

    Sample->update_motion(rinit,vinit,0.0);

    // START NUMERICAL TESTING
    std::string filepath = "Tests/IntegrationTests/Data/out_" + ForceType + "_" + Element + "_Test.txt";
    ForceModel MyModel(filepath,accuracy,ForceTerms,Sample,Pdata);

    double Mass = Sample->get_mass();
    MyModel.UpdateTimeStep();
    double ModelTimeStep = MyModel.get_timestep();
	double TimeMax = 0.1;
	double TotalTime(0);
    for( TotalTime = 0; TotalTime < TimeMax; TotalTime+=ModelTimeStep){
        MyModel.Force();
    }
    threevector ModelVelocity = Sample->get_velocity();
    // END NUMERICAL TESTING
    
	// START ANALYTICAL MODEL
    threevector AnalyticVelocity;
    double qtom = -4*PI*epsilon0*Sample->get_radius()*Kb*Pdata.ElectronTemp*Sample->get_potential()/(echarge*Mass);

    if( ForceType == "Gravity" ){
        AnalyticVelocity = (Pdata.Gravity)*(TotalTime)+vinit;
    }else if( ForceType == "EFieldGravity" ){
        AnalyticVelocity = (Efield*qtom + Pdata.Gravity)*(TotalTime)+vinit;
    }else if( ForceType == "BField"){
        
        double GyroFrequency = qtom*Pdata.MagneticField.mag3();
        double VPerp = sqrt(vinit.getx()*vinit.getx()+vinit.gety()*vinit.gety());
    
        double vx = VPerp*sin(GyroFrequency*TotalTime);
        double vy = VPerp*cos(GyroFrequency*TotalTime);
        double vz = vinit.getz();
        AnalyticVelocity.setx(vx);
        AnalyticVelocity.sety(vy);
        AnalyticVelocity.setz(vz);
    }else if( ForceType == "Lorentz" ){
        // Starting from zero initial velocity with EField = (1.0,0.0,0.0), BField = (0.0,0.0,10.0)
        // We have an ExB drift vE, and gyro-motion with gyro-velocity equal to the max of vE
        double GyroFrequency = qtom*Pdata.MagneticField.mag3(); // Direction of this as vector is parallel to B
        double VPerp = sqrt(vinit.getx()*vinit.getx()+vinit.gety()*vinit.gety());
    
        threevector vE = (Pdata.ElectricField^Pdata.MagneticField)*(1/(Pdata.MagneticField*Pdata.MagneticField));    
    
        double vx = vE.getx()+(Pdata.ElectricField.getx()/Pdata.MagneticField.mag3())*sin(GyroFrequency*TotalTime); // vE.getx()
        double vy = vE.gety()+(Pdata.ElectricField.getx()/Pdata.MagneticField.mag3())*cos(GyroFrequency*TotalTime); // vE.gety()
        double vz = qtom*Efield.getz()*TotalTime+vinit.getz();
    
        AnalyticVelocity.setx(vx);
        AnalyticVelocity.sety(vy);
        AnalyticVelocity.setz(vz);
    }else if( ForceType == "LorentzGravity" ){
        // Starting from zero initial velocity with 
        // EField = (1.0,0.0,0.0), BField = (0.0,0.0,10.0) & Gravity = (0.0,0.0,-9.81)
        // We have an ExB drift vE, and Gravity force drift vG and gyro-motion with gyro-velocity equal to the max of vE+vG
        double GyroFrequency = qtom*Pdata.MagneticField.mag3(); // Direction of this as vector is parallel to B
        double VPerp = sqrt(vinit.getx()*vinit.getx()+vinit.gety()*vinit.gety());
    
        threevector vE = (Pdata.ElectricField^Pdata.MagneticField)*(1/(Pdata.MagneticField*Pdata.MagneticField));    
        threevector vG = qtom*(Pdata.Gravity^Pdata.MagneticField)*(1/(Pdata.MagneticField*Pdata.MagneticField));    
    
        double vdriftMag = (vE +vG).mag3();
    
        double vx = vE.getx()+vdriftMag*sin(GyroFrequency*TotalTime); // vE.getx()
        double vy = vE.gety()+vdriftMag*cos(GyroFrequency*TotalTime); // vE.gety()
        double vz = (qtom*Efield.getz()+Pdata.Gravity.getz())*TotalTime+vinit.getz();
    
        AnalyticVelocity.setx(vx);
        AnalyticVelocity.sety(vy);
        AnalyticVelocity.setz(vz);

    }

    // END ANALYTICAL MODEL
    
    double ReturnVal = 0;

    if( (ModelVelocity - AnalyticVelocity).mag3() == 0 )                 ReturnVal = 1;
    else if( fabs(1.0-ModelVelocity.mag3()/AnalyticVelocity.mag3()) < 0.0001 )    ReturnVal = 3;
    else if( fabs(1.0-ModelVelocity.mag3()/AnalyticVelocity.mag3()) < 0.1 )    ReturnVal = 2;
    else                                        ReturnVal = -1;

    clock_t end = clock();
    double elapsd_secs = double(end-begin)/CLOCKS_PER_SEC;
        
    std::cout << "\n\n*****\n\nIntegrationTest 1 completed in " << elapsd_secs << "s\n";
    std::cout << "\n\n*****\nModelVel = " << ModelVelocity << "m s^-1 : AnalyticVel = " << AnalyticVelocity << "m s^-1";

    if( AnalyticVelocity.mag3() != 0.0 )
        std::cout << "\nPercentage Deviation: Mag = " << 100*fabs(1.0-ModelVelocity.mag3()/AnalyticVelocity.mag3()) <<"%\n*****\n\n";
    if( AnalyticVelocity.getx() != 0.0 )
        std::cout << "\nPercentage Deviation: Xdir = " << 100*fabs(1.0-ModelVelocity.getx()/AnalyticVelocity.getx()) <<"%\n*****\n\n";
    if( AnalyticVelocity.gety() != 0.0 )
    std::cout << "\nPercentage Deviation: Ydir = " << 100*fabs(1.0-ModelVelocity.gety()/AnalyticVelocity.gety()) <<"%\n*****\n\n";
    if( AnalyticVelocity.getz() != 0.0 )
        std::cout << "\nPercentage Deviation: Zdir = " << 100*fabs(1.0-ModelVelocity.getz()/AnalyticVelocity.getz()) <<"%\n*****\n\n";

    return ReturnVal;
}
/*
    threevector Eparr(Bfield.getx()*Efield.getx(),Bfield.gety()*Efield.gety(),Bfield.getz()*Efield.getz());
    threevector Eperp = Efield - Eparr;
    double ConvertKelvsToeV(8.621738e-5);
    // NOTE: this is the charge to mass ratio for a negative grain only...
    double qtom = -3.0*epsilon0*Pdata.ElectronTemp*ConvertKelvsToeV*Sample->get_potential()
            /(pow(Sample->get_radius(),2)*Sample->get_density());
    double AngularVel = qtom*Bfield.mag3();
    double rc = 0.0;
    threevector Vperp(rc*AngularVel*cos(AngularVel*TotalTime),0.0,rc*AngularVel*sin(AngularVel*TotalTime));
//    threevector CrossParts = (GravityForce^Bfield + qtom*Eperp^Bfield)*(1.0/ sqrt(pow(Bfield.mag3(),2)))*TotalTime;
    threevector CrossParts = (GravityForce + qtom*Eperp)*TotalTime;
    threevector EPart = Eparr*qtom*TotalTime;
    
    std::cout << "\nLorentzForce = " << (qtom*(Efield+ModelVelocity^Bfield));

*/
