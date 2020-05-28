/** @file DTOKSU.cpp
 *  @brief Implementation of class for physics models relevant to dust charging
 *  
 *  Implement the member functions of the DTOKSU class
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug bugs, they definitely exist
 */

#include "DTOKSU.h"

DTOKSU::DTOKSU( std::array<float,MN> acclvls, Matter *& sample, PlasmaData 
&pdata, std::vector<HeatTerm*> HeatTerms, std::vector<ForceTerm*> ForceTerms, 
std::vector<CurrentTerm*> CurrentTerms): 
Sample(sample), WallBound(BoundaryDefaults), CoreBound(BoundaryDefaults),
HM("Data/default_hm_0.txt",acclvls[1],HeatTerms,sample,pdata),
FM("Data/default_fm_0.txt",acclvls[2],ForceTerms,sample,pdata),
CM("Data/default_cm_0.txt",acclvls[0],CurrentTerms,sample,pdata){
    D_Debug("\n\nIn DTOKSU::DTOKSU( std::array<float,MN> acclvls, "
        << "Matter *& sample, PlasmaData &pdata, "
        << "std::vector<HeatTerm*> HeatTerms, "
        << "std::vector<ForceTerm*> ForceTerms, "
        << "std::vector<CurrentTerm*> CurrentTerms): "
        << "Sample(sample), WallBound(BoundaryDefaults), "
        << "CoreBound(BoundaryDefaults),"
        << "HM(\"Data/default_hm_0.txt\",acclvls[1],heatmodels,sample,pdata),"
        << "FM(\"Data/default_fm_0.txt\",acclvls[2],forcemodels,sample,pdata),"
        << "CM(\"Data/default_cm_0.txt\",acclvls[0],chargemodels,sample,pdata)"
        << "\n\n");
    D_Debug("\n\n******************* SETUP FINISHED ******************* \n\n");

    MaxTime = 0.5;
    TotalTime = 0;
    ReflectedLastStep = false;
    create_file("Data/df.txt");
}

DTOKSU::DTOKSU( std::array<float,MN> acclvls, Matter *& sample,
PlasmaGrid_Data &pgrid,std::vector<HeatTerm*> HeatTerms, 
std::vector<ForceTerm*> ForceTerms, std::vector<CurrentTerm*> CurrentTerms):
Sample(sample), WallBound(BoundaryDefaults), CoreBound(BoundaryDefaults),
HM("Data/default_hm_0.txt",acclvls[1],HeatTerms,sample,pgrid),
FM("Data/default_fm_0.txt",acclvls[2],ForceTerms,sample,pgrid),
CM("Data/default_cm_0.txt",acclvls[0],CurrentTerms,sample,pgrid){
    D_Debug("\n\nIn DTOKSU::DTOKSU( std::array<float,MN> acclvls, "
        << "Matter *& sample, PlasmaGrid_Data &pgrid, "
        << "std::vector<HeatTerm*> HeatTerms, "
        << "std::vector<ForceTerm*> ForceTerms, "
        << "std::vector<CurrentTerm*> CurrentTerms): "
        << "Sample(sample), WallBound(BoundaryDefaults), "
        << "CoreBound(BoundaryDefaults),"
        << "HM(\"Data/default_hm_0.txt\",acclvls[1],heatmodels,sample,pdata),"
        << "FM(\"Data/default_fm_0.txt\",acclvls[2],forcemodels,sample,pdata),"
        << "CM(\"Data/default_cm_0.txt\",acclvls[0],chargemodels,sample,pdata)"
        << "\n\n");
    D_Debug("\n\n******************* SETUP FINISHED ******************* \n\n");

    MaxTime = 0.5;
    TotalTime = 0;
    ReflectedLastStep = false;
    create_file("Data/df.txt");
}

DTOKSU::DTOKSU( std::array<float,MN> acclvls, Matter *& sample, 
PlasmaGrid_Data &pgrid, PlasmaData &pdata, std::vector<HeatTerm*> HeatTerms, 
std::vector<ForceTerm*> ForceTerms, std::vector<CurrentTerm*> CurrentTerms): 
Sample(sample), WallBound(BoundaryDefaults), CoreBound(BoundaryDefaults),
HM("Data/default_hm_0.txt",acclvls[1],HeatTerms,sample,pgrid,pdata),
FM("Data/default_fm_0.txt",acclvls[2],ForceTerms,sample,pgrid,pdata),
CM("Data/default_cm_0.txt",acclvls[0],CurrentTerms,sample,pgrid,pdata){
    D_Debug("\n\nIn DTOKSU::DTOKSU( std::array<float,MN> acclvls, "
        << "Matter *& sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata,"
        << "std::vector<HeatTerm*> HeatTerms, "
        << "std::vector<ForceTerm*> ForceTerms, "
        << "std::vector<CurrentTerm*> CurrentTerms): "
        << "Sample(sample), WallBound(BoundaryDefaults), "
        << "CoreBound(BoundaryDefaults),"
        << "HM(\"Data/default_hm_0.txt\",acclvls[1],heatmodels,sample,pdata),"
        << "FM(\"Data/default_fm_0.txt\",acclvls[2],forcemodels,sample,pdata),"
        << "CM(\"Data/default_cm_0.txt\",acclvls[0],chargemodels,sample,pdata)"
        << "\n\n");
    D_Debug("\n\n******************* SETUP FINISHED ******************* \n\n");

    MaxTime = 0.5;
    TotalTime = 0;
    ReflectedLastStep = false;
    create_file("Data/df.txt");
}

DTOKSU::DTOKSU( std::array<float,MN> acclvls, Matter *& sample, 
PlasmaGrid_Data &pgrid, PlasmaData &pdata, Boundary_Data &wbound, 
Boundary_Data &cbound, std::vector<HeatTerm*> HeatTerms, 
std::vector<ForceTerm*> ForceTerms, std::vector<CurrentTerm*> CurrentTerms): 
Sample(sample), WallBound(wbound), CoreBound(cbound),
HM("Data/default_hm_0.txt",acclvls[1],HeatTerms,sample,pgrid,pdata),
FM("Data/default_fm_0.txt",acclvls[2],ForceTerms,sample,pgrid,pdata),
CM("Data/default_cm_0.txt",acclvls[0],CurrentTerms,sample,pgrid,pdata){
    D_Debug("\n\nIn DTOKSU::DTOKSU( std::array<float,MN> acclvls, "
        << "Matter *& sample, PlasmaGrid_Data &pgrid, PlasmaData &pdata,"
        << "Boundary_Data &wbound, Boundary_Data &cbound,"
        << "std::vector<HeatTerm*> HeatTerms, "
        << "std::vector<ForceTerm*> ForceTerms, "
        << "std::vector<CurrentTerm*> CurrentTerms): "
        << "Sample(sample), WallBound(BoundaryDefaults), "
        << "CoreBound(BoundaryDefaults),"
        << "HM(\"Data/default_hm_0.txt\",acclvls[1],heatmodels,sample,pdata),"
        << "FM(\"Data/default_fm_0.txt\",acclvls[2],forcemodels,sample,pdata),"
        << "CM(\"Data/default_cm_0.txt\",acclvls[0],chargemodels,sample,pdata)"
        << "\n\n");
    D_Debug("\n\n******************* SETUP FINISHED ******************* \n\n");

    MaxTime = 0.5;
    TotalTime = 0;
    ReflectedLastStep = false;
    create_file("Data/df.txt");
}

void DTOKSU::create_file( std::string filename ){
    D_Debug("\n\nIn DTOKSU::create_file(std::string filename)\n\n");
    MyFile.open(filename);
    MyFile << "TotalTime\n";
}

void DTOKSU::OpenFiles( std::string filename, unsigned int i ){
    D_Debug("\n\nIn DTOKSU::OpenFiles(std::string filename,"
        << " unsigned int i)\n\n");
    create_file(filename + "_df_" + std::to_string(i) + ".txt");
    HM.CreateFile(filename + "_hm_" + std::to_string(i) + ".txt",false);
    FM.CreateFile(filename + "_fm_" + std::to_string(i) + ".txt");
    CM.CreateFile(filename + "_cm_" + std::to_string(i) + ".txt");
}

void DTOKSU::CloseFiles(){
    D_Debug("\n\nIn DTOKSU::CloseFiles()\n\n");
    MyFile.close();
    HM.close_file();
    FM.close_file();
    CM.close_file();
}

void DTOKSU::ResetModelTime(double HMTime, double FMTime, double CMTime){
    D_Debug("\n\nIn DTOKSU::ResetModelTime(double HMTime, double FMTime, "
        << "double CMTime)\n\n");
    HM.AddTime(HMTime);
    FM.AddTime(FMTime);
    CM.AddTime(CMTime);
}

void DTOKSU::print(){
    D_Debug("\tIn DTOKSU::print()\n\n");
    MyFile  << TotalTime;
    MyFile << "\n";
}

void DTOKSU::SpecularReflection(){
    D_Debug("\tIn DTOKSU::SpecularReflection()\n\n");
    double x = Sample->get_position().getx();
    double y = Sample->get_position().getz();

    double MinDist = 100.0;
    int MinIndex(0),SecondMinIndex(0);
    for (int i=0; i<WallBound.Grid_Pos.size(); i++) {
        //!< Calculate distance to point i
        double DistanceToPoint = sqrt((WallBound.Grid_Pos[i].first-x)*
            (WallBound.Grid_Pos[i].first-x)+
            (WallBound.Grid_Pos[i].second-y)*(WallBound.Grid_Pos[i].second-y));
        //!< if this is the smallest distance so far
        if( DistanceToPoint < MinDist ){ 
            MinIndex = i; // save index of nearest point
            MinDist = DistanceToPoint;
        }
//        std::cout << "\n" << i << "\t" << WallBound.Grid_Pos[i].first << "\t" 
//          << WallBound.Grid_Pos[i].first;
    }

    //!< If MinIndex is last element of array
    int MinIndexCopy = MinIndex;
    if( MinIndex == (WallBound.Grid_Pos.size()-1) ) 
        MinIndexCopy = -1; //!< Go back to start of the array
    double Dist1 = sqrt((WallBound.Grid_Pos[MinIndexCopy+1].first-x)*
        (WallBound.Grid_Pos[MinIndexCopy+1].first-x)+
        (WallBound.Grid_Pos[MinIndexCopy+1].second-y)*
        (WallBound.Grid_Pos[MinIndexCopy+1].second-y));
    if( MinIndex == 0 ) //!< If MinIndex is first element of array,
        MinIndexCopy = WallBound.Grid_Pos.size(); //!< Go to end of array
    double Dist2 = sqrt((WallBound.Grid_Pos[MinIndexCopy-1].first-x)*
        (WallBound.Grid_Pos[MinIndexCopy-1].first-x)+
        (WallBound.Grid_Pos[MinIndexCopy-1].second-y)*
        (WallBound.Grid_Pos[MinIndexCopy-1].second-y));

    //!< Determine which of Dist1 or Dist2 is closer  
    if( Dist1 < Dist2 )
        if( MinIndex == (WallBound.Grid_Pos.size()-1) )
            SecondMinIndex = 0;
        else
            SecondMinIndex = MinIndex+1;
    else
        if( MinIndex == 0 )
            SecondMinIndex = WallBound.Grid_Pos.size()-1;
        else
            SecondMinIndex = MinIndex-1;

    double dr = WallBound.Grid_Pos[MinIndex].first-
        WallBound.Grid_Pos[SecondMinIndex].first;
    double dz = WallBound.Grid_Pos[MinIndex].second-
        WallBound.Grid_Pos[SecondMinIndex].second;

    threevector normal(dz,0.0,-1.0*dr);
    threevector Zeroes(0.0,0.0,0.0);
    threevector SpecularReflDir = -2.0*(normal.getunit()*
        Sample->get_velocity())*normal.getunit()+Sample->get_velocity();
    threevector ReflectedVel(SpecularReflDir.getx(),
        Sample->get_velocity().gety(),SpecularReflDir.getz());
    
//    std::cout << "\n" << MinIndex << "\t" << SecondMinIndex;
//    std::cout << "\n" << WallBound.Grid_Pos[MinIndex].first << "\t" 
//        <<  WallBound.Grid_Pos[MinIndex].second;
//    std::cout << "\t" << x << "\t" <<  y << "\t" << normal << "\t" 
//     << Sample->get_velocity() << "\t" 
//     << normal.getunit()*Sample->get_velocity() << "\t" << SpecularReflDir 
//     << "\t" << ReflectedVel;
//    std::cout << "\n" << WallBound.Grid_Pos[SecondMinIndex].first << "\t" 
//     <<  WallBound.Grid_Pos[SecondMinIndex].second;
    std::cout << "\n\nCollision with Wall!";

    Sample->update_motion(Zeroes,ReflectedVel-Sample->get_velocity(),0.0);
}

//!< http://alienryderflex.com/polygon/
bool DTOKSU::Boundary_Check(bool InOrOut){
    D_Debug("\tIn DTOKSU::Boundary_Check(bool InOrOut)\n\n");
    Boundary_Data Edge;
    //!< Determine if it's core or wall boundary
    if( InOrOut ){
        Edge = CoreBound;
    }else{
        Edge = WallBound;
    }
    
    assert(Edge.Grid_Pos.size() > 2);
    int j=Edge.Grid_Pos.size()-1 ;
    bool Inside=false;
    double x = Sample->get_position().getx();
    double y = Sample->get_position().getz();

    for (int i=0; i<Edge.Grid_Pos.size(); i++) {
        if (((Edge.Grid_Pos[i].second < y && Edge.Grid_Pos[j].second >= y)
            ||   (Edge.Grid_Pos[j].second < y && Edge.Grid_Pos[i].second >= y))
            &&  (Edge.Grid_Pos[i].first <= x || Edge.Grid_Pos[j].first <= x)) {

            Inside^=(Edge.Grid_Pos[i].first+(y-Edge.Grid_Pos[i].second)
                    /(Edge.Grid_Pos[j].second-Edge.Grid_Pos[i].second)
                    *(Edge.Grid_Pos[j].first-Edge.Grid_Pos[i].first)<x);
        }
        j=i; 
    }

    //!< In this case, it's a wall and the particle has gone through it,
    //!< Here we implement specular refulection
    if( InOrOut ){
        return Inside;
    }else{
        if( !Inside && !ReflectedLastStep ){
            SpecularReflection();
            ReflectedLastStep = true;
            return Inside; //!< Pretend we're still inside
        }else if( !Inside && ReflectedLastStep ){
            return Inside; //!< Pretend we're still inside
        }else if( Inside ){
            ReflectedLastStep = false;
        }
        return !Inside;
    }
}

void DTOKSU::ImpurityPrint(){
    HM.ImpurityPrint();
}

int DTOKSU::Run(){
    D_Debug("- In DTOKSU::Run()\n\n");

    double HeatTime(0),ForceTime(0),ChargeTime(0);

    //!< Update the plasma data from the plasma grid for all models...
    //!< This has to be done individually because PlasmaData is shared accross 
    //!< models in current design
    D_Debug("\n\n***** DTOKSU::Run() :: update_plasmadata() *****\n\n");
    bool cm_InGrid = CM.update_plasmadata(); 
    bool hm_InGrid = HM.update_plasmadata(); 
    bool fm_InGrid = FM.update_plasmadata(); 
    assert( cm_InGrid == fm_InGrid 
        && fm_InGrid == hm_InGrid 
        && hm_InGrid == cm_InGrid );
    //!< Charge instantaneously as soon as we start, have to add time though...
    CM.Charge(1e-100);
    //!< Need to manually update the first time as first step is not necessarily
    //!< heating      
    Sample->update();
    bool ErrorFlag(false);
    while( cm_InGrid && !Sample->is_split() ){

        // ***** START OF : DETERMINE TIMESCALES OF PROCESSES ***** //  
        //!< Charge instantaneously as soon as we start, have to add a time 
        //!< though...
        CM.Charge(1e-100);
        D_Debug("\n\n***** DTOKSU::Run() :: get Model timescales *****\n\n");
        ChargeTime  = CM.UpdateTimeStep(); //!< Check Time step length is good
        ForceTime   = FM.UpdateTimeStep(); //!< Check Time step length is good
        HeatTime    = HM.UpdateTimeStep(); //!< Check Time step length is good
        if( HeatTime == 1) break; //!< Thermal Equilibrium Reached

        HM.UpdateRERN();
        //!< We will assume Charging Time scale is much faster than either 
        //!< heating or moving, but check for the other case.
        double MaxTimeStep = std::max(ForceTime,HeatTime);
        double MinTimeStep = std::min(ForceTime,HeatTime);

        //!< Check Charging timescale isn't the fastest timescale.
        if( ChargeTime > MinTimeStep && ChargeTime != 1){
            static bool runOnce = true;
            std::string Warning = "*** Charging Time scale is not the shortest";
            Warning += " timescale!! ***\n";
            WarnOnce(runOnce,Warning);
            std::cout << "\nChargeTime = " << ChargeTime 
                << "\t:\tMinTime = " << MinTimeStep;
            ErrorFlag = true;
//          break;
        }

        // ***** END OF : DETERMINE TIMESCALES OF PROCESSES ***** //    
        // ***** START OF : NUMERICAL METHOD BASED ON TIME SCALES ***** //  
    
        //!< Resolve region where the Total Power is zero for a Plasma Grid.
        //!< This typically occurs when plasma parameters are zero in a cell and 
        //!< other models are off (or zero)... 
        //!< Even in No Plasma Region, cooling processes can still occur, but 
        //!< this is specifically for zero power
        
        D_Debug("\n\n***** DTOKSU::Run() :: Begin Global Step *****\n\n");
        if( HeatTime == 10 ){
            D1_Debug("\n\tNo Net Power Region...\n");
            FM.Force(); 
            HM.AddTime(ForceTime);
            CM.Charge(ForceTime);
            TotalTime += ForceTime;
        }else if( MinTimeStep*2.0 > MaxTimeStep){
            D1_Debug("\n\tComparable Timescales, taking time steps through both "
                << "processes at shorter time scale\n");
            CM.Charge(MinTimeStep);
            HM.Heat(MinTimeStep);
            FM.Force(MinTimeStep);
            TotalTime += MinTimeStep;
        }else{ 
            //!< Else, we can take steps through the smaller one til the sum 
            //!< of the steps is the larger.
            D1_Debug("\n\tDifferent Timescales, taking many time steps through "
                << "quicker process at shorter time scale\n");
            unsigned int j(1);
            bool Loop(true);
            for( j =1; (j*MinTimeStep) < MaxTimeStep && Loop; j ++){
                D1_Debug("\n\tIntermediateStep/MaxTimeStep = " << j*MinTimeStep 
                    << "/" << MaxTimeStep << "\n");

                //!< Take the time step in the faster time process
                if( MinTimeStep == HeatTime ){
                    HM.Heat(MinTimeStep);
                    HM.Record_MassLoss(false);
                    if( Sample->is_gas() ){
                        D1_Debug("\n\tSample is gaseous!\n");
                        Loop=false;
                    }
                    //!< This has to go here, Think break; statement
                    CM.Charge(MinTimeStep); 
                    D1_Debug("\n\tHeat Step Taken.\n");
                    //!< Check that time scales haven't changed significantly 
                    //!< whilst looping... 
                }else if( MinTimeStep == ForceTime ){
                    FM.Force(MinTimeStep);
                    //!< This has to go here, Think break; statement
                    CM.Charge(MinTimeStep); 
                    if( FM.new_cell() ){
                        D1_Debug("\n\tWe have stepped into a new cell!\n");
                        Loop=false;
                    }
                    D1_Debug("\n\tForce Step Taken.\n");
                    //!< Check that time scales haven't changed significantly 
                    //!< whilst looping... 
                }else{
                    std::cerr << "\nUnexpected Timescale Behaviour (1)!";
                }


                //!< Check that time scales haven't changed significantly whilst 
                //!< looping...
                //!< NOTE the FM.ProbeTimeStep() command takes a significant 
                //!< amount of time and has been found to be rarely activated.
                //!< This could still be a nice/necessary check in some
                //!< circumstances. However, I can't see a way to make this 
                //!< faster right now
                double HM_ProbeTime = HM.ProbeTimeStep();
                if( ForceTime/FM.ProbeTimeStep() > 2 ){
                    D1_Debug("\n\tForce TimeStep Has Changed Significantly whilst"
                        << " taking small steps...\n");
                    j ++;
                    //!< Can't do this: MaxTimeStep = j*MinTimeStep; 
                    //!< as we change MaxTimeStep...
                    break; 
                }
                if( HeatTime/HM_ProbeTime > 2 ){
                    D1_Debug("\n\tHeat TimeStep Has Changed Significantly whilst"
                        << " taking small steps...\n");

                    j ++;
                    //!< Can't do this: MaxTimeStep = j*MinTimeStep; 
                    //!< as we change MaxTimeStep...
                    break; 
                }

                if( HM_ProbeTime == 1 ) break; //!< Thermal Equilibrium Reached
            }

            // Take a time step in the slower time process
            D1_Debug("\n\t*STEP* = " << (j-1)*MinTimeStep << "\nMaxTimeStep = " 
                << MaxTimeStep << "\nj = " << j << "\n");
            D1_Debug("\n\tHeatTime == " << HeatTime << "\nForceTime == " 
                << ForceTime << "\n");
            
            if( MaxTimeStep == HeatTime ){
                if( j > 1 )
                    HM.Heat((j-1)*MinTimeStep);
                else
                    HM.AddTime(MinTimeStep);
                D_Debug("\nHeat Step Taken."); 
            }else if( MaxTimeStep == ForceTime ){
                if( j > 1 )
                    FM.Force((j-1)*MinTimeStep); 
                else
                    FM.AddTime(MinTimeStep);
                D_Debug("\nForce Step Taken."); 
            }else{  std::cerr << "\nUnexpected Timescale Behaviour! (2)";   }

            TotalTime += (j-1)*MinTimeStep;
            //!< This is effectively 'an extra step' which is necessary because 
            //!< we need the last step to be charging
            //!< So we set the time step to be arbitrarily small... 
            //!< Not the best practice ever but okay.
            CM.Charge(1e-100);
        }
        // ***** END OF : NUMERICAL METHOD BASED ON TIME SCALES ***** //    
        D_Debug("\n\n***** DTOKSU::Run() :: End Global Step *****\n\n");
        D_Debug("\n\tTemperature = " << Sample->get_temperature()
            << "\n\tMinTimeStep = " << MinTimeStep << "\n\tChargeTime = " 
            << ChargeTime << "\n\tForceTime = " << ForceTime
            << "\n\tHeatTime = " << HeatTime << "\n");

        //!< Update the plasma data from the plasma grid for all models...
        D_Debug("\n\n***** DTOKSU::Run() :: Updating Plasma Data *****\n\n");
        cm_InGrid = CM.update_plasmadata();
        hm_InGrid = HM.update_plasmadata();
        fm_InGrid = FM.update_plasmadata();
        
        assert( cm_InGrid == fm_InGrid 
            && fm_InGrid == hm_InGrid 
            && hm_InGrid == cm_InGrid );

        D_Debug("\n\n***** DTOKSU::Run() :: Record Plasma Data and Impurities *****\n\n");
        CM.RecordPlasmadata("pd.txt");
        HM.Record_MassLoss(false);
        //HM.RecordPlasmadata("hm_pd.txt");
        //FM.RecordPlasmadata("fm_pd.txt");
        print();
        Pause();
        // ***** START OF : DETERMINE IF END CONDITION HAS BEEN REACHED ***** //
        D_Debug("\n\n***** DTOKSU::Run() :: Determine If End Conditions Satisfied *****\n\n");
        if( Sample->is_gas() 
            && Sample->get_superboilingtemp() <= Sample->get_temperature() ){
            HM.Record_MassLoss(true);
            std::cout << "\n\nSample has Boiled!";
            break;
        }else if( Sample->is_gas() 
            && Sample->get_superboilingtemp() > Sample->get_temperature() ){
            HM.Record_MassLoss(true);
            std::cout << "\n\nSample has Evaporated!";
            break;
        }else if( Sample->is_split() ){
            std::cout << "\n\nSample has broken up into two large parts";
            //!< Could replace with return 3; and leave off the end check?...
            break; 
        }else if( Sample->is_gas() ){
            HM.Record_MassLoss(true);
            std::cout << "\n\nSample has vapourised";
            break;
        }else if( HeatTime == 1 ){
            std::cout << "\n\nThermal Equilibrium reached!";
            break;
        }else{
            if( CoreBound.Grid_Pos.size() > 2 ){
                if( Boundary_Check(true) ){
                    HM.Record_MassLoss(true);
                    std::cout << "\n\nCollision with Core!";
                    break;
                }
            }
            if( WallBound.Grid_Pos.size() > 2 ){
                if( Boundary_Check(false) ){
                    break;
                }
            }
        }
        // ***** END OF : DETERMINE IF END CONDITION HAS BEEN REACHED ***** //
    }

    D_Debug("\n\n***** DTOKSU::Run() :: Determine Return and Termination Status *****\n\n");
    if( fabs(1 - (HM.get_totaltime()/FM.get_totaltime())) > 0.001  
        || fabs(1 - (FM.get_totaltime()/CM.get_totaltime())) > 0.001 ){
        std::cout << "\nWarning! Total Times recorded by processes don't "
            << "match!\nHM.get_totaltime() = " << HM.get_totaltime() 
            << "\nFM.get_totaltime() = " << FM.get_totaltime() 
            << "\nCM.get_totaltime() = " << CM.get_totaltime() 
            << "\n\nTotalTime = " << TotalTime;
    }
    if( !cm_InGrid ){
        std::cout << "\nSample has left simulation domain";
        return 1;
    }else if( HeatTime == 1 ){
        std::cout << "\nEnd of Run, exiting due to Thermal Equilibrium being "
            << "reached!\n\n";
        return 2;
    }else if( Sample->is_split() && !Sample->is_gas() ){
        return 3; //!< Return status=3 for continue simulating Breakup condition
    }else if( Sample->is_gas() ){
        std::cout << "\nSample has boiled, evaporated or vapourised!\n\n";
        return 4;
    }else if( TotalTime > MaxTime ){
        std::cout << "\nMax Simulation Time Exceeded!\n\n";
        return 5;

    }
    else if( ErrorFlag ){
        std::cout << "\nGeneric run-time error!\n\n";
        return 10;
    }

    std::cout << "\nFinished DTOKS-U run.";
    return 0;
}
