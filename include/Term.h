/** @file Term.h
 *  @brief Contains a class defining the behaviour of a physics term
 *  
 *  Structure which defines the behaviour of physical terms acting 
 *  within a model that actually describe the source terms of different 
 *  physics.
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __TERM_H_INCLUDED__
#define __TERM_H_INCLUDED__

#include "Matter.h"
#include "PlasmaData.h"

#include <memory>
#include <string>

/** @struct ForceTerm
 *  @brief struct defining the behaviour of Force terms within a dynamics model
 *  
 *  Abstract structure which defines the behaviour of force terms acting 
 *  within the force model that describe the source terms of different 
 *  accelerations.
 */
struct ForceTerm{
    virtual threevector Evaluate(Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata)=0;
    virtual std::string PrintName()=0;
};

/** @struct HeatTerm
 *  @brief struct defining the behaviour of heat terms within a heating model
 *  
 *  Abstract structure which defines the behaviour of heat terms acting 
 *  within the heating model that describe the source terms of different 
 *  power fluxes.
 */
struct HeatTerm{
    virtual double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata,
     double Temp)=0;
    virtual std::string PrintName()=0;
};

/** @struct ChargingTerm
 *  @brief struct defining the behaviour of different charging models
 *  
 *  Abstract structure which defines the behaviour of the charging terms acting 
 *  within the charging model.
 */
struct CurrentTerm{
    virtual void Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double& Potential)=0;
    virtual std::string PrintName()=0;
};

#endif /* __TERM_H_INCLUDED__ */