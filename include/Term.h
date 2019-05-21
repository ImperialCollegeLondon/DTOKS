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
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Temp the temperature of the dust at which it is evaluated
 *  @return The acceleration in m/s^2 of \p Sample due to the ForceTerm
 */
struct ForceTerm{
    virtual threevector Evaluate(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, 
        const threevector velocity)=0;
    virtual std::string PrintName()=0;
};

/** @struct HeatTerm
 *  @brief struct defining the behaviour of heat terms within a heating model
 *   
 *  Abstract structure which defines the behaviour of heat terms acting 
 *  within the heating model that describe the source terms of different 
 *  power fluxes.
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Temp the temperature of the dust at which it is evaluated
 *  @return The power in kW to the surface of \p Sample due to the HeatTerm
 */
struct HeatTerm{
    virtual double Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, const double Temp)=0;
    virtual std::string PrintName()=0;
};

/** @struct CurrentTerm
 *  @brief struct defining the behaviour of different charging models
 *  
 *  Abstract structure which defines the behaviour of the charging terms acting 
 *  within the charging model.
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Temp the temperature of the dust at which it is evaluated
 *  @return The current flux to the surface of \p Sample due to the CurrentTerm
 */
struct CurrentTerm{
    virtual double Evaluate(const Matter* Sample, 
        const std::shared_ptr<PlasmaData> Pdata, const double Potential)=0;
    virtual std::string PrintName()=0;
};

#endif /* __TERM_H_INCLUDED__ */