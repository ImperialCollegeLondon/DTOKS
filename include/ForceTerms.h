/** @file ForceTerms.h
 *  @brief Contains functions defining the physics of different
 *  force terms
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __FORCETERMS_H_INCLUDED__
#define __FORCETERMS_H_INCLUDED__

#include "Term.h"
#include "PlasmaFluxes.h"
#include "MathHeader.h"

namespace Term{

/** @name Force Terms
 *  @brief functions for different terms in force equation
 *
 *  These functions all return the acceleration due to this force term
 *  in the equation of motion for a dust grain in a tokamak plasma.
 *  The number of functions is equal to FMN+1
 *  @return threevector acceleration of dust due to force term
 */
///@{
/** @brief The acceleration due to gravity
 *  @return The acceleration in m/s^2 due to gravity
 */
struct Gravity:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "Gravity"; };
};
/** @brief Lorentz force due to electric and magnetic fields
 *  @return The acceleration in m/s^2 due to the lorentz force
 */
struct LorentzForce:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "LorentzForce"; };
};
/** @brief SOML ion drag model due to collisions of dust with ions
 *  @return The acceleration in m/s^2 due to SOML Ion Drag force
 */
struct SOMLIonDrag:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "SOMLIonDrag"; };
};
/** @brief SMOML ion drag model due to collisions of dust with ions
 *  @return The acceleration in m/s^2 due to SMOML Ion Drag force
 */
struct SMOMLIonDrag:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "SMOMLIonDrag"; };
};
/** @brief DTOKS ion drag model due to collisions of dust with ions
 *  @return The acceleration in m/s^2 due to DTOKS Ion Drag force
*/
struct DTOKSIonDrag:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "DTOKSIonDrag"; };
};
/** @brief DUSTT ion drag model due to collisions of dust with ions
 *  @return The acceleration in m/s^2 due to DUSTT Ion Drag force
 */
struct DUSTTIonDrag:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "DUSTTIonDrag"; };
};
/** @brief Hybrid ion drag model due to collisions of dust with ions
 *  @return The acceleration in m/s^2 due to Hybrid Ion Drag force
 */
struct HybridIonDrag:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "HybridIonDrag"; };
};
/** @brief Lloyd ion drag model due to collisions of dust with ions
 *  @return The acceleration in m/s^2 due to Lloyd Ion Drag force
 */
struct LloydIonDrag:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "LloydIonDrag"; };
};

/** @brief Neutral drag due to collisions of dust with neutrals
 *  @return The acceleration in m/s^2 due to Neutral Drag force
 */
struct NeutralDrag:ForceTerm{
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "NeutralDrag"; };
};
/** @brief Neutral drag due to collisions of dust with neutrals
 *  @return The acceleration in m/s^2 due to the Rocket force
 */
struct RocketForce:ForceTerm{
    double OldTemp;
    threevector Evaluate(const Matter* Sample, 
        std::shared_ptr<PlasmaData> Pdata, threevector velocity);
    std::string PrintName(){ return "RocketForce"; };
};
///@}
}

#endif /* __FORCETERMS_H_INCLUDED__ */
