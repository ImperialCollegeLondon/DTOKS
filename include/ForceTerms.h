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
 */
struct Gravity:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "Gravity"; };
};
/** @brief Centrifugal force due to cylindrical coordinate system
  */
struct Centrifugal:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "Centrifugal"; };
};
/** @brief Lorentz force due to electric and magnetic fields
  */
struct LorentzForce:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "LorentzForce"; };
};
/** @brief SOML ion drag model due to collisions of dust with ions
  */
struct SOMLIonDrag:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "SOMLIonDrag"; };
};
/** @brief SMOML ion drag model due to collisions of dust with ions
  */
struct SMOMLIonDrag:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "SMOMLIonDrag"; };
};
/** @brief DTOKS ion drag model due to collisions of dust with ions
  */
struct DTOKSIonDrag:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "DTOKSIonDrag"; };
};
/** @brief DUSTT ion drag model due to collisions of dust with ions
  */
struct DUSTTIonDrag:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "DUSTTIonDrag"; };
};
/** @brief Hybrid ion drag model due to collisions of dust with ions
  */
struct HybridIonDrag:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "HybridIonDrag"; };
};
/** @brief Neutral drag due to collisions of dust with neutrals
  */
struct NeutralDrag:ForceTerm{
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "NeutralDrag"; };
};

struct RocketForce:ForceTerm{
    double OldTemp;
    threevector Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
    std::string PrintName(){ return "RocketForce"; };
};
///@}
}

#endif /* __FORCETERMS_H_INCLUDED__ */
