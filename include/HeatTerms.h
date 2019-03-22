/** @file HeatTerms.h
 *  @brief Contains functions defining the physics of different
 *  heat terms
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __HEATTERMS_H_INCLUDED__
#define __HEATTERMS_H_INCLUDED__

#include "Term.h"
#include "PlasmaFluxes.h"

namespace Term{

/** @name Heating Terms
 *  @brief functions for different terms in heating equation
 *
 *  These functions all return the power due to a particular heat term
 *  in the equation of heating for a dust grain in a tokamak plasma.
 *  The number of functions is equal to HMN+1
 *  @return double the power to the dust grain in kW
 *  @param DustTemperature the temperature of the dust, needed for RK4
 */
///@{
/** @brief Power to surface due to black body radiation
 */
struct EmissivityModel:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "EmissivityModel"; };
};
/** @brief Power to surface due evaporative loss of particles
 */
struct EvaporationModel:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "EvaporationModel"; };
};
/** @brief Power to surface due to contact with air
 */
struct NewtonCooling:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "NewtonCooling"; };
};
/** @brief Power to surface due to neutral bombardment
 */
struct NeutralHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "NeutralHeatFlux"; };
};
/** @brief Power to surface due to SOML ion bombardment 
 */
struct SOMLIonHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "SOMLIonHeatFlux"; };
};
/** @brief Power to surface due to SOML neutral recombination
 */
struct SOMLNeutralRecombination:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "SOMLNeutralRecombination"; };
};
/** @brief Power to surface due to SMOML ion bombardment 
 */
struct SMOMLIonHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "SMOMLIonHeatFlux"; };
};
/** @brief Power to surface due to SMOML neutral recombination
 */
struct SMOMLNeutralRecombination:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "SMOMLNeutralRecombination"; };
};
/** @brief Power to surface due to secondary electron emission
 */
struct SEE:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "SEE"; };
};
/** @brief Power to surface due to thermionic electron emission
 */
struct TEE:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "TEE"; };
};
/** @brief Power to surface due to PHL electron bombardment 
 */
struct PHLElectronHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "PHLElectronHeatFlux"; };
};
/** @brief Power to surface due to OML electron bombardment 
 */
struct OMLElectronHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "OMLElectronHeatFlux"; };
};
/** @brief Power to surface due to DTOKS secondary electron emission
 */
struct DTOKSSEE:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "DTOKSSEE"; };
};
/** @brief Power to surface due to DTOKS secondary electron emission
 */
struct DTOKSTEE:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "DTOKSTEE"; };
};
/** @brief Power to surface due to DTOKS ion bombardment 
 */
struct DTOKSIonHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "DTOKSIonHeatFlux"; };
};
/** @brief Power to surface due to DTOKS neutral recombination
 */
struct DTOKSNeutralRecombination:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "DTOKSNeutralRecombination"; };
};
/** @brief Power to surface due to DTOKS electron bombardment 
 */
struct DTOKSElectronHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "DTOKSElectronHeatFlux"; };
};
/** @brief Power to surface due to DUSTT ion bombardment 
 */
struct DUSTTIonHeatFlux:HeatTerm{
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
    std::string PrintName(){ return "DUSTTIonHeatFlux"; };
};
///@}
}

#endif /* __HEATTERMS_H_INCLUDED__ */