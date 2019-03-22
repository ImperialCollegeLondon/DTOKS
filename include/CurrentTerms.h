/** @file ChargingTerms.h
 *  @brief Contains functions defining the different currents to 
 *  a spherical dust grain
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __CHARGINGTERMS_H_INCLUDED__
#define __CHARGINGTERMS_H_INCLUDED__

#include "Term.h"
#include "PlasmaFluxes.h"
#include "solveMOMLEM.h"

namespace Term{

/** @name Current Terms
 *  @brief different terms in the current balance for a dust grain
 *
 *  These functions all return different components of the currents
 *  to a dust grain. The normalised potential is then calculated by
 *  root-find the current balance to the grain surface assuming 
 *  no net current in steady state.
 *  @return double the normalised surface potential of dust grain
 *  @param a the total electron emission yield
 *  @param guess an initial value guess for root finding algorithm
 */
///@{

struct OMLe:CurrentTerm{
    /** @brief Calculate the OML electron Flux
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the electron flux following OML theory
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "OMLe"; };
};

struct PHLe:CurrentTerm{
    /** @brief Calculate the electron Flux following PHL's theory
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the electron Flux following PHL's theory
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "PHLe"; };
};

struct OMLi:CurrentTerm{
    /** @brief Calculate the OML ion Flux
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the ion flux following OML theory
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "OMLi"; };
};

struct SOMLi:CurrentTerm{
    /** @brief Calculate the SOML ion Flux
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the ion flux following SOML theory
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "SOMLi"; };
};

struct SMOMLi:CurrentTerm{
    /** @brief Calculate the SMOML ion Flux
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the ion flux following SMOML theory
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "SMOMLi"; };
};

struct TEE:CurrentTerm{
    /** @brief Calculate the thermionic electron emission Flux
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the thermionic electron emission Flux
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "TEE"; };
};

struct TEESchottky:CurrentTerm{
    /** @brief Calculate the thermionic electron emission Flux with Schottky
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the thermionic electron emission Flux with Schottky correction
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "TEESchottky"; };
};

struct THSe:CurrentTerm{
    /** @brief Calculate the electron current in magnetic field with THS
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the electron current in magnetic field with THS
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "THSe"; };
};

struct THSi:CurrentTerm{
    /** @brief Calculate the ion current in magnetic field with THS
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the ion current in magnetic field with THS
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "THSi"; };
};

struct DTOKSi:CurrentTerm{
    /** @brief Calculate the ion current following original DTOKS method
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the ion current following original DTOKS method
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "DTOKSi"; };
};

struct DTOKSe:CurrentTerm{
    /** @brief Calculate the electron current following original DTOKS method
     *  @param Sample Pointer to class containing all data about matter
     *  @param Pdata Pointer to data structure with information about plasma
     *  @param Potential the normalised potential on the dust grain
     *  @return the electron current following original DTOKS method
     */
    double Evaluate(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, 
        double Potential);
    std::string PrintName(){ return "DTOKSe"; };
};
///@}
}
#endif /* __CHARGINGTERMS_H_INCLUDED__ */