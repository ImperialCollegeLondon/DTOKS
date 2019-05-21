/** @file PlasmaFluxes.h
 *  @brief Contains functions declaring the default plasma fluxes
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __PLASMAFLUXES_H_INCLUDED__
#define __PLASMAFLUXES_H_INCLUDED__

//#define PLASMAFLUX_DEBUG

#include "Matter.h"
#include <memory>
#include "PlasmaData.h"

namespace Flux{
/** @name Plasma Fluxes
 *  @brief Functions which calculate the plasma flux to the sphere
 *
 *  These functions calculate the approximate flux to a charged sphere
 *  under different regimes relevant to different charging models. These
 *  are used to provide consistency in calculations made by different 
 *  models
 */
///@{
/** @brief Calculate the Ion Flux as specified by OML
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux according to OML theory
 */
double OMLIonFlux(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double Potential);
/** @brief Calculate the Ion Flux as specified by MOML
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux following MOML theory for negative dust
 */
double MOMLIonFlux(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double Potential);
/** @brief Calculate the Ion Flux as specified by SOML
 *  For the case of no flow, to avoid dividing by zero we return OMLIonFlux.
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by SOML
 */
double SOMLIonFlux(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double Potential);
/** @brief Calculate the Ion Flux as specified by SMOML
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by SMOML
 */
double SMOMLIonFlux(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double Potential);
/** @brief Calculate the Ion Flux as specified by PHL
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by PHL
 */
double PHLElectronFlux(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double Potential);
/** @brief Calculate the Ion Flux as specified originally in DTOKS
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by DTOKS
 */
double DTOKSIonFlux(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double Potential);
/** @brief Calculate the electron flux as specified originally in DTOKS
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the electron flux as originally given by DTOKS
 */
double DTOKSElectronFlux(const std::shared_ptr<PlasmaData> Pdata, 
    const double Potential);
/** @brief Calculate the electron flux as specified by OML
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by OML
 */
double OMLElectronFlux(const std::shared_ptr<PlasmaData> Pdata, 
    const double Potential);
/** @brief Calculate the thermal neutral flux 
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @return the  thermal neutral flux
 */
double NeutralFlux(const std::shared_ptr<PlasmaData> Pdata);

/** @brief flux of particles away from liquid surface due to evaporation
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param DustTemperature The temperature of the dust
 *  @return The flux of evaporating particles away from the liquid
 */
double EvaporationFlux(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double DustTemperature);
/** @brief Thermionic electron emission yield
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @return the thermionic electron emission yield coefficient
 */
double DeltaTherm(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata);
/** @brief Thermionic electron emission flux
 *  @param Sample Const pointer to class containing all data about matter
 *  @return the thermionic electron emission flux
 */
double ThermFlux(const Matter* Sample);
/** @brief Thermionic electron emission flux with schottky correction
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @param Potential the potential of the dust grain surface
 *  @return the thermionic electron emission flux
 */
double ThermFluxSchottky(const Matter* Sample, 
    const std::shared_ptr<PlasmaData> Pdata, const double Potential);
/** @brief Secondary electron emission yield
 *  @param Sample Const pointer to class containing all data about matter
 *  @param Pdata Const pointer to data structure with information about plasma
 *  @return the secondary electron emission yield coefficient
 */
double DeltaSec(const Matter* Sample, const std::shared_ptr<PlasmaData> Pdata);
///@}

}
#endif /* __PLASMAFLUXES_H_INCLUDED__ */
