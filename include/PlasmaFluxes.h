/** @file PlasmaFluxes.h
 *  @brief Contains functions declaring the default plasma fluxes
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 */

#ifndef __PLASMAFLUXES_H_INCLUDED__
#define __PLASMAFLUXES_H_INCLUDED__

#include "Matter.h"
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
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by DTOKS
 */
double OMLIonFlux(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Calculate the Ion Flux as specified by SOML
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by SOML
 */
double SOMLIonFlux(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Calculate the Ion Flux as specified by SMOML
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by SMOML
 */
double SMOMLIonFlux(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Calculate the Ion Flux as specified by PHL
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by PHL
 */
double PHLElectronFlux(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Calculate the Ion Flux as specified originally in DTOKS
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by DTOKS
 */
double DTOKSIonFlux(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Calculate the electron flux as specified originally in DTOKS
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the electron flux as originally given by DTOKS
 */
double DTOKSElectronFlux(std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Calculate the electron flux as specified by OML
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the normalised potential on the dust grain
 *  @return the ion flux as originally given by OML
 */
double OMLElectronFlux(std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Calculate the thermal neutral flux 
 *  @param Pdata Pointer to data structure with information about plasma
 *  @return the  thermal neutral flux
 */
double NeutralFlux(std::shared_ptr<PlasmaData> Pdata);

/** @brief flux of particles away from liquid surface due to evaporation
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 */
double EvaporationFlux(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double DustTemperature);
/** @brief Thermionic electron emission yield
 *  @return the thermionic electron emission yield coefficient
 */
double DeltaTherm();
/** @brief Thermionic electron emission flux
 *  @param Sample Pointer to class containing all data about matter
 *  @return the thermionic electron emission flux
 */
double ThermFlux(Matter* Sample);
/** @brief Thermionic electron emission flux with schottky correction
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 *  @param Potential the potential of the dust grain surface
 *  @return the thermionic electron emission flux
 */
double ThermFluxSchottky(Matter* Sample, std::shared_ptr<PlasmaData> Pdata, double Potential);
/** @brief Secondary electron emission yield
 *  @param Sample Pointer to class containing all data about matter
 *  @param Pdata Pointer to data structure with information about plasma
 *  @return the secondary electron emission yield coefficient
 */
double DeltaSec(Matter* Sample, std::shared_ptr<PlasmaData> Pdata);
///@}

}
#endif /* __PLASMAFLUXES_H_INCLUDED__ */