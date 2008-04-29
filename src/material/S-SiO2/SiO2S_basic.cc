/*****************************************************************************/
/*                                                                           */
/*              8888888         88888888         88888888                    */
/*            8                8                8                            */
/*           8                 8                8                            */
/*           8                  88888888         88888888                    */
/*           8      8888                8                8                   */
/*            8       8                 8                8                   */
/*              888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS material database Version 0.4                                        */
/*  Last update: Feb 17, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: SiO2 as semiconductor


#include "PMI.h"

class GSS_SiO2S_BasicParameter : public PMIS_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of SiO2.
  PetscScalar PERMEABI;  // The relative megnetic permeability of SiO2.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.
  void   Basic_Init()
  {
    PERMITTI = 3.900000e+00;
    PERMEABI = 1.0;
    AFFINITY = 0.970000e+00*eV;
    DENSITY  = 2.260000e-03*kg*pow(cm,-3);
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }
  GSS_SiO2S_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_SiO2S_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIS_BasicParameter* PMIS_SiO2S_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_SiO2S_BasicParameter(env);
  }
}
