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
// Material Type: Air


#include "PMI.h"

class GSS_Air_BasicParameter : public PMII_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of Air.
  PetscScalar PERMEABI;  // The relative megnetic permeability of Air.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.

  void   Basic_Init()
  {
    PERMITTI = 1.0;
    PERMEABI = 1.0;
    AFFINITY = 0.0*eV;
    DENSITY  = 1.290000e-06*kg*pow(cm,-3);
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }
  GSS_Air_BasicParameter(const PMII_Environment &env):PMII_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_Air_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMII_BasicParameter* PMII_Air_BasicParameter_Default (const PMII_Environment& env)
  {
    return new GSS_Air_BasicParameter(env);
  }
}
