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
/*  Last update: July 25, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Nitride


#include "PMI.h"

class GSS_Nitride_BasicParameter : public PMII_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of Nitride.
  PetscScalar PERMEABI;  // The relative megnetic permeability of Nitride.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.
  void   Basic_Init()
  {
    PERMITTI = 7.500000E+00;
    PERMEABI = 1.0;
    AFFINITY = 9.700000E-01*eV;
    DENSITY  = 3.440000E-03*kg*pow(cm,-3);
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }
  GSS_Nitride_BasicParameter(const PMII_Environment &env):PMII_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_Nitride_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMII_BasicParameter* PMII_Nitride_BasicParameter_Default (const PMII_Environment& env)
  {
    return new GSS_Nitride_BasicParameter(env);
  }
}
