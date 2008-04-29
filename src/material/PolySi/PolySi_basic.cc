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
// Material Type: Poly Silicon


#include "PMI.h"

class GSS_PolySi_BasicParameter : public PMIC_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of silicon.
  PetscScalar PERMEABI;  // The relative megnetic permeability of silicon.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.
  void   Basic_Init()
  {
    PERMITTI = 1.1800000e+01;
    PERMEABI = 1.0;//---not magnetism material,assigned its permeability equals 1.
    AFFINITY = 4.170000e+00*eV;
    DENSITY  = 2.320000e-03*kg*pow(cm,-3);
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }
  GSS_PolySi_BasicParameter(const PMIC_Environment &env):PMIC_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_PolySi_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIC_BasicParameter* PMIC_PolySi_BasicParameter_Default (const PMIC_Environment& env)
  {
    return new GSS_PolySi_BasicParameter(env);
  }
}
