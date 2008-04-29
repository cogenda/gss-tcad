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
/*  Last update: May 27, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: TiSi2


#include "PMI.h"

class GSS_TiSi2_BasicParameter : public PMIC_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity.
  PetscScalar PERMEABI;  // The relative megnetic permeability.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.
  void   Basic_Init()
  {
    PERMITTI = 1.000000e+00;//TiSi2
    PERMEABI = 1.0;
    AFFINITY = 0.970000e+00*eV;//use value for Al, should be corrected later.
    DENSITY  = 4.040000e-03*kg*pow(cm,-3);//TiSi2
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }
  GSS_TiSi2_BasicParameter(const PMIC_Environment &env):PMIC_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_TiSi2_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIC_BasicParameter* PMIC_TiSi2_BasicParameter_Default (const PMIC_Environment& env)
  {
    return new GSS_TiSi2_BasicParameter(env);
  }
}
