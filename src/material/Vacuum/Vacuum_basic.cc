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
// Material Type: Vacuum


#include "PMI.h"

class GSS_Vacuum_BasicParameter : public PMIV_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of silicon.
  PetscScalar PERMEABI;  // The relative megnetic permeability of silicon.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar DENSITY;   // Specific mass density for the material.
  void   Basic_Init()
  {
    PERMITTI = 1.0;//corrected by zhangxih
    PERMEABI = 1.0;
    AFFINITY = 0.0*eV;
    DENSITY  = 0.0*kg*pow(cm,-3);//corrected by zhangxih
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }
  GSS_Vacuum_BasicParameter(const PMIV_Environment &env):PMIV_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_Vacuum_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIV_BasicParameter* PMIV_Vacuum_BasicParameter_Default (const PMIV_Environment& env)
  {
    return new GSS_Vacuum_BasicParameter(env);
  }
}
