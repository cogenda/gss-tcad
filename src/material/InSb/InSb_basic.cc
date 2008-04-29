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
// Material Type: InSb


#include "PMI.h"

class GSS_InSb_BasicParameter : public PMIS_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of InAs.
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar PERMEABI;  // The relative megnetic permeability of InAs.
  PetscScalar DENSITY;   // Specific mass density for the material.

  void   Basic_Init()
  {
    //Source: Semiconductors on NSM
    PERMITTI =  16.8;
    AFFINITY =  4.59*eV;
    PERMEABI =  1.0;
    DENSITY  =  5.77E-03*kg*pow(cm,-3);
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const { return DENSITY;  }
  PetscScalar Permittivity  ()                      const { return PERMITTI; }
  PetscScalar Permeability  ()                      const { return PERMEABI; }
  PetscScalar Affinity      (const PetscScalar &Tl) const { return AFFINITY; }

  GSS_InSb_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_InSb_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIS_BasicParameter* PMIS_InSb_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_InSb_BasicParameter(env);
  }
}
