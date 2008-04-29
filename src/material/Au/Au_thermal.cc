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
// Material Type: Au


#include "PMI.h"

class GSS_Au_Thermal : public PMIC_Thermal
{
public:
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    return 2.49*J/(K*pow(cm,3));
  }
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 3.1*W/(K*cm);
  }
  GSS_Au_Thermal(const PMIC_Environment &env):PMIC_Thermal(env)
  {
    
  }
  ~GSS_Au_Thermal()
  {
  }
}
;

extern "C"
{
  PMIC_Thermal* PMIC_Au_Thermal_Default (const PMIC_Environment& env)
  {
    return new GSS_Au_Thermal(env);
  }
}
