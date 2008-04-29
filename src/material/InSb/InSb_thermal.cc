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

class GSS_InSb_Thermal : public PMIS_Thermal
{
public:
  //---------------------------------------------------------------------------
  // Heat Capacity
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    return 0.2*J/g/K;    //Source: Semiconductors on NSM
  }
  AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const
  {
    return 0.2*J/g/K;    //Source: Semiconductors on NSM
  }

  
  //---------------------------------------------------------------------------
  // Heat Conduction
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 0.18*W/cm/K;  //Source: Semiconductors on NSM
  }
  AutoDScalar HeatConduction(const AutoDScalar &Tl) const
  {
    return 0.18*W/cm/K;  //Source: Semiconductors on NSM
  }


// constructor and destructor  
public:     
  GSS_InSb_Thermal(const PMIS_Environment &env):PMIS_Thermal(env)
  {
    
  }
  ~GSS_InSb_Thermal()
  {
  }
}
;

extern "C"
{
  PMIS_Thermal* PMIS_InSb_Thermal_Default (const PMIS_Environment& env)
  {
    return new GSS_InSb_Thermal(env);
  }
}
