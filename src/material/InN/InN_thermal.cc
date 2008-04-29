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
// Material Type: InN


#include "PMI.h"

class GSS_InN_Thermal : public PMIS_Thermal
{
public:
  //---------------------------------------------------------------------------
  // Heat Capacity
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    PetscScalar Cv = 2.2*J/(K*pow(cm,3));  //Source: Semiconductors on NSM
    return Cv;
  }
  AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const
  {
    PetscScalar Cv = 2.2*J/(K*pow(cm,3));  //Source: Semiconductors on NSM
    return Cv;
  }
  
  //---------------------------------------------------------------------------
  // Heat Conduction
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 0.45*W/cm/K;  //Source: Semiconductors on NSM
  }
  AutoDScalar HeatConduction(const AutoDScalar &Tl) const
  {
    return 0.45*W/cm/K;  //Source: Semiconductors on NSM
  }


// constructor and destructor  
public:     
  GSS_InN_Thermal(const PMIS_Environment &env):PMIS_Thermal(env)
  {
    
  }
  ~GSS_InN_Thermal()
  {
  }
}
;

extern "C"
{
  PMIS_Thermal* PMIS_InN_Thermal_Default (const PMIS_Environment& env)
  {
    return new GSS_InN_Thermal(env);
  }
}
