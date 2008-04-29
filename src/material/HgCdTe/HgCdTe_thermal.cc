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
// Material Type: HgCdTe


#include "PMI.h"

//Source: Numerical Analys is of the Temperature Field in HgCdTe Detector by Laser Irradiation

class GSS_HgCdTe_Thermal : public PMIS_Thermal
{
public:
  //---------------------------------------------------------------------------
  // Heat Capacity
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    return 150*J/kg/K;
  }
  AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const
  {
    return 150*J/kg/K;
  }

  
  //---------------------------------------------------------------------------
  // Heat Conduction
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 20*W/m/K;
  }
  AutoDScalar HeatConduction(const AutoDScalar &Tl) const
  {
    return 20*W/m/K;
  }

// constructor and destructor  
public:       
  GSS_HgCdTe_Thermal(const PMIS_Environment &env):PMIS_Thermal(env)
  {
  }
  ~GSS_HgCdTe_Thermal()
  {
  }
}
;

extern "C"
{
  PMIS_Thermal* PMIS_HgCdTe_Thermal_Default (const PMIS_Environment& env)
  {
    return new GSS_HgCdTe_Thermal(env);
  }
}
