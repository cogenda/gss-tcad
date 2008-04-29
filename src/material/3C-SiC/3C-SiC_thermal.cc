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
// Material Type: 3C-SiC


#include "PMI.h"

class GSS_SiC3C_Thermal : public PMIS_Thermal
{
private:
  PetscScalar A_SP_HEA;	// First parameter for the specific heat model of the material.
  PetscScalar B_SP_HEA;	// Second parameter for the specific heat model of the material.
  PetscScalar C_SP_HEA;	// Third parameter for the specific heat model of the material.
  PetscScalar D_SP_HEA;	// Fourth parameter for the specific heat model of the material.
  PetscScalar F_SP_HEA;	// Fifth parameter for the specific heat model of the material.
  PetscScalar G_SP_HEA;	// Sixth parameter for the specific heat model of the material.	 
  PetscScalar T300;
  void   Thermal_Init()
  {
     A_SP_HEA =  8.509000e+02*J/kg/K;
     B_SP_HEA =  1.522000e-01*J/kg/pow(K,2);
     C_SP_HEA =  0.000000e+00*J/kg/pow(K,3);
     D_SP_HEA = -1.582000e+07*J/kg*K;
     F_SP_HEA =  0.000000e+00*J/kg/pow(K,4);
     G_SP_HEA =  0.000000e+00*J/kg/pow(K,5);
     T300     =  300.0*K;
  }
public:
  //---------------------------------------------------------------------------
  // Heat Capacity, no exact values, use paramters for Si.
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const
  {
    return A_SP_HEA + B_SP_HEA*Tl + C_SP_HEA*Tl*Tl + D_SP_HEA/Tl/Tl
           + F_SP_HEA*Tl*Tl*Tl + G_SP_HEA*Tl*Tl*Tl*Tl;
  }
  AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const
  {
    return A_SP_HEA + B_SP_HEA*Tl + C_SP_HEA*Tl*Tl + D_SP_HEA/Tl/Tl
           + F_SP_HEA*Tl*Tl*Tl + G_SP_HEA*Tl*Tl*Tl*Tl;
  }

  //---------------------------------------------------------------------------
  // Heat Conduction, source: Semiconductors on NSM
  PetscScalar HeatConduction(const PetscScalar &Tl) const
  {
    return 3.6*W/cm/K;
  }
  AutoDScalar HeatConduction(const AutoDScalar &Tl) const
  {
    return 3.6*W/cm/K;
  }

// constructor and destructor  
public:     
  GSS_SiC3C_Thermal(const PMIS_Environment &env):PMIS_Thermal(env)
  {
    Thermal_Init();
  }
  ~GSS_SiC3C_Thermal()
  {
  }
}
;

extern "C"
{
  PMIS_Thermal* PMIS_SiC3C_Thermal_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC3C_Thermal(env);
  }
}
