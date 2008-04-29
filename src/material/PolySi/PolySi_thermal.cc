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

class GSS_PolySi_Thermal : public PMIC_Thermal
{
private:
  PetscScalar A_SP_HEA;	// First parameter for the specific heat model of the material.
  PetscScalar B_SP_HEA;	// Second parameter for the specific heat model of the material.
  PetscScalar C_SP_HEA;	// Third parameter for the specific heat model of the material.
  PetscScalar D_SP_HEA;	// Fourth parameter for the specific heat model of the material.
  PetscScalar F_SP_HEA;	// Fifth parameter for the specific heat model of the material.
  PetscScalar G_SP_HEA;	// Sixth parameter for the specific heat model of the material.
  PetscScalar A_TH_CON;	// First parameter for the thermal conductivity model of the material.
  PetscScalar B_TH_CON;	// Second parameter for the thermal conductivity model of the material.
  PetscScalar C_TH_CON;	// Third parameter for the thermal conductivity model of the material.
  PetscScalar E_TH_CON;	// Fifth parameter for the thermal conductivity model of the material.
  PetscScalar D_TH_CON;	// Fourth parameter for the thermal conductivity model of the material.
  void   Thermal_Init()
  {
    A_SP_HEA =  7.370000e+02*J/kg/K;
    B_SP_HEA =  4.420000e-01*J/kg/pow(K,2);
    C_SP_HEA =  0.000000e+00*J/kg/pow(K,3);
    D_SP_HEA =  0.000000e+00*J/kg*K;
    F_SP_HEA =  0.000000e+00*J/kg/pow(K,4);
    G_SP_HEA =  0.000000e+00*J/kg/pow(K,5);
    A_TH_CON =  0.437000e+00*cm*K/J*s;
    B_TH_CON =  0.000000e+00*cm/J*s;
    C_TH_CON =  0.000000e+00*cm/J*s/K;
    E_TH_CON =  0.000000e+00;
    D_TH_CON =  0.000000e+00*cm/J*s*pow(K,1-E_TH_CON);
  }
public:
  PetscScalar HeatCapacity  (const PetscScalar &Tl) const 
  {
    return A_SP_HEA + B_SP_HEA*Tl + C_SP_HEA*Tl*Tl + D_SP_HEA/Tl/Tl
           + F_SP_HEA*Tl*Tl*Tl + G_SP_HEA*Tl*Tl*Tl*Tl;
  }
  PetscScalar HeatConduction(const PetscScalar &Tl) const 
  {
    return 1.0/(A_TH_CON + B_TH_CON*Tl + C_TH_CON*Tl*Tl + D_TH_CON*pow(Tl,E_TH_CON));
  }
  GSS_PolySi_Thermal(const PMIC_Environment &env):PMIC_Thermal(env)
  {
    Thermal_Init();
  }
  ~GSS_PolySi_Thermal()
  {
  }
}
;

extern "C"
{
  PMIC_Thermal* PMIC_PolySi_Thermal_Default (const PMIC_Environment& env)
  {
    return new GSS_PolySi_Thermal(env);
  }
}
