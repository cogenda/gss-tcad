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
// Material Type: Si(1-x)Ge(x).


#include "PMI.h"

class GSS_SiGe_BasicParameter : public PMIS_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of Si(1-x)Ge(x).
  PetscScalar EPS_X1;
  PetscScalar EPS_X2;
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar AF_X1;
  PetscScalar AF_X2;
  PetscScalar PERMEABI;  // The relative megnetic permeability of silicon.
  PetscScalar DENSITY;   // Specific mass density for the material.

  void   Basic_Init()
  {
    PERMITTI = 1.180000e+01;
    EPS_X1   = 0.000000E+00;
    EPS_X2   = 0.000000E+00;
    AFFINITY = 4.170000e+00*eV;
    AF_X1    = 0.000000E+00;
    AF_X2    = 0.000000E+00;
    PERMEABI = 1.0;
    DENSITY  = 2.320000e-03*kg*pow(cm,-3);
  }
public:
  PetscScalar Density(const PetscScalar &Tl) const 
  { 
  	return DENSITY;  
  }

  PetscScalar Permittivity() const 
  {
        PetscScalar mole_x = ReadxMoleFraction();
        return PERMITTI + EPS_X1*mole_x + EPS_X2*mole_x*mole_x;
  }

  PetscScalar Permeability() const          
  { 
  	return PERMEABI; 
  }

  PetscScalar Affinity(const PetscScalar &Tl) const 
  {
        PetscScalar mole_x = ReadxMoleFraction();
        return AFFINITY + AF_X1*mole_x + AF_X2*mole_x*mole_x;
  }

public:
  GSS_SiGe_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_SiGe_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIS_BasicParameter* PMIS_SiGe_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_SiGe_BasicParameter(env);
  }
}
