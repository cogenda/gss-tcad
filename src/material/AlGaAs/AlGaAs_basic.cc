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
// Material Type: AlGaAs


#include "PMI.h"

class GSS_AlGaAs_BasicParameter : public PMIS_BasicParameter
{
private:
  PetscScalar PERMITTI;  // The relative dielectric permittivity of AlGaAs.
  PetscScalar EPS_X1;
  PetscScalar EPS_X2;
  PetscScalar AFFINITY;  // The electron affinity for the material.
  PetscScalar AF_X0;
  PetscScalar AF_X1;
  PetscScalar AF_X2;
  PetscScalar AF_X3;
  PetscScalar AF_X4;
  PetscScalar AF_X5;
  PetscScalar AF_XL;
  PetscScalar PERMEABI;  // The relative megnetic permeability of AlGaAs.
  PetscScalar DENSITY;   // Specific mass density for the material.

  void   Basic_Init()
  {
    PERMITTI =  1.310000E+01;
    EPS_X1   =  0.000000E+00;
    EPS_X2   =  0.000000E+00;
    AFFINITY =  4.070000e+00*eV;
    AF_X0    =  0.000000E+00;
    AF_X1    = -1.100000E+00;
    AF_X2    =  0.000000E+00;
    AF_X3    = -4.300000E-01;
    AF_X4    = -1.400000E-01;
    AF_X5    =  0.000000E+00;
    AF_XL    =  4.500000E-01;
    PERMEABI =  1.0;
    DENSITY  =  5.317600E-03*kg*pow(cm,-3);
  }
public:
  PetscScalar Density       (const PetscScalar &Tl) const 
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
  PetscScalar Affinity      (const PetscScalar &Tl) const 
  {
        PetscScalar mole_x = ReadxMoleFraction();
        if(mole_x<AF_XL)
                return AFFINITY + AF_X0 + AF_X1*mole_x + AF_X2*mole_x*mole_x;
        else
                return AFFINITY + AF_X3 + AF_X4*mole_x + AF_X5*mole_x*mole_x;
  }

  GSS_AlGaAs_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
    Basic_Init();
  }
  ~GSS_AlGaAs_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIS_BasicParameter* PMIS_AlGaAs_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_AlGaAs_BasicParameter(env);
  }
}
