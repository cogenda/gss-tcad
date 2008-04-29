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
// Material Type: Hg(1-x)Cd(x)Te


#include "PMI.h"

//Source: Two-Dimensional Analysis of Double-Layer Heterojunction HgCdTe Photodiodes

class GSS_HgCdTe_BasicParameter : public PMIS_BasicParameter
{
public:
  PetscScalar Density(const PetscScalar &Tl) const 
  { 
    return 7.6*g*pow(cm,-3);  
  }
  PetscScalar Permittivity() const 
  {
        PetscScalar mole_x = ReadxMoleFraction(0.17,0.443);
        return 20.5 - 15.5*mole_x + 5.7*mole_x*mole_x;
  }
  PetscScalar Permeability() const          
  { 
  	return 1.0; 
  }
  PetscScalar Affinity(const PetscScalar &Tl) const 
  {
        PetscScalar mole_x = ReadxMoleFraction(0.17,0.443);
        PetscScalar Eg = - 0.302*eV + 1.93*mole_x*eV - 0.810*mole_x*mole_x*eV;
                         + 0.832*pow(mole_x,3)*eV + 5.354e-4*(1-2*mole_x)*eV/K*Tl;
        PetscScalar chi = 4.23*eV - 0.813*(Eg - 0.0083*eV);
        return chi;
  }

public:  
  GSS_HgCdTe_BasicParameter(const PMIS_Environment &env):PMIS_BasicParameter(env)
  {
  }
  ~GSS_HgCdTe_BasicParameter()
  {
  }
}
;

extern "C"
{
  PMIS_BasicParameter* PMIS_HgCdTe_BasicParameter_Default (const PMIS_Environment& env)
  {
    return new GSS_HgCdTe_BasicParameter(env);
  }
}
