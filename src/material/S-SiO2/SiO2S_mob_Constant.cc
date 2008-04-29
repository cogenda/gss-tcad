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
// Material Type: SiO2 as semicondcutor


#include "PMI.h"

class GSS_Mob_Constant : public PMIS_Mobility
{
private:
  // parameters for constant mobility
  PetscScalar mumaxn ;
  PetscScalar mumaxp ;
  PetscScalar Exponentn ;
  PetscScalar Exponentp ;
  PetscScalar T0;
  void Mob_Constant_Init()
  {
    mumaxn = 0.05*cm*cm/V/s;
    mumaxp = 4.7050e+02*cm*cm/V/s;
    Exponentn = 2.5;
    Exponentp = 2.2;
    T0 = 300*K;
  }
  
public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
     return mumaxn*pow(Tl/T0,-Exponentn);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
     return mumaxn*pow(Tl/T0,-Exponentn);
  }
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
     return mumaxp*pow(Tl/T0,-Exponentp);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
     return mumaxp*pow(Tl/T0,-Exponentp);
  }
  
// constructor
public:  
  GSS_Mob_Constant(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Constant_Init();
  }


  ~GSS_Mob_Constant()
  {
  }

}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 */
extern "C"
{
  PMIS_Mobility* PMIS_SiO2S_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_Mob_Constant(env);
  }
}
