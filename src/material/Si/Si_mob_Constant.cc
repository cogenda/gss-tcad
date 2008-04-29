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
// Material Type: Silicon


#include "PMI.h"

class GSS_Mob_Constant : public PMIS_Mobility
{
private:
  // parameters for constant mobility
  PetscScalar MUN ;
  PetscScalar MUP ;
  void Mob_Constant_Init()
  {
    MUN = 1.400E+03*cm*cm/V/s;
    MUP = 4.500E+02*cm*cm/V/s;
  }

  
public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
     return MUN;
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    return MUN;
  }
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
     return MUP;
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    return MUP;
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
  PMIS_Mobility* PMIS_Si_Mob_Constant (const PMIS_Environment& env)
  {
    return new GSS_Mob_Constant(env);
  }
}
