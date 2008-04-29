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
// Material Type: InP


#include "PMI.h"

class GSS_Mob_Analytic : public PMIS_Mobility
{

private:
  // parameters for Analytic mobility
  PetscScalar Ar_mumin_n;
  PetscScalar Ar_mumin_p;
  PetscScalar Ar_alm_n  ;  
  PetscScalar Ar_alm_p  ;   
  PetscScalar Ar_mud_n  ;  
  PetscScalar Ar_mud_p  ;  
  PetscScalar Ar_ald_n  ;  
  PetscScalar Ar_ald_p  ;  
  PetscScalar Ar_N0_n   ;
  PetscScalar Ar_N0_p   ;
  PetscScalar Ar_alN_n  ;  
  PetscScalar Ar_alN_p  ;  
  PetscScalar Ar_a_n    ;  
  PetscScalar Ar_a_p    ;  
  PetscScalar Ar_ala_n  ;  
  PetscScalar Ar_ala_p  ;  
  PetscScalar T300      ;
  // parameters for high field modification
  PetscScalar E0N       ;
  PetscScalar E0P       ;
  void Mob_Analytic_Init()
  {
    //Source: data base of DESSIS
    Ar_mumin_n	= 4.5000e+03*cm*cm/V/s;
    Ar_mumin_p	= 1.5000e+02*cm*cm/V/s;
    Ar_alm_n  	= -1.5000e+00;
    Ar_alm_p    = -1.5000e+00;
    Ar_mud_n  	= 0.0000e+00*cm*cm/V/s;
    Ar_mud_p  	= 0.0000e+00*cm*cm/V/s;
    Ar_ald_n  	= 0.0000e+00;
    Ar_ald_p  	= 0.0000e+00;
    Ar_N0_n   	= 1.0000e+17*pow(cm,-3);
    Ar_N0_p   	= 1.0000e+17*pow(cm,-3);
    Ar_alN_n  	= 0.0000e+00;
    Ar_alN_p  	= 0.0000e+00;
    Ar_a_n    	= 0.0000e+00;
    Ar_a_p    	= 0.0000e+00;
    Ar_ala_n  	= 0.0000e+00;
    Ar_ala_p  	= 0.0000e+00;
    T300        = 300*K;
    E0N         = 4.0000e+03*V/cm;
    E0P         = 4.0000e+03*V/cm;
  }

private:

  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    return Ar_mumin_n*pow(Tl/T300,Ar_alm_n);
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    return Ar_mumin_n*pow(Tl/T300,Ar_alm_n);
  }

  //---------------------------------------------------------------------------
  // Hole low field mobility
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    return Ar_mumin_p*pow(Tl/T300,Ar_alm_p);
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    return Ar_mumin_p*pow(Tl/T300,Ar_alm_p);
  }
  

public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar vsat = 1.0000e+07*cm/s;
    PetscScalar mu0  = ElecMobLowField(Tl);
    PetscScalar E = Ep>0 ? Ep : 0 ;
    return (mu0+vsat*pow(E,3)/pow(E0N,4))/(1+pow(E/E0N,4));
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    PetscScalar vsat = 1.0000e+07*cm/s;
    AutoDScalar mu0  = ElecMobLowField(Tl);
    AutoDScalar E = fmax(Ep, 0.0) ;
    return (mu0+vsat*pow(E,3)/pow(E0N,4))/(1+pow(E/E0N,4));
  }
  
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = 1.0000e+07*cm/s;
    PetscScalar mu0  = HoleMobLowField(Tl);
    PetscScalar E = Ep>0 ? Ep : 0 ;
    return (mu0+vsat*pow(E,3)/pow(E0P,4))/(1+pow(E/E0P,4));
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    PetscScalar vsat = 1.0000e+07*cm/s;
    AutoDScalar mu0  = HoleMobLowField(Tl);
    AutoDScalar E = fmax(Ep, 0.0) ;
    return (mu0+vsat*pow(E,3)/pow(E0P,4))/(1+pow(E/E0P,4));
  }

  
// constructor and destructor
public:    
  GSS_Mob_Analytic(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Analytic_Init();
  }


  ~GSS_Mob_Analytic()
  {}

}
;


/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  use Analytic model as default mobility model
 */
extern "C"
{
  PMIS_Mobility* PMIS_InP_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
/* alias */
extern "C"
{
  PMIS_Mobility* PMIS_InP_Mob_Analytic (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
