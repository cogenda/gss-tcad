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
  PetscScalar E0N     ;
  PetscScalar E0P     ;
  void Mob_Analytic_Init()
  {
    //Source: data base of DESSIS
    Ar_mumin_n	= 88*cm*cm/V/s;
    Ar_mumin_p	= 54.3*cm*cm/V/s;
    Ar_alm_n  	= -6.70E-01;
    Ar_alm_p    = -5.70E-01;
    Ar_mud_n  	= 2.20E+03*cm*cm/V/s;
    Ar_mud_p  	= 4.07E+02*cm*cm/V/s;
    Ar_ald_n  	= -4.00E+00;
    Ar_ald_p  	= -2.23E+00;
    Ar_N0_n   	= 1.25E+17*pow(cm,-3);
    Ar_N0_p   	= 2.35E+17*pow(cm,-3);
    Ar_alN_n  	= 1.9;
    Ar_alN_p  	= 2.4;
    Ar_a_n    	= 0.98;
    Ar_a_p    	= 0.88;
    Ar_ala_n  	= -1.50E-01;
    Ar_ala_p  	= -1.46E-01;
 	   
    T300        = 300*K;
    E0N         = 7.1000e+04*V/cm;
    E0P         = 7.1000e+04*V/cm;
  }

private:
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar AA = Ar_a_n*pow(Tl/T300,Ar_ala_n);
    PetscScalar N00=Ar_N0_n*pow(Tl/T300,Ar_alN_n); 
    PetscScalar muminA=Ar_mumin_n*pow(Tl/T300,Ar_alm_n); 
    PetscScalar mudA = Ar_mud_n*pow(Tl/T300,Ar_ald_n);
    PetscScalar mu_dop = muminA + mudA/(1.0+pow((Na+Nd)/N00,AA));
    return mu_dop;
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    AutoDScalar AA = Ar_a_n*pow(Tl/T300,Ar_ala_n);
    AutoDScalar N00=Ar_N0_n*pow(Tl/T300,Ar_alN_n); 
    AutoDScalar muminA=Ar_mumin_n*pow(Tl/T300,Ar_alm_n); 
    AutoDScalar mudA = Ar_mud_n*pow(Tl/T300,Ar_ald_n);
    AutoDScalar mu_dop = muminA + mudA/(1.0+pow((Na+Nd)/N00,AA));
    return mu_dop;
  }
  
  //---------------------------------------------------------------------------
  // Hole low field mobility
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar AA = Ar_a_p*pow(Tl/T300,Ar_ala_p);
    PetscScalar N00=Ar_N0_p*pow(Tl/T300,Ar_alN_p); 
    PetscScalar muminA=Ar_mumin_p*pow(Tl/T300,Ar_alm_p); 
    PetscScalar mudA = Ar_mud_p*pow(Tl/T300,Ar_ald_p);
    PetscScalar mu_dop = muminA + mudA/(1.0+pow((Na+Nd)/N00,AA));
    return mu_dop;
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    AutoDScalar AA = Ar_a_p*pow(Tl/T300,Ar_ala_p);
    AutoDScalar N00=Ar_N0_p*pow(Tl/T300,Ar_alN_p); 
    AutoDScalar muminA=Ar_mumin_p*pow(Tl/T300,Ar_alm_p); 
    AutoDScalar mudA = Ar_mud_p*pow(Tl/T300,Ar_ald_p);
    AutoDScalar mu_dop = muminA + mudA/(1.0+pow((Na+Nd)/N00,AA));
    return mu_dop;
  }
  

  
  
public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar vsat = 2.600e+07*cm/s;
    PetscScalar mu0  = ElecMobLowField(Tl);
    PetscScalar E = Ep>0 ? Ep : 0 ;
    return (mu0+vsat*pow(E,3)/pow(E0N,4))/(1+pow(E/E0N,4));
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    PetscScalar vsat = 2.600e+07*cm/s;
    AutoDScalar mu0  = ElecMobLowField(Tl);
    AutoDScalar E = fmax(Ep, 0.0) ;
    return (mu0+vsat*pow(E,3)/pow(E0N,4))/(1+pow(E/E0N,4));
  }
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = 2.600e+07*cm/s;
    PetscScalar mu0  = HoleMobLowField(Tl);
    PetscScalar E = Ep>0 ? Ep : 0 ;
    return (mu0+vsat*pow(E,3)/pow(E0P,4))/(1+pow(E/E0P,4));
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    PetscScalar vsat = 2.600e+07*cm/s;
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
  PMIS_Mobility* PMIS_InN_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
/* alias */
extern "C"
{
  PMIS_Mobility* PMIS_InN_Mob_Analytic (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
