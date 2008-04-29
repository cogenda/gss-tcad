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
// Material Type: Si(1-x)Ge(x). (in fact, use parameters for Si instead.)


#include "PMI.h"

class GSS_Mob_Analytic : public PMIS_Mobility
{
private:
  // parameters for Analytic mobility
  PetscScalar MUN_MIN ;
  PetscScalar MUN_MAX ;
  PetscScalar NREFN   ;
  PetscScalar NUN     ;
  PetscScalar XIN     ;
  PetscScalar ALPHAN  ;
  PetscScalar MUP_MIN ;
  PetscScalar MUP_MAX ;
  PetscScalar NREFP   ;
  PetscScalar NUP     ;
  PetscScalar XIP     ;
  PetscScalar ALPHAP  ;
  PetscScalar T300    ;
  // parameters for high field modification
  PetscScalar BETAN;
  PetscScalar BETAP;
  void Mob_Analytic_Init()
  {
    MUN_MIN = 5.524000E+01*cm*cm/V/s;
    MUN_MAX = 1.429230E+03*cm*cm/V/s;
    NREFN   = 1.072000E+17*pow(cm,-3);
    NUN     = -2.300000E+00;
    XIN     = -3.800000E+00;
    ALPHAN  = 7.300000E-01;
    MUP_MIN = 4.970000E+01*cm*cm/V/s;
    MUP_MAX = 4.793700E+02*cm*cm/V/s;
    NREFP   = 1.606000E+17*pow(cm,-3);
    NUP     = -2.200000E+00;
    XIP     = -3.700000E+00;
    ALPHAP  = 7.000000E-01;
    T300    = 300.0*K;
    BETAN   =  2.000000E+00;
    BETAP   =  1.000000E+00;
  }

public:
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUN_MIN+(MUN_MAX*pow(Tl/T300,NUN)-MUN_MIN)/ \
           (1+pow(Tl/T300,XIN)*pow((Na+Nd)/NREFN,ALPHAN));
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUN_MIN+(MUN_MAX*pow(Tl/T300,NUN)-MUN_MIN)/ \
           (1+pow(Tl/T300,XIN)*pow((Na+Nd)/NREFN,ALPHAN));
  }
  
  //---------------------------------------------------------------------------
  // Hole low field mobility
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUP_MIN+(MUP_MAX*pow(Tl/T300,NUP)-MUP_MIN)/ \
           (1+pow(Tl/T300,XIP)*pow((Na+Nd)/NREFP,ALPHAP));
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUP_MIN+(MUP_MAX*pow(Tl/T300,NUP)-MUP_MIN)/ \
           (1+pow(Tl/T300,XIP)*pow((Na+Nd)/NREFP,ALPHAP));
  }


public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    PetscScalar mu0  = ElecMobLowField(Tl);
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    AutoDScalar mu0  = ElecMobLowField(Tl);
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    PetscScalar mu0  = HoleMobLowField(Tl);
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    AutoDScalar mu0  = HoleMobLowField(Tl);
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }

  
// constructor and destructor
public:  
  GSS_Mob_Analytic(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Analytic_Init();
  }


  ~GSS_Mob_Analytic()
  {
  }

}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  use Analytic model as default mobility model
 */
extern "C"
{
  PMIS_Mobility* PMIS_SiGe_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
/* alias */
extern "C"
{
  PMIS_Mobility* PMIS_SiGe_Mob_Analytic (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
