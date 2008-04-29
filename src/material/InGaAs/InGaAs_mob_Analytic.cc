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
// Material Type: InGaAs


#include "PMI.h"

class GSS_Mob_Analytic : public PMIS_Mobility
{
private:
  // parameters for Analytic mobility
  PetscScalar MUN_MIN;
  PetscScalar MIN_X1;
  PetscScalar MIN_X2;
  PetscScalar MUN_MAX;
  PetscScalar MAN_X1;
  PetscScalar MAN_X2;
  PetscScalar NREFN;
  PetscScalar NREFN2;
  PetscScalar NUN;
  PetscScalar XIN;
  PetscScalar ALPHAN;
  PetscScalar MUP_MIN;
  PetscScalar MIP_X1;
  PetscScalar MIP_X2;
  PetscScalar MUP_MAX;
  PetscScalar MAP_X1;
  PetscScalar MAP_X2;
  PetscScalar NREFP;
  PetscScalar NREFP2;
  PetscScalar NUP;
  PetscScalar XIP;
  PetscScalar ALPHAP;
  PetscScalar T300;
  // parameters for high field modification
  PetscScalar VSATN;
  PetscScalar VSN1;
  PetscScalar VSN2;
  PetscScalar VSATP;
  PetscScalar VSP1;
  PetscScalar VSP2;
  PetscScalar E0N;
  PetscScalar EN1;
  PetscScalar EN2;
  PetscScalar E0P;

  void Mob_Analytic_Init()
  {
    MUN_MIN =  3.996900E+03*cm*cm/V/s;
    MIN_X1  = -2.229000E-01;
    MIN_X2  = -1.852000E-01;
    MUN_MAX =  2.725200E+04*cm*cm/V/s;
    MAN_X1  = -1.760000E+00;
    MAN_X2  =  1.122500E+00;
    NREFN   =  3.629700E+17*pow(cm,-3);
    NREFN2  =  1.745100E+18*pow(cm,-3);
    NUN     =  0.000000E+00;
    XIN     =  0.000000E+00;
    ALPHAN  =  1.000000E+00;

    MUP_MIN =  0.000000E+00*cm*cm/V/s;
    MIP_X1  =  0.000000E+00;
    MIP_X2  =  0.000000E+00;
    MUP_MAX =  4.800000E+02*cm*cm/V/s;
    MAP_X1  =  0.000000E+00;
    MAP_X2  =  0.000000E+00;
    NREFP   =  1.000000E+30*pow(cm,-3);
    NREFP2  =  1.000000E+30*pow(cm,-3);
    NUP     =  0.000000E+00;
    XIP     =  0.000000E+00;
    ALPHAP  =  3.950000E-01;
    T300    =  300.0*K;

    VSATN   =  1.000000E+07*cm/s;
    VSN1    = -1.019000E+00;
    VSN2    =  7.091700E-01;
    VSATP   =  0.000000E+00*cm/s;
    VSP1    =  0.000000E+00;
    VSP2    =  0.000000E+00;
    E0N     =  6.000000E+03*V/cm;
    EN1     = -7.956000E-01;
    EN2     =  6.299800E-01;
    E0P     =  1.000000E+06*V/cm;

  }

private:
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUN_MIN*(1+MIN_X1*x+MIN_X2*x*x);
    PetscScalar mu_max = MUN_MAX*(1+MAN_X1*x+MAN_X2*x*x);
    return mu_min+(mu_max*pow(Tl/T300,NUN)-mu_min)/ \
           (1+pow(Tl/T300,XIN)*(pow((Na+Nd)/NREFN,ALPHAN)+pow((Na+Nd)/NREFN2,3)));
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUN_MIN*(1+MIN_X1*x+MIN_X2*x*x);
    PetscScalar mu_max = MUN_MAX*(1+MAN_X1*x+MAN_X2*x*x);
    return mu_min+(mu_max*pow(Tl/T300,NUN)-mu_min)/ \
           (1+pow(Tl/T300,XIN)*(pow((Na+Nd)/NREFN,ALPHAN)+pow((Na+Nd)/NREFN2,3)));
  }
  
  
  //---------------------------------------------------------------------------
  // Hole low field mobility
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUP_MIN*(1+MIP_X1*x+MIP_X2*x*x);
    PetscScalar mu_max = MUP_MAX*(1+MAP_X1*x+MAP_X2*x*x);
    return mu_min+(mu_max*pow(Tl/T300,NUP)-mu_min)/ \
           (1+pow(Tl/T300,XIP)*(pow((Na+Nd)/NREFP,ALPHAP)+pow((Na+Nd)/NREFP2,3)));
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar x = ReadxMoleFraction();
    PetscScalar mu_min = MUP_MIN*(1+MIP_X1*x+MIP_X2*x*x);
    PetscScalar mu_max = MUP_MAX*(1+MAP_X1*x+MAP_X2*x*x);
    return mu_min+(mu_max*pow(Tl/T300,NUP)-mu_min)/ \
           (1+pow(Tl/T300,XIP)*(pow((Na+Nd)/NREFP,ALPHAP)+pow((Na+Nd)/NREFP2,3)));
  }
  

  
public:
  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar x = ReadxMoleFraction();
    PetscScalar vsat = VSATN*(1+ VSN1*x + VSN2*x*x);
    PetscScalar E0   = E0N*(1+EN1*x+EN2*x*x);
    PetscScalar mu0  = ElecMobLowField(Tl);
    return (mu0+vsat*pow(Ep,3)/pow(E0,4))/(1+pow(Ep/E0,4));
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    PetscScalar x = ReadxMoleFraction();
    PetscScalar vsat = VSATN*(1+ VSN1*x + VSN2*x*x);
    PetscScalar E0   = E0N*(1+EN1*x+EN2*x*x);
    AutoDScalar mu0  = ElecMobLowField(Tl);
    return (mu0+vsat*pow(Ep,3)/pow(E0,4))/(1+pow(Ep/E0,4));
  }
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    return HoleMobLowField(Tl);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    return HoleMobLowField(Tl);
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
  PMIS_Mobility* PMIS_InGaAs_Mob_Default (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
/* alias */
extern "C"
{
  PMIS_Mobility* PMIS_InGaAs_Mob_Analytic (const PMIS_Environment& env)
  {
    return new GSS_Mob_Analytic(env);
  }
}
