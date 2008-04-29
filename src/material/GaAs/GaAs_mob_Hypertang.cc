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
// Material Type: GaAs


#include "PMI.h"

class GSS_Mob_Hypertang : public PMIS_Mobility
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
  PetscScalar E0N     ;
  PetscScalar E0P     ;
  void Mob_Hypertang_Init()
  {
    MUN_MIN =  0.000000E+00*cm*cm/V/s;
    MUN_MAX =  8.500000E+03*cm*cm/V/s;
    NREFN   =  1.690000E+17*pow(cm,-3);
    NUN     = -1.000000E+00;
    XIN     =  0.000000E+00;
    ALPHAN  =  4.360000E-01;
    MUP_MIN =  0.000000E+00*cm*cm/V/s;
    MUP_MAX =  4.000000E+02*cm*cm/V/s;
    NREFP   =  2.750000E+17*pow(cm,-3);
    NUP     = -2.100000E+00;
    XIP     =  0.000000E+00;
    ALPHAP  =  3.950000E-01;
    T300    =  300.0*K;
    E0N     =  4.000000E+03*V/cm;
    E0P     =  1.000000E+06*V/cm;
  }

public:

  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUN_MIN+(MUN_MAX*pow(Tl/T300,NUN)-MUN_MIN)/ \
           (1+pow(Tl/T300,XIN)*pow((Na+Nd)/NREFN,ALPHAN));
  }

  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUP_MIN+(MUP_MAX*pow(Tl/T300,NUP)-MUP_MIN)/ \
           (1+pow(Tl/T300,XIP)*pow((Na+Nd)/NREFP,ALPHAP));
  }
  
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return MUN_MIN+(MUN_MAX*pow(Tl/T300,NUN)-MUN_MIN)/ \
           (1+pow(Tl/T300,XIN)*pow((Na+Nd)/NREFN,ALPHAN));
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
    PetscScalar vsat = 11.3E6*cm/s - 1.2E4*cm/s*Tl;
    PetscScalar mu0  = ElecMobLowField(Tl);
    if(Ep < 1e2*V/cm)    return mu0;
    return vsat/Ep*tanh(mu0*Ep/vsat);
  }

  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = 11.3E6*cm/s - 1.2E4*cm/s*Tl;
    AutoDScalar mu0  = ElecMobLowField(Tl);
    if(Ep < 1e2*V/cm)    return mu0;
    return vsat/Ep*tanh(mu0*Ep/vsat);
  }
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = 11.3E6*cm/s - 1.2E4*cm/s*Tl;
    PetscScalar mu0  = HoleMobLowField(Tl);
    if(Ep < 1e2*V/cm)    return mu0;
    return vsat/Ep*tanh(mu0*Ep/vsat);
  }

  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = 11.3E6*cm/s - 1.2E4*cm/s*Tl;
    AutoDScalar mu0  = HoleMobLowField(Tl);
    if(Ep < 1e2*V/cm)    return mu0;
    return vsat/Ep*tanh(mu0*Ep/vsat);
  }
 


// constructor and destructor
public:  
  GSS_Mob_Hypertang(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Hypertang_Init();
  }


  ~GSS_Mob_Hypertang()
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
  PMIS_Mobility* PMIS_GaAs_Mob_Hypertang (const PMIS_Environment& env)
  {
    return new GSS_Mob_Hypertang(env);
  }
}
