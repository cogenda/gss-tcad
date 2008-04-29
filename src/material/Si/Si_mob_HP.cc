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
/*  Last update: July 23, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Silicon


#include "PMI.h"

class GSS_Mob_HP : public PMIS_Mobility
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
  // parameters for HP high field modification
  PetscScalar   MUN0_HP;    
  PetscScalar   ECN_HP;     
  PetscScalar   VSN_HP;    
  PetscScalar   VCN_HP;  
  PetscScalar   GN_HP;   
  PetscScalar   MUP0_HP;
  PetscScalar   ECP_HP;
  PetscScalar   VSP_HP;
  PetscScalar   VCP_HP;
  PetscScalar   GP_HP;
  PetscScalar   NREF_HP;
    
  void Mob_HP_Init()
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
    
    MUN0_HP = 7.740000E+02*cm*cm/V/s;      
    ECN_HP  = 5.500000E+05*V/cm;      
    VSN_HP  = 1.036000E+07*cm/s;         
    VCN_HP  = 4.900000E+06*cm/s;       
    GN_HP   = 8.800000E+00;        
    MUP0_HP = 2.500000E+02*cm*cm/V/s;      
    ECP_HP  = 2.780000E+05*V/cm;         
    VSP_HP  = 1.200000E+07*cm/s;       
    VCP_HP  = 2.928000E+06*cm/s;       
    GP_HP   = 1.600000E+00;
    NREF_HP = 5.000000E+17*pow(cm,-3);
  }

public:
  //---------------------------------------------------------------------------
  // Electron low field mobility, Analytic model
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
  // Hole low field mobility, Analytic model
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
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    PetscScalar mu0;
    if(Ntotal > NREF_HP) mu0=ElecMobLowField(Tl);
    else mu0 = MUN0_HP/(1+Et/ECN_HP);
    PetscScalar alpha = mu0*Ep/VCN_HP;
    PetscScalar beta  = mu0*Ep/VSN_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GN_HP)+beta*beta);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    AutoDScalar mu0;
    if(Ntotal > NREF_HP) mu0=ElecMobLowField(Tl);
    else mu0 = MUN0_HP/(1+Et/ECN_HP);
    AutoDScalar alpha = mu0*Ep/VCN_HP;
    AutoDScalar beta  = mu0*Ep/VSN_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GN_HP)+beta*beta);
  }
  
  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, 
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    PetscScalar mu0;
    if(Ntotal > NREF_HP) mu0=HoleMobLowField(Tl);
    else mu0 = MUP0_HP/(1+Et/ECP_HP);
    PetscScalar alpha = mu0*Ep/VCP_HP;
    PetscScalar beta  = mu0*Ep/VSP_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GP_HP)+beta*beta);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, 
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    PetscScalar Ntotal = ReadDopingNa()+ReadDopingNd();
    AutoDScalar mu0;
    if(Ntotal > NREF_HP) mu0=HoleMobLowField(Tl);
    else mu0 = MUP0_HP/(1+Et/ECP_HP);
    AutoDScalar alpha = mu0*Ep/VCP_HP;
    AutoDScalar beta  = mu0*Ep/VSP_HP;
    return mu0/sqrt(1+alpha*alpha/(alpha+GP_HP)+beta*beta);
  }
  
// constructor
public:  
  GSS_Mob_HP(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_HP_Init();
  }


  ~GSS_Mob_HP()
  {
  }

}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  use Analytic model as default mobility model
 */
/* alias */
extern "C"
{
  PMIS_Mobility* PMIS_Si_Mob_HP (const PMIS_Environment& env)
  {
    return new GSS_Mob_HP(env);
  }
}
