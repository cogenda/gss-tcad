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
/*  Last update: July 26, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//
// Material Type: Si(1-x)Ge(x). (in fact, use parameters for Si instead.)


#include "PMI.h"

class GSS_Mob_Lombardi : public PMIS_Mobility
{
private:
  PetscScalar EXN1_LSM;
  PetscScalar EXN2_LSM;
  PetscScalar EXN3_LSM;
  PetscScalar EXN4_LSM;
  PetscScalar EXN8_LSM;
  PetscScalar MUN0_LSM;
  PetscScalar MUN1_LSM;
  PetscScalar MUN2_LSM;
  PetscScalar CRN_LSM;
  PetscScalar CSN_LSM;
  PetscScalar BN_LSM;
  PetscScalar CN_LSM;
  PetscScalar DN_LSM;

  PetscScalar EXP1_LSM;
  PetscScalar EXP2_LSM;
  PetscScalar EXP3_LSM;
  PetscScalar EXP4_LSM;
  PetscScalar EXP8_LSM;
  PetscScalar MUP0_LSM;
  PetscScalar MUP1_LSM;
  PetscScalar MUP2_LSM;
  PetscScalar CRP_LSM;
  PetscScalar CSP_LSM;
  PetscScalar BP_LSM;
  PetscScalar CP_LSM;
  PetscScalar DP_LSM;
  PetscScalar PC_LSM;

  // parameters for parallel field modification
  PetscScalar BETAN;
  PetscScalar BETAP;
  PetscScalar T300;

  void Mob_Lombardi_Init()
  {
    EXN1_LSM=6.800000E-01;
    EXN2_LSM=2.000000E+00;
    EXN3_LSM=2.500000E+00;
    EXN4_LSM=1.250000E-01;
    EXN8_LSM=2.000000E+00;
    MUN0_LSM=5.220000E+01*cm*cm/V/s;
    MUN1_LSM=4.340000E+01*cm*cm/V/s;
    MUN2_LSM=1.417000E+03*cm*cm/V/s;
    CRN_LSM=9.680000E+16*pow(cm,-3);
    CSN_LSM=3.430000E+20*pow(cm,-3);
    BN_LSM=4.750000E+07*cm/s;
    CN_LSM=1.740000E+05*K*cm/s*pow(V/cm,PetscScalar(-2.0/3.0))*pow(cm,3*EXN4_LSM);
    DN_LSM=5.820000E+14*cm*cm/V/s*pow(V/cm,EXN8_LSM);

    EXP1_LSM=7.190000E-01;
    EXP2_LSM=2.000000E+00;
    EXP3_LSM=2.200000E+00;
    EXP4_LSM=3.170000E-02;
    EXP8_LSM=2.000000E+00;
    MUP0_LSM=4.490000E+01*cm*cm/V/s;
    MUP1_LSM=2.900000E+01*cm*cm/V/s;
    MUP2_LSM=4.705000E+02*cm*cm/V/s;
    CRP_LSM=2.230000E+17*pow(cm,-3);
    CSP_LSM=6.100000E+20*pow(cm,-3);
    BP_LSM=9.930000E+06*cm/s;
    CP_LSM=8.840000E+05*K*cm/s*pow(V/cm,PetscScalar(-2.0/3.0))*pow(cm,3*EXP4_LSM);;
    DP_LSM=2.050000E+14*cm*cm/V/s*pow(V/cm,EXN8_LSM);
    PC_LSM=9.230000E+16*pow(cm,-3);

    BETAN = 2.000000E+00;
    BETAP = 1.000000E+00;
    T300  = 300.0*K;
  }
  //---------------------------------------------------------------------------
  // Electron low field mobility
  PetscScalar ElecMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    PetscScalar mu_max = MUN2_LSM*pow(Tl/T300,-EXN3_LSM);
    return MUN0_LSM+(mu_max-MUN0_LSM)/(1+pow(N_total/CRN_LSM,EXN1_LSM))-MUN1_LSM/(1+pow(CSN_LSM/N_total,EXN2_LSM));
  }
  AutoDScalar ElecMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    AutoDScalar mu_max = MUN2_LSM*pow(Tl/T300,-EXN3_LSM);
    return MUN0_LSM+(mu_max-MUN0_LSM)/(1+pow(N_total/CRN_LSM,EXN1_LSM))-MUN1_LSM/(1+pow(CSN_LSM/N_total,EXN2_LSM));
  }

  //---------------------------------------------------------------------------
  // Hole low field mobility, Analytic model
  PetscScalar HoleMobLowField(const PetscScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    PetscScalar mu_max = MUP2_LSM*pow(Tl/T300,-EXP3_LSM);
    return MUP0_LSM*exp(-PC_LSM/N_total)+mu_max/(1+pow(N_total/CRP_LSM,EXP1_LSM))-MUP1_LSM/(1+pow(CSP_LSM/N_total,EXP2_LSM));
  }
  AutoDScalar HoleMobLowField(const AutoDScalar &Tl) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    AutoDScalar mu_max = MUP2_LSM*pow(Tl/T300,-EXP3_LSM);
    return MUP0_LSM*exp(-PC_LSM/N_total)+mu_max/(1+pow(N_total/CRP_LSM,EXP1_LSM))-MUP1_LSM/(1+pow(CSP_LSM/N_total,EXP2_LSM));
  }

  //---------------------------------------------------------------------------
  // Electron surface mobility, acoustical phono scattering and roughness scattering
  PetscScalar ElecMobSurface(const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BN_LSM/ET + CN_LSM*pow(N_total,EXN4_LSM)/Tl*pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar mu_sr = DN_LSM*pow(ET,-EXN8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar ElecMobSurface(const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BN_LSM/ET + CN_LSM*pow(N_total,EXN4_LSM)/Tl*pow(ET,PetscScalar(-1.0/3.0));
    AutoDScalar mu_sr = DN_LSM*pow(ET,-EXN8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }

  //---------------------------------------------------------------------------
  // Hole surface mobility, acoustical phono scattering and roughness scattering
  PetscScalar HoleMobSurface(const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BP_LSM/ET + CP_LSM*pow(N_total,EXP4_LSM)/Tl*pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar mu_sr = DN_LSM*pow(ET,-EXP8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar HoleMobSurface(const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BP_LSM/ET + CP_LSM*pow(N_total,EXP4_LSM)/Tl*pow(ET+1.0*V/cm,PetscScalar(-1.0/3.0));
    AutoDScalar mu_sr = DN_LSM*pow(ET,-EXP8_LSM);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }


public:

  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/ElecMobLowField(Tl)+1.0/ElecMobSurface(Tl,Et));
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/ElecMobLowField(Tl)+1.0/ElecMobSurface(Tl,Et));
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN);
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/HoleMobLowField(Tl)+1.0/HoleMobSurface(Tl,Et));
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/HoleMobLowField(Tl)+1.0/HoleMobSurface(Tl,Et));
    return mu0/pow(1+pow(mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP);
  }

// constructor
public:
  GSS_Mob_Lombardi(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Lombardi_Init();
  }

  ~GSS_Mob_Lombardi(){}
}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  and it setup Lucent mobility model
 */
extern "C"
{
  PMIS_Mobility* PMIS_SiGe_Mob_Lombardi (const PMIS_Environment& env)
  {
    return new GSS_Mob_Lombardi(env);
  }
}
