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

class GSS_Mob_Lucent : public PMIS_Mobility
{
private:
  // parameters for Lucent mobility
  PetscScalar MMNN_UM;
  PetscScalar MMXN_UM;
  PetscScalar NRFN_UM;
  PetscScalar ALPN_UM;
  PetscScalar TETN_UM;
  PetscScalar NRFD_UM;
  PetscScalar CRFD_UM;

  PetscScalar MMNP_UM;
  PetscScalar MMXP_UM;
  PetscScalar NRFP_UM;
  PetscScalar ALPP_UM;
  PetscScalar TETP_UM;
  PetscScalar NRFA_UM;
  PetscScalar CRFA_UM;
  PetscScalar NSC_REF;
  PetscScalar CAR_REF;
  PetscScalar me_over_m0;
  PetscScalar mh_over_m0;
  PetscScalar me_over_mh;
  // temperature
  PetscScalar T300;
  // parameters for transvers field modification
  PetscScalar AN_LUC;
  PetscScalar AP_LUC;
  PetscScalar BN_LUC;
  PetscScalar BP_LUC;
  PetscScalar CN_LUC;
  PetscScalar CP_LUC;
  //PetscScalar DN_LUC;
  //PetscScalar DP_LUC;
  PetscScalar FN_LUC;
  PetscScalar FP_LUC;
  PetscScalar KN_LUC;
  PetscScalar KP_LUC;
  PetscScalar EXN4_LUC;
  PetscScalar EXP4_LUC;
  PetscScalar EXN9_LUC;
  PetscScalar EXP9_LUC;
  // parameters for parallel field modification
  PetscScalar BETAN;
  PetscScalar BETAP;

  void Mob_Lucent_Init()
  {
    MMNN_UM = 5.220000E+01*cm*cm/V/s;
    MMXN_UM = 1.417000E+03*cm*cm/V/s;
    NRFN_UM = 9.680000E+16*pow(cm,-3);
    ALPN_UM = 6.800000E-01;
    TETN_UM = 2.285000E+00;
    NRFD_UM = 4.000000E+20*pow(cm,-3);
    CRFD_UM = 2.100000E-01;

    MMNP_UM = 4.490000E+01*cm*cm/V/s;
    MMXP_UM = 4.705000E+02*cm*cm/V/s;
    NRFP_UM = 2.230000E+17*pow(cm,-3);
    ALPP_UM = 7.190000E-01;
    TETP_UM = 2.247000E+00;
    NRFA_UM = 7.200000E+20*pow(cm,-3);
    CRFA_UM = 5.000000E-01;

    NSC_REF = 3.97e13*pow(cm,-2);
    CAR_REF = 1.36e20*pow(cm,-3);
    me_over_m0 = 1.0;
    mh_over_m0 = 1.258;
    me_over_mh = 1.0/1.258;
    T300    = 300.0*K;

    EXN4_LUC = 2.330000E-02;
    EXP4_LUC = 1.190000E-02;
    EXN9_LUC = 7.670000E-02;
    EXP9_LUC = 1.230000E-01;
    AN_LUC   = 2.580000E+00;
    AP_LUC   = 2.180000E+00;
    BN_LUC   = 3.610000E+07*cm/s;
    BP_LUC   = 1.510000E+07*cm/s;
    CN_LUC   = 1.700000E+04*cm*cm/V/s*pow(V/cm,PetscScalar(1.0/3.0))*pow(cm,3*EXN4_LUC);
    CP_LUC   = 4.180000E+03*cm*cm/V/s*pow(V/cm,PetscScalar(1.0/3.0))*pow(cm,3*EXP4_LUC);
    //DN_LUC=3.580000E+18*cm*cm/V/s*pow(V/cm,rn);
    //DP_LUC=4.100000E+15*cm*cm/V/s*pow(V/cm,rp);
    FN_LUC   = 6.850000E-21*pow(cm,3*(1-EXN9_LUC));
    FP_LUC   = 7.820000E-21*pow(cm,3*(1-EXP9_LUC));
    KN_LUC   = 1.700000E+00;
    KP_LUC   = 9.000000E-01;

    BETAN = 2.000000E+00;
    BETAP = 1.000000E+00;

  }
  //---------------------------------------------------------------------------
  // Electron low field bulk mobility
  PetscScalar ElecMobPhilips(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl) const
  {
    PetscScalar mu_lattice = MMXN_UM*pow(Tl/T300,-TETN_UM);
    PetscScalar mu1 = MMXN_UM*MMXN_UM/(MMXN_UM-MMNN_UM)*pow(Tl/T300,3*ALPN_UM-1.5);
    PetscScalar mu2 = MMXN_UM*MMNN_UM/(MMXN_UM-MMNN_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    PetscScalar Nsc = Nds+Nas+fabs(p);

    PetscScalar P   = 1.0/(2.459/(NSC_REF/pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*me_over_m0))*(Tl/T300)*(Tl/T300);
    PetscScalar pp1 = pow(P,PetscScalar(0.6478));
    PetscScalar F   = (0.7643*pp1+2.2999+6.5502*me_over_mh)/(pp1+2.3670-0.8552*me_over_mh);
    //PetscScalar G   = 1-0.89233/pow(0.41372+P*pow(Tl/T300/me_over_m0,PetscScalar(0.28227)),PetscScalar(0.19778))+0.005978/pow(P*pow(Tl/T300/me_over_m0,PetscScalar(0.72169)),PetscScalar(1.80618));
    PetscScalar G   = 1-4.41804/pow(39.9014+P*pow(Tl/T300/me_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/pow(P*pow(T300/Tl*me_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    PetscScalar Nsce = Nds+Nas*G+fabs(p)/F;
    PetscScalar mu_scatt = mu1*(Nsc/Nsce)*pow(NRFN_UM/Nsc,ALPN_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }
  AutoDScalar ElecMobPhilips(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl) const
  {
    AutoDScalar mu_lattice = MMXN_UM*pow(Tl/T300,-TETN_UM);
    AutoDScalar mu1 = MMXN_UM*MMXN_UM/(MMXN_UM-MMNN_UM)*pow(Tl/T300,3*ALPN_UM-1.5);
    AutoDScalar mu2 = MMXN_UM*MMNN_UM/(MMXN_UM-MMNN_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    AutoDScalar Nsc = Nds+Nas+fabs(p);

    AutoDScalar P   = 1.0/(2.459/(NSC_REF/pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*me_over_m0))*(Tl/T300)*(Tl/T300);
    AutoDScalar pp1 = pow(P,PetscScalar(0.6478));
    AutoDScalar F   = (0.7643*pp1+2.2999+6.5502*me_over_mh)/(pp1+2.3670-0.8552*me_over_mh);
    //PetscScalar G   = 1-0.89233/pow(0.41372+P*pow(Tl/T300/me_over_m0,PetscScalar(0.28227)),PetscScalar(0.19778))+0.005978/pow(P*pow(Tl/T300/me_over_m0,PetscScalar(0.72169)),PetscScalar(1.80618));
    AutoDScalar G   = 1-4.41804/pow(39.9014+P*pow(Tl/T300/me_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/pow(P*pow(T300/Tl*me_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    AutoDScalar Nsce = Nds+Nas*G+fabs(p)/F;
    AutoDScalar mu_scatt = mu1*(Nsc/Nsce)*pow(NRFN_UM/Nsc,ALPN_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }

  //---------------------------------------------------------------------------
  // Hole low field bulk mobility
  PetscScalar HoleMobPhilips(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl) const
  {
    PetscScalar mu_lattice = MMXP_UM*pow(Tl/T300,-TETP_UM);
    PetscScalar mu1 = MMXP_UM*MMXP_UM/(MMXP_UM-MMNP_UM)*pow(Tl/T300,3*ALPP_UM-1.5);
    PetscScalar mu2 = MMXP_UM*MMNP_UM/(MMXP_UM-MMNP_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    PetscScalar Nsc = Nds+Nas+fabs(n);

    PetscScalar P   = 1.0/(2.459/(NSC_REF/pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*mh_over_m0))*(Tl/T300)*(Tl/T300);
    PetscScalar pp1 = pow(P,PetscScalar(0.6478));
    PetscScalar F   = (0.7643*pp1+2.2999+6.5502/me_over_mh)/(pp1+2.3670-0.8552/me_over_mh);
    //PetscScalar G   = 1-0.89233/pow(0.41372+P*pow(Tl/T300/mh_over_m0,0.28227),0.19778)+0.005978/pow(P*pow(Tl/T300/mh_over_m0,0.72169),1.80618);
    PetscScalar G   = 1-4.41804/pow(39.9014+P*pow(Tl/T300/mh_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/pow(P*pow(T300/Tl*mh_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    PetscScalar Nsce = Nas+Nds*G+fabs(n)/F;
    PetscScalar mu_scatt = mu1*(Nsc/Nsce)*pow(NRFP_UM/Nsc,ALPP_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }
  AutoDScalar HoleMobPhilips(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl) const
  {
    AutoDScalar mu_lattice = MMXP_UM*pow(Tl/T300,-TETP_UM);
    AutoDScalar mu1 = MMXP_UM*MMXP_UM/(MMXP_UM-MMNP_UM)*pow(Tl/T300,3*ALPP_UM-1.5);
    AutoDScalar mu2 = MMXP_UM*MMNP_UM/(MMXP_UM-MMNP_UM)*sqrt(T300/Tl);
    PetscScalar Na  = ReadDopingNa()+1e0*pow(cm,-3);
    PetscScalar Nd  = ReadDopingNd()+1e0*pow(cm,-3);
    PetscScalar Nds = Nd*(1.0+1.0/(CRFD_UM+(NRFD_UM/Nd)*(NRFD_UM/Nd)));
    PetscScalar Nas = Na*(1.0+1.0/(CRFA_UM+(NRFA_UM/Na)*(NRFA_UM/Na)));
    AutoDScalar Nsc = Nds+Nas+fabs(n);

    AutoDScalar P   = 1.0/(2.459/(NSC_REF/pow(Nsc,PetscScalar(2.0/3.0)))+3.828/(CAR_REF/fabs(n+p)*mh_over_m0))*(Tl/T300)*(Tl/T300);
    AutoDScalar pp1 = pow(P,PetscScalar(0.6478));
    AutoDScalar F   = (0.7643*pp1+2.2999+6.5502/me_over_mh)/(pp1+2.3670-0.8552/me_over_mh);
    //PetscScalar G   = 1-0.89233/pow(0.41372+P*pow(Tl/T300/mh_over_m0,0.28227),0.19778)+0.005978/pow(P*pow(Tl/T300/mh_over_m0,0.72169),1.80618);
    AutoDScalar G   = 1-4.41804/pow(39.9014+P*pow(Tl/T300/mh_over_m0,PetscScalar(0.0001)),PetscScalar(0.38297))+0.52896/pow(P*pow(T300/Tl*mh_over_m0,PetscScalar(1.595787)),PetscScalar(0.25948));
    AutoDScalar Nsce = Nas+Nds*G+fabs(n)/F;
    AutoDScalar mu_scatt = mu1*(Nsc/Nsce)*pow(NRFP_UM/Nsc,ALPP_UM)+mu2*(fabs(n+p)/Nsce);
    return 1.0/(1.0/mu_lattice+1.0/mu_scatt);
  }

  //---------------------------------------------------------------------------
  //the surface acoustical phono scattering and roughness scattering for electron
  PetscScalar ElecMobLombardi(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BN_LUC/ET + CN_LUC*pow(N_total,EXN4_LUC)*pow(Tl/T300,-KN_LUC)*pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar r = AN_LUC + FN_LUC*fabs(n+p)/pow(N_total,EXN9_LUC);
    PetscScalar DN_LUC=3.580000E+18*cm*cm/V/s*pow(V/cm,r);
    PetscScalar mu_sr = DN_LUC*pow(ET,-r);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar ElecMobLombardi(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BN_LUC/ET + CN_LUC*pow(N_total,EXN4_LUC)*pow(Tl/T300,-KN_LUC)*pow(ET,PetscScalar(-1.0/3.0));
    AutoDScalar r = AN_LUC + FN_LUC*fabs(n+p)/pow(N_total,EXN9_LUC);
    AutoDScalar DN_LUC=3.580000E+18*cm*cm/V/s*pow(V/cm,r);
    AutoDScalar mu_sr = DN_LUC*pow(ET,-r);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }

  //---------------------------------------------------------------------------
  //the surface acoustical phono scattering and roughness scattering for hole
  PetscScalar HoleMobLombardi(const PetscScalar &p,const PetscScalar &n,const PetscScalar &Tl,const PetscScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    PetscScalar ET = Et+1.0*V/cm;
    PetscScalar mu_ac = BP_LUC/ET + CP_LUC*pow(N_total,EXP4_LUC)*pow(Tl/T300,-KP_LUC)*pow(ET,PetscScalar(-1.0/3.0));
    PetscScalar r = AP_LUC + FP_LUC*fabs(n+p)/pow(N_total,EXP9_LUC);
    PetscScalar DP_LUC=4.100000E+15*cm*cm/V/s*pow(V/cm,r);
    PetscScalar mu_sr = DP_LUC*pow(ET,-r);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }
  AutoDScalar HoleMobLombardi(const AutoDScalar &p,const AutoDScalar &n,const AutoDScalar &Tl,const AutoDScalar &Et) const
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N_total = Na+Nd+1e0*pow(cm,-3);
    AutoDScalar ET = Et+1.0*V/cm;
    AutoDScalar mu_ac = BP_LUC/ET + CP_LUC*pow(N_total,EXP4_LUC)*pow(Tl/T300,-KP_LUC)*pow(ET,PetscScalar(-1.0/3.0));
    AutoDScalar r = AP_LUC + FP_LUC*fabs(n+p)/pow(N_total,EXP9_LUC);
    AutoDScalar DP_LUC=4.100000E+15*cm*cm/V/s*pow(V/cm,r);
    AutoDScalar mu_sr = DP_LUC*pow(ET,-r);
    return 1.0/(1.0/mu_ac+1.0/mu_sr);
  }

public:

  //---------------------------------------------------------------------------
  // Electron mobility
  PetscScalar ElecMob(const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                      const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const
  {
    PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/ElecMobPhilips(p,n,Tl)+1.0/ElecMobLombardi(p,n,Tl,Et));
    PetscScalar mu = 2*mu0/(1+pow(1+pow(2*mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN));
    return mu;
  }
  AutoDScalar ElecMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const
  {
    AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/ElecMobPhilips(p,n,Tl)+1.0/ElecMobLombardi(p,n,Tl,Et));
    AutoDScalar mu = 2*mu0/(1+pow(1+pow(2*mu0*fabs(Ep)/vsat,BETAN),1.0/BETAN));
    return mu;
  }

  //---------------------------------------------------------------------------
  // Hole mobility
  PetscScalar HoleMob (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl,
                       const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const
  {
    PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    PetscScalar mu0  = 1.0/(1.0/HoleMobPhilips(p,n,Tl)+1.0/HoleMobLombardi(p,n,Tl,Et));
    PetscScalar mu = 2*mu0/(1+pow(1+pow(2*mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP));
    return mu;
  }
  AutoDScalar HoleMob(const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl,
                      const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const
  {
    AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
    AutoDScalar mu0  = 1.0/(1.0/HoleMobPhilips(p,n,Tl)+1.0/HoleMobLombardi(p,n,Tl,Et));
    AutoDScalar mu = 2*mu0/(1+pow(1+pow(2*mu0*fabs(Ep)/vsat,BETAP),1.0/BETAP));
    return mu;
  }


  // constructor
public:
  GSS_Mob_Lucent(const PMIS_Environment &env):PMIS_Mobility(env)
  {
    Mob_Lucent_Init();
  }

  ~GSS_Mob_Lucent(){}}
;

/*---------------------------------------------------------------
 *  the interface function called by material databse controller
 *  and it setup Lucent mobility model
 */
extern "C"
{
  PMIS_Mobility* PMIS_SiGe_Mob_Lucent (const PMIS_Environment& env)
  {
    return new GSS_Mob_Lucent(env);
  }
}
