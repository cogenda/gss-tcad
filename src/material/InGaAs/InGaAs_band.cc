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


class GSS_InGaAs_BandStructure : public PMIS_BandStructure
{
private:
  PetscScalar T300;
  //[Bandgap]
  // Bandgap and Effective Density of States
  PetscScalar EG300;     // The energy bandgap of the material at 300 K.
  PetscScalar EG_X0;
  PetscScalar EG_X1;
  PetscScalar EG_X2;
  PetscScalar EG_X3;
  PetscScalar EG_X4;
  PetscScalar EG_X5;
  PetscScalar EG_X6;
  PetscScalar EG_X7;
  PetscScalar EG_X8;
  PetscScalar EG_X9;
  PetscScalar EG_X10;
  PetscScalar EG_X11;
  PetscScalar EG_X12;
  PetscScalar EG_X13;
  PetscScalar EG_X14;
  PetscScalar pm;
  PetscScalar EGGAMM;
  PetscScalar EGGAX;
  PetscScalar EGGAL;

  PetscScalar EGBETA;
  PetscScalar EGBEX;
  PetscScalar EGBEL;

  PetscScalar EGALPH;
  PetscScalar EGALX;
  PetscScalar EGALL;

  PetscScalar NC300;     // The effective density of states in the conduction band at 300K.
  PetscScalar NV300;     // The effective density of states in the valence band at 300K.
  PetscScalar NC_F;      // The parameter for temperature depended effective density of states in the conduction band.
  PetscScalar NV_F;      // The parameter for temperature depended effective density of states in the valence band.
  PetscScalar MEG;
  PetscScalar MEG_X1;
  PetscScalar MEX;
  PetscScalar MEX_X1;
  PetscScalar MEL;
  PetscScalar MEL_X1;
  PetscScalar MH0;
  PetscScalar MH0_X1;
  PetscScalar ML0;
  PetscScalar ML0_X1;
  // Model of Bandgap Narrowing due to Heavy Doping
  PetscScalar N0_BGN;    // The concentration parameter used in Slotboom's band-gap narrowing model.
  PetscScalar V0_BGN;    // The voltage parameter used in Slotboom's band-gap narrowing model.
  PetscScalar CON_BGN;   // The const parameter used in Slotboom's band-gap narrowing model.

  // Init value
  void Eg_Init()
  {
    EG300     = 3.540000E-01*eV;
    EG_X0     = 0.000000E+00*eV;
    EG_X1     = 2.865000E-01*eV;
    EG_X2     = 1.828000E+00*eV;
    EG_X3     = -2.948000E+00*eV;
    EG_X4     = 1.903500E+00*eV;
    EG_X5     = 1.786000E+00*eV;
    EG_X6     = -2.400000E-01*eV;
    EG_X7     = 0.000000E+00*eV;
    EG_X8     = 0.000000E+00*eV;
    EG_X9     = 0.000000E+00*eV;
    EG_X10    = 1.096000E+00*eV;
    EG_X11    = 2.600000E-01*eV;
    EG_X12    = 0.000000E+00*eV;
    EG_X13    = 0.000000E+00*eV;
    EG_X14    = 0.000000E+00*eV;
    pm        = 1.5;

    EGGAMM    = 1.545000E-04*eV/K;
    EGGAX     = 0.000000E+00*eV/K;
    EGGAL     = 0.000000E+00*eV/K;

    EGBETA    = -1.329900E+01*K;
    EGBEX     = 0.000000E+00*K;
    EGBEL     = 0.000000E+00*K;

    EGALPH    = 4.660000E-04*eV/K;
    EGALX     = 0.000000E+00*eV/K;
    EGALL     = 0.000000E+00*eV/K;

    NC300     = 8.130000E+16*pow(cm,-3);
    NV300     = 1.520000E+19*pow(cm,-3);
    NC_F      = 1.500000e+00;
    NV_F      = 1.500000e+00;

    MEG       = 2.190000E-02*me;
    MEG_X1    = 4.510000E-02*me;
    MEX       = 0.000000E+00*me;
    MEX_X1    = 8.500000E-01*me;
    MEL       = 0.000000E+00*me;
    MEL_X1    = 5.600000E-01*me;
    MH0       = 6.000000E-01*me;
    MH0_X1    = 2.000000E-02*me;
    ML0       = 2.700000E-02*me;
    ML0_X1    = 6.000000E-02*me;

    N0_BGN    = 1.000000e+17*pow(cm,-3);
    V0_BGN    = 0.000000e+00*V;
    CON_BGN   = 0.000000e+00*eV;
  }
public:
  //---------------------------------------------------------------------------
  // we need calculate all the bandgap valley and choose lowest
  PetscScalar E_Gamma(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X0 + EG_X1*x + EG_X2*x*x + EG_X3*x*x*x + EG_X4*x*x*x*x
               + (T300*T300/(T300+EGBETA)-Tl*Tl/(Tl+EGBETA))*(EGALPH+EGGAMM*x);
  }
  AutoDScalar E_Gamma(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X0 + EG_X1*x + EG_X2*x*x + EG_X3*x*x*x + EG_X4*x*x*x*x
               + (T300*T300/(T300+EGBETA)-Tl*Tl/(Tl+EGBETA))*(EGALPH+EGGAMM*x);
  }
  
  PetscScalar E_X(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X5 + EG_X6*x + EG_X7*x*x + EG_X8*x*x*x + EG_X9*x*x*x*x
               + (T300*T300/(T300+EGBEX)-Tl*Tl/(Tl+EGBEX))*(EGALX+EGGAX*x);
  }
  AutoDScalar E_X(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X5 + EG_X6*x + EG_X7*x*x + EG_X8*x*x*x + EG_X9*x*x*x*x
               + (T300*T300/(T300+EGBEX)-Tl*Tl/(Tl+EGBEX))*(EGALX+EGGAX*x);
  }
  
  PetscScalar E_L(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X10 + EG_X11*x + EG_X12*x*x + EG_X13*x*x*x + EG_X14*x*x*x*x
               + (T300*T300/(T300+EGBEL)-Tl*Tl/(Tl+EGBEL))*(EGALL+EGGAL*x);
  }
  AutoDScalar E_L(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return EG300 + EG_X10 + EG_X11*x + EG_X12*x*x + EG_X13*x*x*x + EG_X14*x*x*x*x
               + (T300*T300/(T300+EGBEL)-Tl*Tl/(Tl+EGBEL))*(EGALL+EGGAL*x);
  }
  
  //---------------------------------------------------------------------------
  // procedure of Bandgap, return the lowest valley
  PetscScalar Eg (const PetscScalar &Tl)
  {
    return E_Gamma(Tl);
  }
  AutoDScalar Eg (const AutoDScalar &Tl)
  {
    return E_Gamma(Tl);
  }

  //---------------------------------------------------------------------------
  // procedure of Bandgap Narrowing due to Heavy Doping
  PetscScalar EgNarrow(const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*pow(cm,-3);
    PetscScalar x = log(N/N0_BGN);
    return V0_BGN*(x+sqrt(x*x+CON_BGN));
  }
  PetscScalar EgNarrowToEc   (const PetscScalar &Tl){return 0.5*EgNarrow(Tl);}
  PetscScalar EgNarrowToEv   (const PetscScalar &Tl){return 0.5*EgNarrow(Tl);}
  
  AutoDScalar EgNarrow(const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    PetscScalar N = Na+Nd+1.0*pow(cm,-3);
    PetscScalar x = log(N/N0_BGN);
    return V0_BGN*(x+sqrt(x*x+CON_BGN));
  }
  AutoDScalar EgNarrowToEc   (const AutoDScalar &Tl){return 0.5*EgNarrow(Tl);}
  AutoDScalar EgNarrowToEv   (const AutoDScalar &Tl){return 0.5*EgNarrow(Tl);}
  
  //---------------------------------------------------------------------------
  //electron and hole effect mass
  PetscScalar EffecElecMass(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        PetscScalar bandgap = Eg(Tl);

        PetscScalar m_Gamma = pow(pow(MEG+MEG_X1*x,pm)*exp((bandgap-E_Gamma(Tl))/(kb*Tl)),1.0/pm);
        PetscScalar m_X = pow(pow(MEX+MEX_X1*x,pm)*exp((bandgap-E_X(Tl))/(kb*Tl)),1.0/pm);
        PetscScalar m_L = pow(pow(MEL+MEL_X1*x,pm)*exp((bandgap-E_L(Tl))/(kb*Tl)),1.0/pm);
        return pow(pow(m_Gamma,pm)+pow(m_X,pm)+pow(m_L,pm),1.0/pm);
  }
  AutoDScalar EffecElecMass(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        AutoDScalar bandgap = Eg(Tl);

        AutoDScalar m_Gamma = pow(pow(MEG+MEG_X1*x,pm)*exp((bandgap-E_Gamma(Tl))/(kb*Tl)),1.0/pm);
        AutoDScalar m_X = pow(pow(MEX+MEX_X1*x,pm)*exp((bandgap-E_X(Tl))/(kb*Tl)),1.0/pm);
        AutoDScalar m_L = pow(pow(MEL+MEL_X1*x,pm)*exp((bandgap-E_L(Tl))/(kb*Tl)),1.0/pm);
        return pow(pow(m_Gamma,pm)+pow(m_X,pm)+pow(m_L,pm),1.0/pm);
  }
  
  PetscScalar EffecHoleMass(const PetscScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return pow(pow(MH0+MH0_X1*x,pm)+pow(ML0+ML0_X1*x,pm),1.0/pm);
  }
  AutoDScalar EffecHoleMass(const AutoDScalar &Tl)
  {
        PetscScalar x = ReadxMoleFraction();
        return pow(pow(MH0+MH0_X1*x,pm)+pow(ML0+ML0_X1*x,pm),1.0/pm);
  }
  
  //---------------------------------------------------------------------------
  // Nc and Nv calculated from effective mass
  PetscScalar Nc (const PetscScalar &Tl)
  {
    //return NC300*pow(Tl/T300,NC_F);
    return 2*pow(2*3.14159265359*EffecElecMass(Tl)*kb*Tl/(h*h),pm);
  }
  AutoDScalar Nc (const AutoDScalar &Tl)
  {
    //return NC300*pow(Tl/T300,NC_F);
    return 2*pow(2*3.14159265359*EffecElecMass(Tl)*kb*Tl/(h*h),pm);
  }
  
  PetscScalar Nv (const PetscScalar &Tl)
  {
    //return NV300*pow(Tl/T300,NV_F);
    return 2*pow(2*3.14159265359*EffecHoleMass(Tl)*kb*Tl/(h*h),pm);
  }
  AutoDScalar Nv (const AutoDScalar &Tl)
  {
    //return NV300*pow(Tl/T300,NV_F);
    return 2*pow(2*3.14159265359*EffecHoleMass(Tl)*kb*Tl/(h*h),pm);
  }
  
  //---------------------------------------------------------------------------
  // nie
  PetscScalar nie (const PetscScalar &Tl)
  {
    PetscScalar bandgap = Eg(Tl);
    PetscScalar Nc = NC300*pow(Tl/T300,NC_F);
    PetscScalar Nv = NV300*pow(Tl/T300,NV_F);
    return sqrt(Nc*Nv)*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(Tl));
  }
  AutoDScalar nie (const AutoDScalar &Tl)
  {
    AutoDScalar bandgap = Eg(Tl);
    AutoDScalar Nc = NC300*pow(Tl/T300,NC_F);
    AutoDScalar Nv = NV300*pow(Tl/T300,NV_F);
    return sqrt(Nc*Nv)*exp(-bandgap/(2*kb*Tl))*exp(EgNarrow(Tl));
  }

  //end of Bandgap

private:
  //[Lifetime]
  //Lifetimes
  PetscScalar TAUN0;         // The Shockley-Read-Hall electron lifetime.
  PetscScalar TAUP0;         // The Shockley-Read-Hall hole lifetime.
  PetscScalar SurfTauN;      // The electron surface recombination velocity.
  PetscScalar SurfTauP;      // The hole surface recombination velocity.
  //Concentration-Dependent Lifetimes
  PetscScalar NSRHN;         // The Shockley-Read-Hall concentration parameter for electrons.
  PetscScalar AN;            // The constant term in the concentration-dependent expression for electron lifetime.
  PetscScalar BN;            // The linear term coefficient in the concentration-dependent expression for electron lifetime.
  PetscScalar CN;            // The exponential term coefficient in the concentration-dependent expression for electron lifetime.
  PetscScalar EN;            // The exponent in the concentration-dependent expression for electron lifetime.
  PetscScalar NSRHP;            // The Shockley-Read-Hall concentration parameter for holes.
  PetscScalar AP;            // The constant term in the concentration-dependent expression for hole lifetime.
  PetscScalar BP;            // The linear term coefficient in the concentration-dependent expression for hole lifetime.
  PetscScalar CP;            // The exponential term coefficient in the concentration-dependent expression for hole lifetime.
  PetscScalar EP;            // The exponent in the concentration-dependent expression for hole lifetime.
  // Lattice Temperature-Dependent Lifetimes
  PetscScalar EXN_TAU;       // The exponent of lattice temperature dependent electron lifetime.
  PetscScalar EXP_TAU;       // The exponent of lattice temperature dependent hole lifetime.

  //Init value
  void Lifetime_Init()
  {
    TAUN0     = 1.000000e-09*s;
    TAUP0     = 1.000000e-09*s;
    SurfTauN  = 0.000000e+00*cm/s;
    SurfTauP  = 0.000000e+00*cm/s;
    NSRHN     = 5.000000e+16*pow(cm,-3);
    AN        = 1.000000e+00;
    BN        = 0.000000e+00;
    CN        = 0.000000e+00;
    EN        = 2.000000e+00;
    NSRHP     = 5.000000e+16*pow(cm,-3);
    AP        = 1.000000e+00;
    BP        = 0.000000e+00;
    CP        = 0.000000e+00;
    EP        = 2.000000e+00;
    EXN_TAU   = 0.000000e+00;
    EXP_TAU   = 0.000000e+00;
  }

public:
  //---------------------------------------------------------------------------
  // electron lift time for SHR Recombination
  PetscScalar TAUN (const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUN0/(1+(Na+Nd)/NSRHN)*pow(Tl/T300,EXN_TAU);
  }
  AutoDScalar TAUN (const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUN0/(1+(Na+Nd)/NSRHN)*pow(Tl/T300,EXN_TAU);
  }
  
  //---------------------------------------------------------------------------
  // hole lift time for SHR Recombination
  PetscScalar TAUP (const PetscScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUP0/(1+(Na+Nd)/NSRHP)*pow(Tl/T300,EXP_TAU);
  }
  AutoDScalar TAUP (const AutoDScalar &Tl)
  {
    PetscScalar Na = ReadDopingNa();
    PetscScalar Nd = ReadDopingNd();
    return TAUP0/(1+(Na+Nd)/NSRHP)*pow(Tl/T300,EXP_TAU);
  }
  // End of Lifetime

private:
  //[Recombination]
  // SRH, Auger, and Direct Recombination
  PetscScalar ETRAP;         // The trap level (Et - Ei) used in determining the Shockley-Read-Hall recombination rate.
  PetscScalar AUGN;          // The Auger coefficient for electrons.
  PetscScalar AUGP;          // The Auger coefficient for holes.
  PetscScalar C_DIRECT;      // The band-to-band recombination coefficient.
  // Recombination Including Tunneling
  PetscScalar M_RTUN;        // The trap-assisted tunneling effective mass. *free electron rest mass m0
  PetscScalar S_RTUN;        // Band-to-band field power ratio.
  PetscScalar B_RTUN;        // Band-to-band tunneling rate proportionality factor.
  PetscScalar E_RTUN;        // Band-to-band reference electric field.

  // Init value
  void Recomb_Init()
  {
    //Source: Semiconductors on NSM
    ETRAP   =  0.000000e+00*eV;
    AUGN    =  7.000000e-29*pow(cm,6)/s;
    AUGP    =  7.000000e-29*pow(cm,6)/s;
    C_DIRECT = 0.960000e-10*pow(cm,3)/s;
    M_RTUN   = 2.500000e-01;
    S_RTUN   = 0.000000e+00;
    B_RTUN   = 0.000000e+00*pow(cm,S_RTUN -3)*pow(V,S_RTUN*-1)/s;
    E_RTUN   = 0.000000e+00*V/cm;
  }

public:
  //---------------------------------------------------------------------------
  // Direct Recombination
  PetscScalar R_Direct     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(Tl);
    return C_DIRECT*(n*p-ni*ni);
  }
  AutoDScalar R_Direct     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(Tl);
    return C_DIRECT*(n*p-ni*ni);
  }
 
  //---------------------------------------------------------------------------
  // Auger Recombination
  PetscScalar R_Auger     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  { 
    PetscScalar ni =   nie(Tl);
    return AUGN*(p*n*n-n*ni*ni)+AUGP*(n*p*p-p*ni*ni);
  }  
  AutoDScalar R_Auger     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  { 
    AutoDScalar ni =   nie(Tl);
    return AUGN*(p*n*n-n*ni*ni)+AUGP*(n*p*p-p*ni*ni);
  } 

  //---------------------------------------------------------------------------
  // SHR Recombination  
  PetscScalar R_SHR     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(Tl);
    PetscScalar taun = TAUN(Tl);
    PetscScalar taup = TAUP(Tl);
    return (p*n-ni*ni)/(taup*(n+ni)+taun*(p+ni));
  }
  AutoDScalar R_SHR     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(Tl);
    AutoDScalar taun = TAUN(Tl);
    AutoDScalar taup = TAUP(Tl);
    return (p*n-ni*ni)/(taup*(n+ni)+taun*(p+ni));
  }

  //---------------------------------------------------------------------------
  // Surface SHR Recombination  
  PetscScalar R_Surf     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, const PetscScalar &reciprocal_len)
  {
    PetscScalar ni =   nie(Tl);
    PetscScalar taun = TAUN(Tl);
    PetscScalar taup = TAUP(Tl);
    taun = 1.0/(SurfTauN*reciprocal_len + 1.0/taun);
    taup = 1.0/(SurfTauP*reciprocal_len + 1.0/taup);
    return (p*n-ni*ni)/(taup*(n+ni)+taun*(p+ni));
  }
  AutoDScalar R_Surf     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, const PetscScalar &reciprocal_len)
  {
    AutoDScalar ni =   nie(Tl);
    AutoDScalar taun = TAUN(Tl);
    AutoDScalar taup = TAUP(Tl);
    taun = 1.0/(SurfTauN*reciprocal_len + 1.0/taun);
    taup = 1.0/(SurfTauP*reciprocal_len + 1.0/taup);
    return (p*n-ni*ni)/(taup*(n+ni)+taun*(p+ni));
  }

  //---------------------------------------------------------------------------
  // total Recombination  
  PetscScalar Recomb (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl)
  {
    PetscScalar ni =   nie(Tl);
    PetscScalar taun = TAUN(Tl);
    PetscScalar taup = TAUP(Tl);
    PetscScalar dn   = p*n-ni*ni; 
    PetscScalar Rshr = dn/(taup*(n+ni)+taun*(p+ni));
    PetscScalar Rdir = C_DIRECT*dn;
    PetscScalar Raug = (AUGN*n+AUGP*p)*dn;
    return Rshr+Rdir+Raug;
  }
  AutoDScalar Recomb (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl)
  {
    AutoDScalar ni =   nie(Tl);
    AutoDScalar taun = TAUN(Tl);
    AutoDScalar taup = TAUP(Tl);
    AutoDScalar dn   = p*n-ni*ni; 
    AutoDScalar Rshr = dn/(taup*(n+ni)+taun*(p+ni));
    AutoDScalar Rdir = C_DIRECT*dn;
    AutoDScalar Raug = (AUGN*n+AUGP*p)*dn;
    return Rshr+Rdir+Raug;
  }

  // End of Recombination
  
private:
  //[energy relax time]
  PetscScalar  WTN0;
  PetscScalar  WTN1;
  PetscScalar  WTN2;
  PetscScalar  WTN3;
  PetscScalar  WTN4;
  PetscScalar  WTN5;
  PetscScalar  WTNL;
  PetscScalar  TNL;
  PetscScalar  WTP0;
  PetscScalar  WTP1;
  PetscScalar  WTP2;
  PetscScalar  WTP3;
  PetscScalar  WTP4;
  PetscScalar  WTP5;
  PetscScalar  WTPL;
  PetscScalar  TPL;
  // Init value
  void RelaxTime_Init()
  {
   WTN0 =  2.400000E-12*s;         
   WTN1 =  4.000000E-13*s;            
   WTN2 =  0.000000E+00*s;        
   WTN3 =  0.000000E+00*s;          
   WTN4 =  0.000000E+00*s;          
   WTN5 =  0.000000E+00*s;            
   WTNL =  6.800000E-13*s;          
   TNL  =  1.866270E+03*K;         
   WTP0 =  0.000000E+00*s;       
   WTP1 =  0.000000E+00*s;            
   WTP2 =  0.000000E+00*s;          
   WTP3 =  0.000000E+00*s;        
   WTP4 =  0.000000E+00*s;       
   WTP5 =  0.000000E+00*s;            
   WTPL =  1.000000E-12*s;          
   TPL  =  0.000000E+00*K;
  }
public:
  //---------------------------------------------------------------------------
  // Electron relaxation time for EBM 
  PetscScalar ElecEnergyRelaxTime(const PetscScalar &Tn,const PetscScalar &Tl)
  {
    PetscScalar r = (Tn-Tl)/TNL;
    return WTN1+(WTN0-WTN1)*r*r*exp(2-2*r);
  }
  AutoDScalar ElecEnergyRelaxTime(const AutoDScalar &Tn,const AutoDScalar &Tl)
  {
    AutoDScalar r = (Tn-Tl)/TNL;
    return WTN1+(WTN0-WTN1)*r*r*exp(2-2*r);
  }

  //---------------------------------------------------------------------------
  // Hole relaxation time for EBM 
  PetscScalar HoleEnergyRelaxTime(const PetscScalar &Tp,const PetscScalar &Tl)
  {
    return WTPL;
  }
  AutoDScalar HoleEnergyRelaxTime(const AutoDScalar &Tp,const AutoDScalar &Tl)
  {
    return WTPL;
  }
  // end of energy relax time
  
private:
  // [Schottky]
  PetscScalar ARICHN;
  PetscScalar ARICHP;
  PetscScalar VSURFN;   // Thermionic emission velocity of electron
  PetscScalar VSURFP;

  void   Schottky_Init()
  {
    ARICHN = 2.631800E+00*A/(K*cm)/(K*cm);
    ARICHP = 7.258600E+01*A/(K*cm)/(K*cm);
  }
public:
  PetscScalar SchottyJsn (PetscScalar n,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFN = ARICHN*Tl*Tl/(e*Nc(Tl));
    PetscScalar nb = Nc(Tl)*exp(-e*Vb/(kb*Tl));
    return -e*VSURFN*(n-nb);
  }
  PetscScalar SchottyJsp (PetscScalar p,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFP = ARICHP*Tl*Tl/(e*Nv(Tl));
    PetscScalar pb = Nv(Tl)*exp((-Eg(Tl)+e*Vb)/(kb*Tl));
    return e*VSURFP*(p-pb);
  }
  PetscScalar SchottyBarrierLowerring (PetscScalar eps, PetscScalar E)
  {
    return sqrt(e/(4*3.1415926535*eps)*E);
  }
  PetscScalar pdSchottyJsn_pdn(PetscScalar n,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFN = ARICHN*Tl*Tl/(e*Nc(Tl));
    return -e*VSURFN;
  }
  PetscScalar pdSchottyJsp_pdp(PetscScalar p,PetscScalar Tl,PetscScalar Vb)
  {
    PetscScalar VSURFP = ARICHP*Tl*Tl/(e*Nv(Tl));
    return e*VSURFP;
  }
  PetscScalar pdSchottyJsn_pdTl(PetscScalar n,PetscScalar Tl,PetscScalar Vb)
  {
    //use finite difference approximate
    PetscScalar dJ = SchottyJsn(n,Tl,Vb)-SchottyJsn(n,(1-1e-10)*Tl,Vb);
    return dJ/(1e-10*Tl);
  }
  PetscScalar pdSchottyJsp_pdTl(PetscScalar p,PetscScalar Tl,PetscScalar Vb)
  {
    //use finite difference approximate
    PetscScalar dJ = SchottyJsp(p,Tl,Vb)-SchottyJsp(p,(1-1e-10)*Tl,Vb);
    return dJ/(1e-10*Tl);
  }

  PetscScalar ThermalVn (PetscScalar Tl)
  {
        return sqrt(kb*Tl/(2*3.14159265359*EffecElecMass(Tl)));
  }
  PetscScalar ThermalVp (PetscScalar Tl)
  {
        return sqrt(kb*Tl/(2*3.14159265359*EffecHoleMass(Tl)));
  }
  PetscScalar pdThermalVn_pdTl (PetscScalar Tl)
  {
        return 0;
  }
  PetscScalar pdThermalVp_pdTl (PetscScalar Tl)
  {
        return 0;
  }

private:
  // [band to band Tunneling]
  PetscScalar  A_BTBT;
  PetscScalar  B_BTBT;
  void   BBTunneling_Init()
  {
    A_BTBT = 0*e*sqrt(V)/cm/s/V/V;
    B_BTBT = 0*V/cm/pow(e*V,PetscScalar(1.5));
  }
public:
  //----------------------------------------------------------------
  // band to band Tunneling
  PetscScalar BB_Tunneling(const PetscScalar &Tl,const  PetscScalar &E)
  {
     return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }
  AutoDScalar BB_Tunneling(const AutoDScalar &Tl,const  AutoDScalar &E)
  {
     return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }
 

// constructor and destructor
public:
  GSS_InGaAs_BandStructure(const PMIS_Environment &env):PMIS_BandStructure(env)
  {
    T300 = 300.0*K;
    Eg_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    BBTunneling_Init();
  }

  ~GSS_InGaAs_BandStructure()
  {}
}
;


extern "C"
{
  PMIS_BandStructure*  PMIS_InGaAs_BandStructure_Default (const PMIS_Environment& env)
  {
    return new GSS_InGaAs_BandStructure(env);
  }
}
