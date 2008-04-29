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
// Material Type: Ge


#include "PMI.h"


class GSS_Ge_BandStructure : public PMIS_BandStructure
{
private:
  PetscScalar T300;
  //[Bandgap]
  // Bandgap and Effective Density of States
  PetscScalar EG300;     // The energy bandgap of the material at 300 K.
  PetscScalar EGALPH;    // The value of alpha used in calculating the temperature depended energy bandgap.
  PetscScalar EGBETA;    // The value of beta  used in calculating the temperature depended energy bandgap.
  PetscScalar ELECMASS;  // The relative effective mass of electron
  PetscScalar HOLEMASS;  // The relative effective mass of hole
  PetscScalar NC300;     // The effective density of states in the conduction band at 300K.
  PetscScalar NV300;     // The effective density of states in the valence band at 300K.
  PetscScalar NC_F;      // The parameter for temperature depended effective density of states in the conduction band.
  PetscScalar NV_F;      // The parameter for temperature depended effective density of states in the valence band.
  // Model of Bandgap Narrowing due to Heavy Doping
  PetscScalar N0_BGN;    // The concentration parameter used in Slotboom's band-gap narrowing model.
  PetscScalar V0_BGN;    // The voltage parameter used in Slotboom's band-gap narrowing model.
  PetscScalar CON_BGN;   // The const parameter used in Slotboom's band-gap narrowing model.

  // Init value
  void Eg_Init()
  {
    EG300     = 6.600000e-01*eV;
    EGALPH    = 4.774000e-04*eV/K;
    EGBETA    = 2.350000e+02*K;
    ELECMASS  = 3.100000E-01*me;
    HOLEMASS  = 5.000000E-01*me;
    NC300     = 1.040000e+19*pow(cm,-3);
    NV300     = 6.000000e+18*pow(cm,-3);
    NC_F      = 1.500000e+00;
    NV_F      = 1.500000e+00;
    N0_BGN    = 1.000000e+17*pow(cm,-3);
    V0_BGN    = 9.000000e-03*V;
    CON_BGN   = 5.000000e-01*eV;
  }
public:
  // function
  PetscScalar Eg (const PetscScalar &Tl)
  {
    return EG300+EGALPH*(T300*T300/(T300+EGBETA) - Tl*Tl/(Tl+EGBETA));
  }
  AutoDScalar Eg (const AutoDScalar &Tl)
  {
    return EG300+EGALPH*(T300*T300/(T300+EGBETA) - Tl*Tl/(Tl+EGBETA));
  }


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
  

  PetscScalar EffecElecMass (const PetscScalar &Tl)
  {
        return ELECMASS;
  }
  PetscScalar EffecHoleMass (const PetscScalar &Tl)
  {
        return HOLEMASS;
  }
  // interface procedure
  PetscScalar Nc (const PetscScalar &Tl)
  {
    return NC300*pow(Tl/T300,NC_F);
  }
  AutoDScalar Nc (const AutoDScalar &Tl)
  {
    return NC300*pow(Tl/T300,NC_F);
  }
  PetscScalar Nv (const PetscScalar &Tl)
  {
    return NV300*pow(Tl/T300,NV_F);
  }
  AutoDScalar Nv (const AutoDScalar &Tl)
  {
    return NV300*pow(Tl/T300,NV_F);
  }
  
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
    TAUN0     = 1.000000e-07*s;
    NSRHN     = 5.000000e+16*pow(cm,-3);
    AN        = 1.000000e+00;
    BN        = 1.000000e+00;
    CN        = 0.000000e+00;
    EN        = 2.000000e+00;
    TAUP0     = 1.000000e-07*s;
    NSRHP     = 5.000000e+16*pow(cm,-3);
    AP        = 1.000000e+00;
    BP        = 1.000000e+00;
    CP        = 0.000000e+00;
    EP        = 2.000000e+00;
    EXN_TAU   = 0.000000e+00;
    EXP_TAU   = 0.000000e+00;
    SurfTauN  = 0.000000e+00*cm/s;
    SurfTauP  = 0.000000e+00*cm/s;
  }

public:
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
  
  //the fit parameter for density-gradient solver
  PetscScalar Gamman         () {return 3.6;}
  PetscScalar Gammap         () {return 5.6;}
  
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
    ETRAP   = 0.000000e+00*eV;
    AUGN    =  2.800000e-31*pow(cm,6)/s;
    AUGP    =  9.900000e-32*pow(cm,6)/s;
    C_DIRECT = 0.000000e+00*pow(cm,3)/s;
    M_RTUN   = 2.500000e-01;
    S_RTUN   = 2.500000e+00;
    B_RTUN   = 4.000000e+14*pow(cm,S_RTUN -3)*pow(V,S_RTUN*-1)/s;
    E_RTUN   = 1.900000e+07*V/cm;
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
    PetscScalar Rshr = (p*n-ni*ni)/(taup*(n+ni)+taun*(p+ni));
    PetscScalar Rdir = C_DIRECT*(n*p-ni*ni);
    PetscScalar Raug = AUGN*(p*n*n-n*ni*ni)+AUGP*(n*p*p-p*ni*ni);
    if(Raug<0) Raug = 0.0;
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
   WTN0 =  1.685200E-13*s;         
   WTN1 =  1.029900E-13*s;            
   WTN2 = -5.184500E-15*s;        
   WTN3 =  0.000000E+00*s;          
   WTN4 =  0.000000E+00*s;          
   WTN5 =  0.000000E+00*s;            
   WTNL =  2.000000E-13*s;          
   TNL  =  2.979800E+03*K;         
   WTP0 = -1.560000E-14*s;         
   WTP1 =  1.380000E-13*s;            
   WTP2 = -2.500000E-14*s;         
   WTP3 =  2.310000E-15*s;          
   WTP4 = -1.050000E-16*s;         
   WTP5 =  1.820000E-18*s;            
   WTPL =  2.000000E-13*s;          
   TPL  =  1.000000E+05*K;
  }
public:
  //---------------------------------------------------------------------------
  // Electron relaxation time for EBM
  PetscScalar ElecEnergyRelaxTime(const PetscScalar &Tn,const PetscScalar &Tl)
  {
    if(Tn>TNL)     return WTNL;
    PetscScalar x = 1+(Tn-Tl)/T300;
    return WTN0+ WTN1*x + WTN2*x*x;
  }
  AutoDScalar ElecEnergyRelaxTime(const AutoDScalar &Tn,const AutoDScalar &Tl)
  {
    if(Tn>TNL)     return WTNL;
    AutoDScalar x = 1+(Tn-Tl)/T300;
    return WTN0+ WTN1*x + WTN2*x*x;
  }

  
  //---------------------------------------------------------------------------
  // Hole relaxation time for EBM 
  PetscScalar HoleEnergyRelaxTime(const PetscScalar &Tp,const PetscScalar &Tl)
  {
    if(Tp>TPL)     return WTPL;
    PetscScalar x = 1+(Tp-Tl)/T300;
    return WTP0+ WTP1*x + WTP2*x*x + WTP3*x*x*x + WTP4*pow(x,4) + WTP5*pow(x,5);
  }
  AutoDScalar HoleEnergyRelaxTime(const AutoDScalar &Tp,const AutoDScalar &Tl)
  {
    if(Tp>TPL)     return WTPL;
    AutoDScalar x = 1+(Tp-Tl)/T300;
    return WTP0+ WTP1*x + WTP2*x*x + WTP3*x*x*x + WTP4*pow(x,4) + WTP5*pow(x,5);
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
    ARICHN = 1.100000e+02*A/(K*cm)/(K*cm);
    ARICHP = 3.000000e+01*A/(K*cm)/(K*cm);
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
    A_BTBT = 3.500000E+21*e*sqrt(V)/cm/s/V/V;
    B_BTBT = 2.250000E+07*V/cm/pow(e*V,PetscScalar(1.5));
  }
public:
  //----------------------------------------------------------------
  // band to band Tunneling
  PetscScalar BB_Tunneling(const PetscScalar &Tl, const PetscScalar &E)
  {
     return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }
  AutoDScalar BB_Tunneling(const AutoDScalar &Tl, const AutoDScalar &E)
  {
     return A_BTBT*E*E/sqrt(Eg(Tl))*exp(-B_BTBT*pow(Eg(Tl),PetscScalar(1.5))/(E+1*V/cm));
  }


// constructor and destructor
public:
  GSS_Ge_BandStructure(const PMIS_Environment &env):PMIS_BandStructure(env)
  {
    T300 = 300.0*K;
    Eg_Init();
    Lifetime_Init();
    Recomb_Init();
    RelaxTime_Init();
    Schottky_Init();
    BBTunneling_Init();
  }

  ~GSS_Ge_BandStructure()
  {}
}
;


extern "C"
{
  PMIS_BandStructure*  PMIS_Ge_BandStructure_Default (const PMIS_Environment& env)
  {
    return new GSS_Ge_BandStructure(env);
  }
}
