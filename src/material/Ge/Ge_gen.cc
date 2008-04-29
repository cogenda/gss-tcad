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

class GSS_Ge_Avalanche : public PMIS_Avalanche
{
private:
  PetscScalar N_IONIZA; // The constant term in the multiplicative prefactor of the electron ionization coefficient.
  PetscScalar ECN_II;   // The critical electric field used in the exponential factor of the electron ionization coefficient.
  PetscScalar EXN_II;   // The exponent of the ratio of the critical electrical field to the local electric field.
  PetscScalar P_IONIZA; // The constant term in the multiplicative prefactor of the hole ionization coefficient.
  PetscScalar ECP_II;   // The critical electric field used in the exponential factor of the hole ionization coefficient.
  PetscScalar EXP_II;   // The exponent of the ratio of the critical electrical field to the local electric field.
  PetscScalar N_ION_1;  // The coefficient multiplying T in the multiplicative prefactor of the electron ionization coefficient.
  PetscScalar N_ION_2;  // The coefficient multiplying T^2 in the multiplicative prefactor of the electron ionization coefficient.
  PetscScalar P_ION_1;  // The coefficient multiplying T in the multiplicative prefactor of the hole ionization coefficient.
  PetscScalar P_ION_2;  // The coefficient multiplying T^2 in the multiplicative prefactor of the hole ionization coefficient.
  //// Impact Ionization Model Depending on Lattice Temperature.
  PetscScalar LAN300;   // Energy free path for electrons at 300 K, used for the impact ionization model depending on lattice temperature.
  PetscScalar LAP300;   // Energy free path for holes at 300 K, used for the impact ionization model depending on lattice temperature.
  PetscScalar OP_PH_EN; // Mean optical phonon energy used for the impact ionization model depending on lattice temperature.
  //
  PetscScalar ElecTauw;
  PetscScalar HoleTauw;
  PetscScalar T300    ;
  
  void 	Avalanche_Init()
  {
    N_IONIZA  =  1.550000e+07/cm;
    ECN_II    =  0.000000e+00*V/cm;
    EXN_II    =  1.000000e+00;
    P_IONIZA  =  1.000000e+07/cm;
    ECP_II    =  0.000000e+00*V/cm;
    EXP_II    =  1.000000e+00;
    N_ION_1   =  0.000000e+00/cm*K;
    N_ION_2   =  0.000000e+00/cm*K*K;
    P_ION_1   =  0.000000e+00/cm*K;
    P_ION_2   =  0.000000e+00/cm*K*K;
    LAN300    =  6.888250e-07*cm;
    LAP300    =  8.395050e-07*cm;
    OP_PH_EN  =  3.700000e-02*eV;
    
    ElecTauw  = 2e-13*s;
    HoleTauw  = 2e-13*s;
    T300      = 300.0*K;
  }
public:
  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for DDM 
  PetscScalar ElecGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    if (Ep < 1e4*V/cm)
    {
      return 0;
    }
    else
    {
      PetscScalar alpha = N_IONIZA + N_ION_1*Tl + N_ION_2*Tl*Tl;
      PetscScalar L=LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      return alpha*exp(-pow(Eg/(e*L)/Ep,EXN_II));
    }
  }
  AutoDScalar ElecGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    if (Ep < 1e4*V/cm)
    {
      return 0;
    }
    else
    {
      AutoDScalar alpha = N_IONIZA + N_ION_1*Tl + N_ION_2*Tl*Tl;
      AutoDScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      return alpha*exp(-pow(Ecrit/Ep,EXN_II));
    }
  }
  
  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for DDM
  PetscScalar HoleGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    if (Ep < 1e4*V/cm)
    {
      return 0;
    }
    else
    {
      PetscScalar alpha = P_IONIZA+P_ION_1*Tl+P_ION_2*Tl*Tl;
      PetscScalar L = LAP300*tanh(OP_PH_EN/(2*kb*Tl));
      return alpha*exp(-pow(Eg/(e*L)/Ep,EXP_II));
    }
  }
  AutoDScalar HoleGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    if (Ep < 1e4*V/cm)
    {
      return 0;
    }
    else
    {
      AutoDScalar alpha = P_IONIZA+P_ION_1*Tl+P_ION_2*Tl*Tl;
      AutoDScalar L = LAP300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      return alpha*exp(-pow(Ecrit/Ep,EXP_II));
    }
  }
  
  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for EBM
  PetscScalar ElecGenRateEBM (const PetscScalar &Tn,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    if (fabs(Tn - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
      PetscScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      PetscScalar Ecrit = Eg/(e*L);
      PetscScalar uc   = 2*vsat*ElecTauw/3*Ecrit;
      PetscScalar ut   = kb/e*(Tn-Tl);
      return N_IONIZA/e*exp(-pow(uc/ut,EXN_II));
    }
  }
  AutoDScalar ElecGenRateEBM (const AutoDScalar &Tn,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    if (fabs(Tn - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
      AutoDScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      AutoDScalar uc   = 2*vsat*ElecTauw/3*Ecrit;
      AutoDScalar ut   = kb/e*(Tn-Tl);
      return N_IONIZA/e*exp(-pow(uc/ut,EXN_II));
    }
  }
  
  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for EBM
  PetscScalar HoleGenRateEBM (const PetscScalar &Tp,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    if (fabs(Tp - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      PetscScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
      PetscScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      PetscScalar Ecrit = Eg/(e*L);
      PetscScalar uc   = 2*vsat*HoleTauw/3*Ecrit;
      PetscScalar ut   = kb/e*(Tp-Tl);
      return P_IONIZA/e*exp(-pow(uc/ut,EXP_II));
    }
  }
  AutoDScalar HoleGenRateEBM (const AutoDScalar &Tp,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    if (fabs(Tp - Tl)<1*K)
    {
      return 0;
    }
    else
    {
      AutoDScalar vsat = (2.4e7*cm/s)/(1+0.8*exp(Tl/(2*T300)));
      AutoDScalar L = LAN300*tanh(OP_PH_EN/(2*kb*Tl));
      AutoDScalar Ecrit = Eg/(e*L);
      AutoDScalar uc   = 2*vsat*HoleTauw/3*Ecrit;
      AutoDScalar ut   = kb/e*(Tp-Tl);
      return P_IONIZA/e*exp(-pow(uc/ut,EXP_II));
    }
  }

  
//----------------------------------------------------------------
// constructor and destructor
public:   
  GSS_Ge_Avalanche(const PMIS_Environment &env):PMIS_Avalanche(env)
  {
    Avalanche_Init();
  }
  ~GSS_Ge_Avalanche()
  {
  }

}
;

extern "C"
{
  PMIS_Avalanche* PMIS_Ge_Avalanche_Default (const PMIS_Environment& env)
  {
    return new GSS_Ge_Avalanche(env);
  }
}
