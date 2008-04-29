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
// Material Type: Hg(1-x)Cd(x)Te

#include "PMI.h"

//sorry, no impact ionization data can be found at present. 
class GSS_HgCdTe_Avalanche : public PMIS_Avalanche
{
public:
  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for DDM 
  PetscScalar ElecGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    return 0;
  }
  AutoDScalar ElecGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    return 0;
  }
  
  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for DDM 
  PetscScalar HoleGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const
  {
    return 0;
  }
  AutoDScalar HoleGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const
  {
    return 0;
  }
  
 
  
  //---------------------------------------------------------------------------
  // Electron Impact Ionization rate for EBM
  PetscScalar ElecGenRateEBM (const PetscScalar &Tn,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    return 0;

  }
  AutoDScalar ElecGenRateEBM (const AutoDScalar &Tn,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    return 0;
  }
  
  //---------------------------------------------------------------------------
  // Hole Impact Ionization rate for EBM
  PetscScalar HoleGenRateEBM (const PetscScalar &Tp,const PetscScalar &Tl,const PetscScalar &Eg) const
  {
    return 0;
  }
  AutoDScalar HoleGenRateEBM (const AutoDScalar &Tp,const AutoDScalar &Tl,const AutoDScalar &Eg) const
  {
    return 0;
  }
  

//----------------------------------------------------------------
// constructor and destructor
public:  
  GSS_HgCdTe_Avalanche(const PMIS_Environment &env):PMIS_Avalanche(env)
  {
  }
  ~GSS_HgCdTe_Avalanche()
  {
  }

}
;

extern "C"
{
  PMIS_Avalanche* PMIS_HgCdTe_Avalanche_Default (const PMIS_Environment& env)
  {
    return new GSS_HgCdTe_Avalanche(env);
  }
}
