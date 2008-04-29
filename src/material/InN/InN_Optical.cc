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
/*       A Two-Dimensional General Purpose Semiconductor GaAsmulator.        */
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
// Material Type: InN


#include "PMI.h"

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;


static table  InN[] =  {
//i can't find any reflection data for InN... 
};

class GSS_InN_Optical : public PMIS_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=1.12, PetscScalar Tl=1.0) const
  {
    complex<PetscScalar> n(1.0,0.0);
    if(!table_size) return n;
    if(lamda<InN[0].wavelength) 
      return complex<PetscScalar> (InN[0].RefractionIndexRe,InN[0].RefractionIndexIm);
    if(lamda>InN[table_size-1].wavelength) 
      return complex<PetscScalar> (InN[table_size-1].RefractionIndexRe,InN[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=InN[i].wavelength && lamda<=InN[i+1].wavelength)
       {
          complex<PetscScalar> n1(InN[i].RefractionIndexRe,InN[i].RefractionIndexIm);
          complex<PetscScalar> n2(InN[i+1].RefractionIndexRe,InN[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-InN[i].wavelength;
          PetscScalar d2 = InN[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2);
          break;  
       }
    }
    return n;   
  }           

// constructions
public:
  GSS_InN_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    table_size = sizeof(InN)/sizeof(table);
    //the wave length should be scaled
    for(int i=0;i<table_size;i++)
    {
        InN[i].wavelength*=um;    
    }
  }

  ~GSS_InN_Optical()
  {}
}
;

extern "C"
{
  PMIS_Optical*  PMIS_InN_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_InN_Optical(env);
  }
}

