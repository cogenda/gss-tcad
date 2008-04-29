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
// Material Type: 3C-SiC. 


#include "PMI.h"

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;

//Source: FreeSnell
static table  SiC3C[] = {
{0.1500,  4.021 , 1.720},
{0.2000,  3.960 , 1.06},
{0.2100,  3.759 , 0.702},
{0.2200,  3.565 , 0.487},
{0.2400,  3.252 , 0.307},
{0.2600,  3.044 , 0.206},
{0.2800,  2.990 , 0.195},
{0.3000,  2.953 , 0.18},
{0.3300,  2.896 , 0.156},
{0.3500,  2.844 , 0.119},
{0.4000,  2.764 , 0.1},
{0.4200,  2.744 , 0.0443},
{0.4500,  2.720 , 7.59e-4},
{0.5000,  2.686 , 0.000},
{0.6000,  2.644 , 0.000},
{0.7000,  2.620 , 0.000},
{0.8000,  2.601 , 0.000},
{0.9000,  2.589 , 0.000},
{1.0000 , 2.573,  0.000},
};

class GSS_SiC3C_Optical : public PMIS_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=1.12, PetscScalar Tl=1.0) const
  {
    complex<PetscScalar> n(1.0,0.0);
    if(lamda<SiC3C[0].wavelength) 
      return complex<PetscScalar> (SiC3C[0].RefractionIndexRe,SiC3C[0].RefractionIndexIm);
    if(lamda>SiC3C[table_size-1].wavelength) 
      return complex<PetscScalar> (SiC3C[table_size-1].RefractionIndexRe,SiC3C[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=SiC3C[i].wavelength && lamda<=SiC3C[i+1].wavelength)
       {
          complex<PetscScalar> n1(SiC3C[i].RefractionIndexRe,SiC3C[i].RefractionIndexIm);
          complex<PetscScalar> n2(SiC3C[i+1].RefractionIndexRe,SiC3C[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-SiC3C[i].wavelength;
          PetscScalar d2 = SiC3C[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2); 
	  break; 
       }
    }
    return n;   
  }           

// constructions
public:
  GSS_SiC3C_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    table_size = sizeof(SiC3C)/sizeof(table);
    //the wave length should be scaled
    for(int i=0;i<table_size;i++)
    {
        SiC3C[i].wavelength*=um;    
    }
  } 
    
  ~GSS_SiC3C_Optical()
  {}
}   
;

extern "C"
{
  PMIS_Optical*  PMIS_SiC3C_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_SiC3C_Optical(env);
  }
}
