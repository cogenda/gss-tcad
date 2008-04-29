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
// Material Type: Cu


#include "PMI.h"

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;

//Source: FreeSnell
static table  Cu[] = {
{0.0200,  0.97,    0.11},
{0.0400,  0.8799,  0.2401},
{0.0600,  0.8935,  0.423},
{0.0800,  0.9801,  0.6901},
{0.0900,  1.064,   0.72},
{0.1000,  1.086,   0.726},	
{0.1500,  1.03,    1.00},
{0.2000,  0.988,   1.504},
{0.3000,  1.393,   1.670},
{0.4000,  1.175,   2.130},
{0.5000,  1.134,   2.570 },
{0.6000,  0.747,   3.362 },
{0.7000,  0.412,   4.203 },
{0.8000,  0.454,   4.978 },
{0.9000,  0.496,   5.754 },
{1.0000,  0.538,   6.530 },
{1.5000,  0.733,   10.017},    
{2.0000,  0.879,   13.400},
{3.0000,  1.562,   20.055},
{4.0000,  2.250,   26.500},
{5.0000,  3.260,   33.000},
{6.0000,  3.944,   39.078},
{7.0000,  4.942,   45.336},
{8.0000,  6.074,   51.292},
{9.0000,  7.118,   57.167},
{10.000,  8.310,   63.000},
};        

class GSS_Cu_Optical : public PMIC_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0.0, PetscScalar Tl=1.0) const
  {
    complex<PetscScalar> n(1.0,0.0);
    if(lamda<Cu[0].wavelength) 
      return complex<PetscScalar> (Cu[0].RefractionIndexRe,Cu[0].RefractionIndexIm);
    if(lamda>Cu[table_size-1].wavelength) 
      return complex<PetscScalar> (Cu[table_size-1].RefractionIndexRe,Cu[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=Cu[i].wavelength && lamda<=Cu[i+1].wavelength)
       {
          complex<PetscScalar> n1(Cu[i].RefractionIndexRe,Cu[i].RefractionIndexIm);
          complex<PetscScalar> n2(Cu[i+1].RefractionIndexRe,Cu[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-Cu[i].wavelength;
          PetscScalar d2 = Cu[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2);
          break;  
       }
    }
    return n;   
  }                                            
  
  // constructions
public:
  GSS_Cu_Optical(const PMIC_Environment &env):PMIC_Optical(env)
  {
     table_size = sizeof(Cu)/sizeof(table);
     //the wave length should be scaled
     for(int i=0;i<table_size;i++)
     {
        Cu[i].wavelength*=um;    
     }
  }

  ~GSS_Cu_Optical(){}
}
;

extern "C"
{
  PMIC_Optical*  PMIC_Cu_Optical_Default (const PMIC_Environment& env)
  {
    return new GSS_Cu_Optical(env);
  }
}
