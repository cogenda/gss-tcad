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
// Material Type: Metal (Al)


#include "PMI.h"

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;


static table  Al[] = {
{0.0165,  1.010,  0.024},
{0.0172,  1.020,  0.0035},
{0.0376,  0.902,  0.0102},
{0.0620,  0.668,  0.0268},
{0.0800,  0.258,  0.0777},
{0.1033,  0.0328, 0.791},
{0.1600,  0.0765, 1.740},
{0.2066,  0.130,  2.39},
{0.3000,  0.276,  3.61},
{0.4000,  0.490,  4.86},
{0.5000,  0.769,  6.08},
{0.6000,  1.20,   7.26},
{0.7000,  1.83,   8.31},
{0.8000,  2.80,   8.45},
{0.9000,  2.06,   8.30},
{1.0000,  1.35,   9.58},
{1.2000,  1.21,   12.0},
{1.5000,  1.38,   15.4},
{2.0000,  2.15,   20.7},
{4.0000,  6.43,   39.8},
{10.000,  25.3,   89.8}, 
{15.000,  44.01,  1.2e2},  
};

class GSS_Elec_Optical : public PMIC_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0.0, PetscScalar Tl=1.0) const
  {
    complex<PetscScalar> n(1.0,0.0);
    if(lamda<Al[0].wavelength) 
      return complex<PetscScalar> (Al[0].RefractionIndexRe,Al[0].RefractionIndexIm);
    if(lamda>Al[table_size-1].wavelength) 
      return complex<PetscScalar> (Al[table_size-1].RefractionIndexRe,Al[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=Al[i].wavelength && lamda<=Al[i+1].wavelength)
       {
          complex<PetscScalar> n1(Al[i].RefractionIndexRe,Al[i].RefractionIndexIm);
          complex<PetscScalar> n2(Al[i+1].RefractionIndexRe,Al[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-Al[i].wavelength;
          PetscScalar d2 = Al[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2);
	  break;  
       }
    }
    return n;   
  }                                            
  
  // constructions
public:
  GSS_Elec_Optical(const PMIC_Environment &env):PMIC_Optical(env)
  {
     table_size = sizeof(Al)/sizeof(table);
     //the wave length should be scaled
     for(int i=0;i<table_size;i++)
     {
        Al[i].wavelength*=um;    
     }
  }

  ~GSS_Elec_Optical(){}
}
;

extern "C"
{
  PMIC_Optical*  PMIC_Elec_Optical_Default (const PMIC_Environment& env)
  {
    return new GSS_Elec_Optical(env);
  }
}
