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
// Material Type: Amorphous Silicon


#include "PMI.h"

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;


static table  PolySi[] = {
{0.1240,  0.466,  1.090},
{0.1378,  0.499,  1.330},
{0.1550,  0.554,  1.660},
{0.1771,  0.670,  2.080},
{0.2066,  0.961,  2.650},
{0.2480,  1.660,  3.380},
{0.2755,  2.310,  3.710},
{0.3100,  3.360,  3.920},
{0.3543,  4.590,  3.380},
{0.4133,  5.430,  2.190},
{0.4960,  5.250,  0.992},
{0.6199,  4.710,  0.217},
{0.8190,  3.900,  2.23e-2},
{0.9613,  3.660,  1.1e-2},
};

class GSS_PolySi_Optical : public PMIC_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0.0, PetscScalar Tl=1.0) const 
  {
    complex<PetscScalar> n(1.0,0.0);
    if(lamda<PolySi[0].wavelength) 
      return complex<PetscScalar> (PolySi[0].RefractionIndexRe,PolySi[0].RefractionIndexIm);
    if(lamda>PolySi[table_size-1].wavelength) 
      return complex<PetscScalar> (PolySi[table_size-1].RefractionIndexRe,PolySi[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=PolySi[i].wavelength && lamda<=PolySi[i+1].wavelength)
       {
          complex<PetscScalar> n1(PolySi[i].RefractionIndexRe,PolySi[i].RefractionIndexIm);
          complex<PetscScalar> n2(PolySi[i+1].RefractionIndexRe,PolySi[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-PolySi[i].wavelength;
          PetscScalar d2 = PolySi[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2);
	  break;  
       }
    }
    return n;   
  }                                            
  
  // constructions
public:
  GSS_PolySi_Optical(const PMIC_Environment &env):PMIC_Optical(env)
  {
     table_size = sizeof(PolySi)/sizeof(table);
     //the wave length should be scaled
     for(int i=0;i<table_size;i++)
     {
        PolySi[i].wavelength*=um;    
     }
  }

  ~GSS_PolySi_Optical(){}
}
;

extern "C"
{
  PMIC_Optical*  PMIC_PolySi_Optical_Default (const PMIC_Environment& env)
  {
    return new GSS_PolySi_Optical(env);
  }
}