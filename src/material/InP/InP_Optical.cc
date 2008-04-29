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
/*       A Two-Dimensional General Purpose Semiconductor InPmulator.         */
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
// Material Type: InP


#include "PMI.h"

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;


static table  InP[] =  {
{0.0620,  0.793,  0.494 },//20eV
{0.0827,  0.695,  0.574 },
{0.1240,  0.797,  1.154 },
{0.1378,  0.872,  1.304 },
{0.1550,  0.960,  1.566 },
{0.1771,  1.215,  1.915 },
{0.2066,  1.500,  2.063 },
{0.2254,  1.426,  2.562 },
{0.2480,  2.131,  3.495 },
{0.2755,  3.697,  2.186 },
{0.3100,  3.141,  1.730 },
{0.3542,  3.193,  1.948 },
{0.4133,  4.395,  1.247 },
{0.4275,  4.256,  0.964 },
{0.4428,  4.121,  0.786 },
{0.4592,  4.004,  0.667 },
{0.4769,  3.903,  0.579 },
{0.4959,  3.818,  0.511 },
{0.5166,  3.745,  0.457 },
{0.5391,  3.682,  0.416 },
{0.5636,  3.629,  0.380 },
{0.5904,  3.585,  0.347 },
{0.6199,  3.549,  0.317 },
{0.6526,  3.517,  0.293 },
{0.6888,  3.492,  0.270 },
{0.7293,  3.476,  0.242 },
{0.7749,  3.467,  0.218 },
{0.8266,  3.456,  0.203 },
{0.8856,  3.440,  0.130 },
{0.9534,  3.400,  7.90e-2},
{0.9759,  3.362,  2.81e-4},//1.27eV
};

class GSS_InP_Optical : public PMIS_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=1.12, PetscScalar Tl=1.0) const
  {
    complex<PetscScalar> n(1.0,0.0);
    if(lamda<InP[0].wavelength) 
      return complex<PetscScalar> (InP[0].RefractionIndexRe,InP[0].RefractionIndexIm);
    if(lamda>InP[table_size-1].wavelength) 
      return complex<PetscScalar> (InP[table_size-1].RefractionIndexRe,InP[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=InP[i].wavelength && lamda<=InP[i+1].wavelength)
       {
          complex<PetscScalar> n1(InP[i].RefractionIndexRe,InP[i].RefractionIndexIm);
          complex<PetscScalar> n2(InP[i+1].RefractionIndexRe,InP[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-InP[i].wavelength;
          PetscScalar d2 = InP[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2);
          break;  
       }
    }
    return n;   
  }           

// constructions
public:
  GSS_InP_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    table_size = sizeof(InP)/sizeof(table);
    //the wave length should be scaled
    for(int i=0;i<table_size;i++)
    {
        InP[i].wavelength*=um;    
    }
  }

  ~GSS_InP_Optical()
  {}
}
;

extern "C"
{
  PMIS_Optical*  PMIS_InP_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_InP_Optical(env);
  }
}

