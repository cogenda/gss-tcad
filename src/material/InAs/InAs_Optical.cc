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
// Material Type: InAs


#include "PMI.h"

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;

//Source: FreeSnell
static table  InAs[] =  {
{1.00, 3.548,  0.278},       // 1.24eV
{1.05, 3.558,  0.294},       // 1.18eV
{1.10, 3.572,  0.309},       // 1.13eV
{1.15, 3.591,  0.324},       // 1.08eV
{1.20, 3.613,  0.340},       // 1.03eV
{1.30, 3.658,  0.371},       // 0.954eV
{1.4,  3.696,  0.410},       // 0.8856eV
{1.5,  3.714,  0.432},       // 0.8266eV
{1.6,  3.755,  0.463},       // 0.7749eV
{1.7,  3.798,  0.493},       // 0.7293eV   
{1.8,  3.851,  0.530},       // 0.6888eV
{1.9,  3.917,  0.572},       // 0.6526eV
{2.0,  3.995,  0.634},       // 0.6199eV
{2.1,  4.088,  0.712},       // 0.5904eV
{2.2,  4.199,  0.822},       // 0.5635eV
{2.3,  4.331,  0.991},       // 0.5391eV
{2.4,  4.466,  1.283},       // 0.5166eV
{2.5,  4.364,  1.786},       // 0.4959eV
{2.6,  4.021,  1.885},       // 0.4769eV
{2.7,  3.911,  2.016},       // 0.4592eV
{2.8,  3.626,  2.208},       // 0.4428eV
{2.9,  3.337,  2.129},       // 0.4275eV
{3.0,  3.197,  2.034},       // 0.4133eV
{3.5,  3.008,  1.754},       // 0.3542eV
{4.0,  3.313,  1.799},       // 0.3100eV
{4.5,  3.194,  3.445},       // 0.2755eV
{5.0,  1.524,  2.871},       // 0.2480eV
{5.5,  1.282,  2.344},       // 0.2254eV
{6.0,  1.434,  2.112},       // 0.2066eV
{7.0,  1.282,  1.754},       // 0.1771eV
{8.0,  0.984,  1.515},       // 0.1550eV
{9.0,  0.897,  1.231},       // 0.1378eV
{10.0,  0.835,  1.071},      // 0.1240eV
{15.0,  0.894,  0.336},      // 0.0827eV
{20.0,  1.125,  0.225},      // 0.0620eV
};

class GSS_InAs_Optical : public PMIS_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=1.12, PetscScalar Tl=1.0) const
  {
    complex<PetscScalar> n(1.0,0.0);
    if(lamda<InAs[0].wavelength) 
      return complex<PetscScalar> (InAs[0].RefractionIndexRe,InAs[0].RefractionIndexIm);
    if(lamda>InAs[table_size-1].wavelength) 
      return complex<PetscScalar> (InAs[table_size-1].RefractionIndexRe,InAs[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=InAs[i].wavelength && lamda<=InAs[i+1].wavelength)
       {
          complex<PetscScalar> n1(InAs[i].RefractionIndexRe,InAs[i].RefractionIndexIm);
          complex<PetscScalar> n2(InAs[i+1].RefractionIndexRe,InAs[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-InAs[i].wavelength;
          PetscScalar d2 = InAs[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2);
          break;  
       }
    }
    return n;   
  }           

// constructions
public:
  GSS_InAs_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    table_size = sizeof(InAs)/sizeof(table);
    //the wave length should be scaled
    for(int i=0;i<table_size;i++)
    {
        InAs[i].wavelength*=um;    
    }
  }

  ~GSS_InAs_Optical()
  {}
}
;

extern "C"
{
  PMIS_Optical*  PMIS_InAs_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_InAs_Optical(env);
  }
}

