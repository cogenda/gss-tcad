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
/*       A Two-Dimensional General Purpose Semiconductor HgCdTemulator.        */
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

typedef struct
{
   PetscScalar  wavelength;
   PetscScalar  RefractionIndexRe;
   PetscScalar  RefractionIndexIm;
} table;

//Source: FreeSnell hgcdte2
static table  HgCdTe[] =  {
{0.2066,  1.446,  1.70  }, //6eV
{0.2254,  1.446,  1.92  },
{0.2480,  1.868,  2.16  },
{0.2755,  2.452,  1.98  },
{0.3100,  2.505,  1.64  },
{0.3542,  2.53,   1.64  },
{0.4133,  2.834,  1.83  },
{0.4275,  3.007,  1.82  },
{0.4428,  3.135,  1.71  },
{0.4592,  3.192,  1.62  },
{0.4769,  3.23,   1.56  },
{0.4959,  3.275,  1.54  },
{0.5166,  3.361,  1.56  },
{0.5391,  3.577,  1.56  },
{0.5636,  3.798,  1.31  },
{0.5904,  3.81 ,  1.06  },
{0.6199,  3.763,  0.88  },
{0.6526,  3.704,  0.756 },
{0.6888,  3.648,  0.667 },
{0.7293,  3.601,  0.592 },
{0.7749,  3.558,  0.5391},
{0.8266,  3.518,  0.499 },
{0.8856,  3.504,  0.42  },
{0.9534,  3.5,    0.4   },
};

class GSS_HgCdTe_Optical : public PMIS_Optical
{
private:
  int   table_size;
public:
  complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=1.12, PetscScalar Tl=1.0) const
  {
    complex<PetscScalar> n(1.0,0.0);
    if(lamda<HgCdTe[0].wavelength) 
      return complex<PetscScalar> (HgCdTe[0].RefractionIndexRe,HgCdTe[0].RefractionIndexIm);
    if(lamda>HgCdTe[table_size-1].wavelength) 
      return complex<PetscScalar> (HgCdTe[table_size-1].RefractionIndexRe,HgCdTe[table_size-1].RefractionIndexIm);
    for(int i=0;i<table_size-1;i++)
    { 
       if(lamda>=HgCdTe[i].wavelength && lamda<=HgCdTe[i+1].wavelength)
       {
          complex<PetscScalar> n1(HgCdTe[i].RefractionIndexRe,HgCdTe[i].RefractionIndexIm);
          complex<PetscScalar> n2(HgCdTe[i+1].RefractionIndexRe,HgCdTe[i+1].RefractionIndexIm);
          PetscScalar d1 = lamda-HgCdTe[i].wavelength;
          PetscScalar d2 = HgCdTe[i+1].wavelength-lamda;
          n = (n1*d2 + n2*d1)/(d1+d2);
          break;  
       }
    }
    return n;   
  }           

// constructions
public:
  GSS_HgCdTe_Optical(const PMIS_Environment &env):PMIS_Optical(env)
  {
    table_size = sizeof(HgCdTe)/sizeof(table);
    //the wave length should be scaled
    for(int i=0;i<table_size;i++)
    {
        HgCdTe[i].wavelength*=um;    
    }
  }

  ~GSS_HgCdTe_Optical()
  {}
}
;

extern "C"
{
  PMIS_Optical*  PMIS_HgCdTe_Optical_Default (const PMIS_Environment& env)
  {
    return new GSS_HgCdTe_Optical(env);
  }
}

