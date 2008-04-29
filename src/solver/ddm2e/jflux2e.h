/*****************************************************************************/
/*   	        8888888         88888888         88888888                    */
/*  	      8                8                8                            */
/* 	     8                 8                8                            */
/*  	     8                  88888888         88888888                    */
/* 	     8      8888                8                8                   */
/* 	      8       8                 8                8                   */
/* 	        888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS 0.4x                                                                 */
/*  Last update: March 22, 2006                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _jflux2_h_
#define _jflux2_h_
#include "petsc.h"


inline PetscScalar nmid(PetscScalar kb,PetscScalar e,PetscScalar dV,PetscScalar n1,PetscScalar n2,PetscScalar T,PetscScalar dT)
{
  PetscScalar Vt = kb*T/e;
  PetscScalar alpha = -dV/(2*Vt)+ dT/(2*T);
  return n1*aux2(alpha) + n2*aux2(-alpha);;
}

inline AutoDScalar nmid(PetscScalar kb,PetscScalar e, const AutoDScalar &dV, const AutoDScalar &n1, const AutoDScalar &n2, 
                        const AutoDScalar &T, const AutoDScalar &dT)
{
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)+ dT/(2*T);
  return n1*aux2(alpha) + n2*aux2(-alpha);;
}



//-----------------------------------------------------------------------------
inline PetscScalar pmid(PetscScalar kb,PetscScalar e,PetscScalar dV,PetscScalar p1,PetscScalar p2,PetscScalar T,PetscScalar dT)
{
  PetscScalar Vt = kb*T/e;
  PetscScalar alpha = -dV/(2*Vt)- dT/(2*T);
  return p1*aux2(-alpha) + p2*aux2(alpha);
}

inline AutoDScalar pmid(PetscScalar kb,PetscScalar e, const AutoDScalar &dV, const AutoDScalar &p1, const AutoDScalar &p2,
                        const AutoDScalar &T, const AutoDScalar &dT)
{
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)- dT/(2*T);
  return p1*aux2(-alpha) + p2*aux2(alpha);
}



//-----------------------------------------------------------------------------

inline PetscScalar In(PetscScalar kb,PetscScalar e,PetscScalar dV,PetscScalar n1,PetscScalar n2,PetscScalar T, PetscScalar dT,PetscScalar h)
{
  PetscScalar E  = -dV/h;
  PetscScalar Vt = kb*T/e;
  PetscScalar alpha = -dV/(2*Vt)+ dT/(2*T);
  PetscScalar n  = n1*aux2(alpha) + n2*aux2(-alpha);
  PetscScalar dndx = aux1(alpha)*(n2-n1)/h;
  return (E*n + Vt*dndx + kb*n/e*dT/h);
}

inline AutoDScalar In(PetscScalar kb,PetscScalar e, const AutoDScalar &dV, const AutoDScalar &n1, const AutoDScalar &n2,
                      const AutoDScalar &T, const AutoDScalar &dT, PetscScalar h)
{
  AutoDScalar E  = -dV/h;
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)+ dT/(2*T);
  AutoDScalar n  = n1*aux2(alpha) + n2*aux2(-alpha);
  AutoDScalar dndx = aux1(alpha)*(n2-n1)/h;
  return (E*n + Vt*dndx + kb*n/e*dT/h);
}



//-----------------------------------------------------------------------------

inline PetscScalar Ip(PetscScalar kb,PetscScalar e,PetscScalar dV,PetscScalar p1,PetscScalar p2,PetscScalar T, PetscScalar dT,PetscScalar h)
{
  PetscScalar E  = -dV/h;
  PetscScalar Vt = kb*T/e;
  PetscScalar alpha = -dV/(2*Vt)- dT/(2*T);
  PetscScalar p  = p1*aux2(-alpha) + p2*aux2(alpha);
  PetscScalar dpdx = aux1(alpha)*(p2-p1)/h;
  return (E*p-Vt*dpdx - kb*p/e*dT/h);
}

inline AutoDScalar Ip(PetscScalar kb,PetscScalar e, const AutoDScalar &dV, const AutoDScalar &p1, const AutoDScalar &p2,
                      const AutoDScalar &T, const AutoDScalar &dT,PetscScalar h)
{
  AutoDScalar E  = -dV/h;
  AutoDScalar Vt = kb*T/e;
  AutoDScalar alpha = -dV/(2*Vt)- dT/(2*T);
  AutoDScalar p  = p1*aux2(-alpha) + p2*aux2(alpha);
  AutoDScalar dpdx = aux1(alpha)*(p2-p1)/h;
  return (E*p-Vt*dpdx - kb*p/e*dT/h);
}


#endif
