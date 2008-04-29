/*****************************************************************************/
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
/*  GSS 0.4x                                                                 */
/*  Last update: May 07, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _flux3e_h_
#define _flux3e_h_
#include "petsc.h"
#include "mathfunc.h"


/* ----------------------------------------------------------------------------
 * Theta function
 *
 *                                           T2-T1
 *                         Theta(T1,T2) = -------------
 *                                         log(T2/T1)
 */
inline PetscScalar Theta(PetscScalar T1, PetscScalar T2)
{
        PetscScalar x = T2/T1-1;
        if(fabs(x)>1e-6)
                return (T2-T1)/log(fabs(T2/T1)); 
        else
                return T1/(1-0.5*x);    
}       
inline AutoDScalar Theta(const AutoDScalar &T1, const AutoDScalar &T2)
{
        AutoDScalar x = T2/T1-1;
        if(fabs(x)>1e-6)
                return (T2-T1)/log(fabs(T2/T1)); 
        else
                return T1/(1-0.5*x);    
}




//-----------------------------------------------------------------------------

inline PetscScalar In(PetscScalar kb, PetscScalar e, PetscScalar V1, PetscScalar V2,  
                      PetscScalar n1,PetscScalar n2, PetscScalar Tn1,PetscScalar Tn2, PetscScalar h)
{
  PetscScalar theta = Theta(Tn1,Tn2);
  PetscScalar alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  return kb*0.5*(Tn1+Tn2)*theta*(bern(alpha)*n2/Tn2 - bern(-alpha)*n1/Tn1)/h;
}

inline AutoDScalar In(PetscScalar kb, PetscScalar e, const AutoDScalar &V1, const AutoDScalar &V2,
                      const AutoDScalar &n1,const AutoDScalar &n2, const AutoDScalar &Tn1,const AutoDScalar &Tn2, PetscScalar h)
{
  AutoDScalar theta = Theta(Tn1,Tn2);
  AutoDScalar alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  return kb*0.5*(Tn1+Tn2)*theta*(bern(alpha)*n2/Tn2 - bern(-alpha)*n1/Tn1)/h;
}



//-----------------------------------------------------------------------------

inline PetscScalar Ip(PetscScalar kb, PetscScalar e, PetscScalar V1, PetscScalar V2,  
                      PetscScalar p1,PetscScalar p2, PetscScalar Tp1,PetscScalar Tp2, PetscScalar h)
{
  PetscScalar theta = Theta(Tp1,Tp2);
  PetscScalar alpha = (e/kb*(V2-V1)+2*(Tp2-Tp1))/theta;
  return kb*0.5*(Tp1+Tp2)*theta*(bern(alpha)*p1/Tp1 - bern(-alpha)*p2/Tp2)/h;  
}

inline AutoDScalar Ip(PetscScalar kb, PetscScalar e, const AutoDScalar &V1, const AutoDScalar &V2,  
                      const AutoDScalar &p1,const AutoDScalar &p2, const AutoDScalar &Tp1,const AutoDScalar &Tp2, PetscScalar h)
{
  AutoDScalar theta = Theta(Tp1,Tp2);
  AutoDScalar alpha = (e/kb*(V2-V1)+2*(Tp2-Tp1))/theta;
  return kb*0.5*(Tp1+Tp2)*theta*(bern(alpha)*p1/Tp1 - bern(-alpha)*p2/Tp2)/h;  
}




//-----------------------------------------------------------------------------

inline PetscScalar Sn(PetscScalar kb, PetscScalar e, PetscScalar V1, PetscScalar V2,
                      PetscScalar n1,PetscScalar n2, PetscScalar Tn1,PetscScalar Tn2, PetscScalar h)
{
  PetscScalar theta = Theta(Tn1,Tn2);
  PetscScalar alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  PetscScalar phi   = (e/kb*(V2-V1)-(Tn2-Tn1))/theta-log(fabs(n2/n1));
  PetscScalar Dn    = kb*0.5*(Tn1+Tn2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dn/h*theta*( - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
  return -2.0*kb*Dn/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*n2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
}

inline AutoDScalar Sn(PetscScalar kb, PetscScalar e, const  AutoDScalar &V1,const  AutoDScalar &V2,
                      const AutoDScalar &n1, const AutoDScalar &n2, const AutoDScalar &Tn1,const AutoDScalar &Tn2, PetscScalar h)
{
  AutoDScalar theta = Theta(Tn1,Tn2);
  AutoDScalar alpha = (e/kb*(V2-V1)-2*(Tn2-Tn1))/theta;
  AutoDScalar phi   = (e/kb*(V2-V1)-(Tn2-Tn1))/theta-log(fabs(n2/n1));
  AutoDScalar Dn    = kb*0.5*(Tn1+Tn2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dn/h*theta*( - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
  return -2.0*kb*Dn/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*n2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*n1);
}



//-----------------------------------------------------------------------------

inline PetscScalar Sp(PetscScalar kb, PetscScalar e, PetscScalar V1, PetscScalar V2,  
                      PetscScalar p1,PetscScalar p2, PetscScalar Tp1,PetscScalar Tp2, PetscScalar h)
{
  PetscScalar theta = Theta(Tp1,Tp2);
  PetscScalar alpha = (-e/kb*(V2-V1)-2*(Tp2-Tp1))/theta;
  PetscScalar phi   = (-e/kb*(V2-V1)-(Tp2-Tp1))/theta-log(fabs(p2/p1));
  PetscScalar Dp    = kb*0.5*(Tp1+Tp2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dp/h*theta*(- bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
  return  -2.0*kb*Dp/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*p2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
}

inline AutoDScalar Sp(PetscScalar kb, PetscScalar e, const  AutoDScalar &V1,const  AutoDScalar &V2,  
                      const AutoDScalar &p1, const AutoDScalar &p2, const AutoDScalar &Tp1,const AutoDScalar &Tp2, PetscScalar h)
{
  AutoDScalar theta = Theta(Tp1,Tp2);
  AutoDScalar alpha = (-e/kb*(V2-V1)-2*(Tp2-Tp1))/theta;
  AutoDScalar phi   = (-e/kb*(V2-V1)-(Tp2-Tp1))/theta-log(fabs(p2/p1));
  AutoDScalar Dp    = kb*0.5*(Tp1+Tp2)/e;
  if(alpha > BP4_BERN || 1.25*phi > BP4_BERN)
    return -2.0*kb*Dp/h*theta*(- bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
  return   -2.0*kb*Dp/h*theta*(bern(alpha)*bern(1.25*phi)/bern(phi)*p2 - bern(-alpha)*bern(-1.25*phi)/bern(-phi)*p1);
}



#endif
