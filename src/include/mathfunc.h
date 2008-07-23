#ifndef __mathfunc_h
#define __mathfunc_h
#include "petsc.h"
#include "brkpnts.h"

#ifdef HAVE_GSL 
#include <gsl/gsl_sf_fermi_dirac.h>
#endif

#ifdef LINUX 
#include <fenv.h>
#endif

#include "adolc.h" 
#include <math.h>
using namespace adtl;

/* define the constant */
const PetscScalar PI     = 3.14159265358979323846;
const PetscScalar SQRTPI = 1.772453850905516;
const PetscScalar LOG10e = 0.434294481903251827651;
const PetscScalar LOGe10 = 2.30258509299404568402;
const PetscScalar VerySmallNumericValue = 1.0e-30;
const PetscScalar VeryLargeNumericValue = 1.0e30;
const PetscScalar MaximumExponent = 76.0;
const PetscScalar MinimumLogarithmArgument = 1.0e-40;
const PetscScalar MinimumExponent = 1.0e-5;

inline PetscScalar dmax(PetscScalar a,PetscScalar b)
{ return a>b? a:b;}

inline PetscScalar dmin(PetscScalar a,PetscScalar b)
{ return a<b? a:b;}

inline PetscScalar dsign(PetscScalar a)
{return a>=0 ? 1.0 : -1.0;}

inline int isign(int a)
{return a>0? 1:-1;}

/* ----------------------------------------------------------------------------
 * NOTE:  the range limits for the Bernoulli and Aux functions were generated
 *        by the "brkpnts" program and stored in brkpnts.h.
 * --------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------
 * bern:  This function returns the Bernoulli function of the argument.  To
 * avoid under and overflows this function is defined by equivalent or approx-
 * imate functions depending upon the value of the argument.
 *
 *                                        x
 *                              B(x) = -------
 *                                     e^x - 1
 *
 */
inline PetscScalar bern ( PetscScalar x )
{
  PetscScalar y;

  if (x <= BP0_BERN)
  { return(-x); }
  else if (x <  BP1_BERN)
  { return(x / (exp(x) - 1.0)); }
  else if (x <= BP2_BERN)
  { return(1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0))); }
  else if (x <  BP3_BERN)
  { y = exp(-x);   return((x * y) / (1.0 - y)); }
  else if (x <  BP4_BERN)
  { return(x * exp(-x)); }
  else { return 0; }

} /* bern */

inline AutoDScalar bern ( const AutoDScalar &x )
{
  AutoDScalar y;

  if (x <= BP0_BERN)
  { return(-x); }
  else if (x <  BP1_BERN)
  { return(x / (exp(x) - 1.0)); }
  else if (x <= BP2_BERN)
  { return(1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0))); }
  else if (x <  BP3_BERN)
  { y = exp(-x);   return((x * y) / (1.0 - y)); }
  else if (x <  BP4_BERN)
  { return(x * exp(-x)); }
  else { return 0; }

} /* bern */


/* ----------------------------------------------------------------------------
 * pd1bern:  This function returns the total derivative of the Bernoulli
 * function with respect to the argument.  To avoid under and overflows this
 * function is defined by equivalent or approximate functions depending upon
 * the value of the argument.
 *
 *                         d        (1-x)*e^x - 1
 *                         --B(x) = -------------
 *                         dx        (e^x - 1)^2
 */
inline PetscScalar pd1bern ( PetscScalar x )
{
  PetscScalar y, z;

  if (x <= BP0_DBERN)
  { return(-1.0); }
  else if (x <= BP1_DBERN)
  { return((1.0 - x)*exp(x) - 1.0); }
  else if (x <  BP2_DBERN)
  { y = exp(x); z = y - 1.0; return(((1.0 - x)*y - 1.0)/(z*z)); }
  else if (x <= BP3_DBERN)
  { return(-0.5 + x/6.0 * (1.0 - x*x/30.0)); }
  else if (x <  BP4_DBERN)
  { y = exp(-x); z = 1 - y;  return(((1.0 - x)*y - y*y )/(z*z)); }
  else if (x <  BP5_DBERN)
  { y = exp(-x); return((1.0 - x)*y - y*y); }
  else { return(0.0); }

} /* pd1bern */


/* ----------------------------------------------------------------------------
 * aux1:  This function returns the aux1 function.  To avoid under and over-
 * flows this function is defined by equivalent or approximate functions
 * depending upon the value of the argument.
 *
 *                                       x
 *                         Aux1(x) =  -------
 *                                    sinh(x)
 */
inline PetscScalar aux1 ( PetscScalar x )
{
  PetscScalar y=x;
  if      (x < -BP0_MISC) y = -BP0_MISC;
  else if (x >  BP0_MISC) y =  BP0_MISC;

  if (y <= BP0_AUX1)
  { return(y / sinh(y)); }
  else if (y <= BP1_AUX1)
  { return(1 - y*y/6.0*(1.0 - 7.0*y*y/60.0)); }
  else { return(y / sinh(y)); }
} /* aux1 */

inline AutoDScalar aux1 ( const AutoDScalar &x )
{
  AutoDScalar y=x;
  if      (x < -BP0_MISC) y = -BP0_MISC;
  else if (x >  BP0_MISC) y =  BP0_MISC;

  if (y <= BP0_AUX1)
  { return(y / sinh(y)); }
  else if (y <= BP1_AUX1)
  { return(1 - y*y/6.0*(1.0 - 7.0*y*y/60.0)); }
  else { return(y / sinh(y)); }
} /* aux1 */


/* ----------------------------------------------------------------------------
 * pd1aux1:  This function returns the total derivative of the aux1 function.
 * To avoid under and overflows this function is defined by equivalent or
 * approximate functions depending upon the value of the argument.
 *
 *                    d           sinh(x) - x*cosh(x)
 *                    --Aux1(x) = -------------------
 *                    dx              (sinh(x))^2
 */
inline PetscScalar pd1aux1 ( PetscScalar x )
{
  PetscScalar y=x;
  PetscScalar z;
  if      (x < -BP0_MISC) y = -BP0_MISC;
  else if (x >  BP0_MISC) y =  BP0_MISC;

  if (y <= BP0_DAUX1)
  { z = sinh(y); return((z - y*cosh(y))/(z*z)); }
  else if (y <= BP1_DAUX1)
  { return(-y/3.0*(1.0 - 7.0*y*y/30.0)); }
  else { z = sinh(y); return((z - y*cosh(y))/(z*z)); }
} /* pd1aux1 */


/* ----------------------------------------------------------------------------
 * aux2:  This function returns the aux2 function.  To avoid under and over-
 * flows this function is defined by equivalent or approximate functions
 * depending upon the value of the argument.
 *
 *                                           1
 *                              Aux2(x) = -------
 *                                        1 + e^x
 */
inline PetscScalar aux2 ( PetscScalar x )
{
  if (x <= BP0_AUX2)
  { return(1.0); }
  else if (x <= BP1_AUX2)
  { return(1.0 / (1.0 + exp(x))); }
  else if (x <= BP2_AUX2)
  { return(exp(-x)); }
  else { return(0.0); }

} /* aux2 */

inline AutoDScalar aux2 ( const AutoDScalar &x )
{
  if (x <= BP0_AUX2)
  { return(1.0); }
  else if (x <= BP1_AUX2)
  { return(1.0 / (1.0 + exp(x))); }
  else if (x <= BP2_AUX2)
  { return(exp(-x)); }
  else { return(0.0); }

} /* aux2 */


/* ----------------------------------------------------------------------------
 * pd1aux2:  This function returns the total derivative of the aux2 function.
 * To avoid under and overflows this function is defined by equivalent or
 * approximate functions depending upon the value of the argument.
 *
 *                         d             - e^x
 *                         --Aux2(x) = -----------
 *                         dx          (1 + e^x)^2
 */
inline PetscScalar pd1aux2 ( PetscScalar x )
{
  PetscScalar y,z;

  if (x <= BP0_DAUX2)
  { return(0.0); }
  else if (x <= BP1_DAUX2)
  { return(-exp(x)); }
  else if (x <= BP2_DAUX2)
  { y = exp(x); z = y + 1.0; return(-y/(z*z)); }
  else if (x <= BP3_DAUX2)
  { return(-exp(-x)); }
  else { return(0.0); }

} /* pd1aux2 */


/* ----------------------------------------------------------------------------
 * pd1erf:  This function returns the derivative of the error function with
 * respect to the first variable.
 */
inline PetscScalar pd1erf ( PetscScalar x )
{
  return 2.0 / sqrt(PI) * exp(-x*x);
}


/* ----------------------------------------------------------------------------
 * fermi_half:  This function returns value of 1/2 order Fermi-Dirac Integral
 */
inline PetscScalar fermi_half(PetscScalar x)
{
#ifdef HAVE_GSL 
  return PetscScalar(gsl_sf_fermi_dirac_half(double(x)));
#else
  /* use an analytic expression. The result achieves within 0.4% error in all ranges.*/
  if(x<-4.5)
  {
    //for small arguments, fhfp and exp are almost identical.
    return 1.0/exp(-x);
  }
  else if(x<0.0)
  {
    PetscScalar v = std::pow(x,4) + 50 + 33.6*x*(1-0.68*exp(-0.17*(x+1)*(x+1)));
    PetscScalar p = 1.329340388179*std::pow(v,PetscScalar(-0.375));
    return 1.0/(exp(-x) + p);
  }
  else
  {
    PetscScalar v = std::pow(x,4) + 50 + 33.6*x*(1-0.68*exp(-0.17*(x+1)*(x+1)));
    PetscScalar p = 1.329340388179*std::pow(v,PetscScalar(-0.375));
    return 1.0/(1.0/exp(x) + p);
  }
#endif  
}




/*-----------------------------------------------------------------------
 *
 *     fhfm evaluates the fermi-dirac integral of minus one-half order 
 *     f-1/2(x) from x
 */
inline PetscScalar fermi_mhalf(PetscScalar x)
{
#ifdef HAVE_GSL 
  return PetscScalar(gsl_sf_fermi_dirac_mhalf(double(x)));
#else 
/*     the formulae used are after j.s. blakemore, ''semiconductor
 *     statistics,'' new york: pergamon press, appendix c, 1962, for
 *     x greater than or equal to 5.5; when x is less than 5.5, 
 *     f-1/2=1/(1/f+a+2*b*f+3*c*f*f+4*d*f*f*f) where f=f1/2 and it
 *     is the derivative of approximate form of f1/2 (see function
 *     fermi_half)
 *     the maximum relative error is about 1.2%
 */  
  PetscScalar a=3.53553e-1,b=4.95009e-3,c=1.48386e-4,d=4.42563e-6;
  if(x<-4.5)
  {
    //for small argument, fhfm1 and exp are almost identical.
    return exp(x);
  }
  else
  {
    if(x<5.5)
    {
      PetscScalar f=fermi_half(x);
      return f/(1.0+f*(a+f*(-2.0*b+f*(3*c-4*d*f))));
    }
    else
      return 2.0*x/(SQRTPI*std::pow((x*x+0.6),PetscScalar(0.25)));
  }
#endif  
}


/*-----------------------------------------------------------------------
 *
 *     fhfm2 evaluates the fermi-dirac integral of minus one-half order 
 *     f-1/2(eta) from x=f1/2
 *     when x<10, the formula 
 *       f-1/2=1/(1/x+a+2*b*x+3*c*x*x+4*d*x*x*x), where x=f1/2, 
 *     is used.
 *     when x>=10, first from x=f1/2(eta) to evaluate eta by
 *       eta=sqrt((3/4/sqrt(pi)*x)**(4/3)-pi*pi/6),
 *     then to evaluate f-1/2 using
 *       f-1/2 = 2*eta/sqrt(pi)/(eta*eta+0.6)**(1/4)
 *     the maximum relative error is about 1.2%
 */
inline PetscScalar fermi_mhalf_f(PetscScalar x)
{
  PetscScalar a=3.53553e-1,b=4.95009e-3,c=1.48386e-4,d=4.42563e-6;
  if(x<1.0e1)
    return x/(1.e0+x*(a+x*(-2.e0*b+x*(3.e0*c-4.e0*d*x))));
  else
  {
    PetscScalar eta=sqrt(std::pow(0.75e0*SQRTPI*x,PetscScalar(4.e0/3.e0))-PI*PI/6.e0);
    return 2.e0*eta/(SQRTPI*std::pow(eta*eta+0.6e0,PetscScalar(0.25e0)));
  }
}



/* ----------------------------------------------------------------------------
 * inv_fermi_half:  This function returns inverse value of Fermi-Dirac Integral
 */
inline PetscScalar inv_fermi_half(PetscScalar x)
{
  if(x<8.463)
    return log(x) + x*(3.5355339059327379e-001 - x*(4.9500897298752622e-003
                       - x*(1.4838577128872821e-004 - x*4.4256301190009895e-006)));
  else
    return sqrt(std::pow(0.75*SQRTPI*x,PetscScalar(4.0/3.0)) - PI*PI/6.0);
}


/*-----------------------------------------------------------------------
 *   GAMMA calculates f1/2(eta)/exp(eta) according to the approximate 
 *   formula in casey's book,dummy arguement x=f1/2(eta).
 */
inline  PetscScalar gamma_f(PetscScalar x)
{
  const PetscScalar a=3.53553e-1,b=4.95009e-3,c=1.48386e-4;
  const PetscScalar d=4.42563e-6,pi1=1.772453851e0,pi2=9.869604401e0;
  PetscScalar temx;
  if(x>1.0e1)
  {
    temx=sqrt(std::pow(7.5e-1*pi1*x,PetscScalar(4.e0/3.e0))-pi2/6.e0);
    if(x > MaximumExponent) 
      return VerySmallNumericValue;
    else
      return x/exp(temx);
  }
  else if(x>0.0)
  {
    temx=x*(a+x*(-b+x*(c-x*d)));
    return 1.0/exp(temx);
  }
  else
    return 1.0;
}


inline  AutoDScalar gamma_f(const AutoDScalar &x)
{
  const PetscScalar a=3.53553e-1,b=4.95009e-3,c=1.48386e-4;
  const PetscScalar d=4.42563e-6,pi1=1.772453851e0,pi2=9.869604401e0;
  AutoDScalar temx;
  if(x>1.0e1)
  {
    temx=sqrt(pow(7.5e-1*pi1*x,PetscScalar(4.e0/3.e0))-pi2/6.e0);
    if(x > MaximumExponent) 
      return VerySmallNumericValue;
    else
      return x/exp(temx);
  }
  else if(x>0.0)
  {
    temx=x*(a+x*(-b+x*(c-x*d)));
    return 1.0/exp(temx);
  }
  else
    return 1.0;
}

#endif

