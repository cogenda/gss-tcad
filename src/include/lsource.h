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
/*  Last update: Nov 29, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//sevel types of light source are defined here
//reference: spice vsource model

#ifndef _lsource_h_
#define _lsource_h_
#include "mathfunc.h"
#include <math.h>
#include <dlfcn.h>
#include "typedef.h"

enum LSOURCE_CARD_ERROR
{
  LSOURCE_NO_ERROR,
  LSOURCE_UNKNOW_PARAMETER   
};

//virtual base class
class LSource
{
public:  
  virtual double currentft(double t)=0;
  virtual ~LSource(){}
};


class LSource_UNIFORM: public LSource
{
private:
  double td;  
  double power;
public:
  LSource_UNIFORM(double t1,double p1): td(t1), power(p1){};
  double currentft(double t)
  { return t>=td? power:0;};
};


class LSource_PULSE: public LSource
{
private:
  double td;
  double tr;
  double tf;
  double pw;
  double pr;  
  double powerlo,powerhi;
public:
  LSource_PULSE(double t1,double t2,double t3,double t4, double t5,double plo,double phi):
  td(t1),tr(t2),tf(t3),pw(t4),pr(t5),powerlo(plo),powerhi(phi){};
  double currentft(double t)
  {
    if(t<td)
      return powerlo;
    else
    {
      t-=td;
      while(t>pr) t-=pr;
      if(t<tr)
        return powerlo+t*(powerhi-powerlo)/tr;
      else if(t<tr+pw)
        return powerhi;
      else if(t<tr+pw+tf)
        return powerhi-(t-tr-pw)*(powerhi-powerlo)/tf;
      else    return powerlo;
    }

  }
};

//just for single pulse
class LSource_GAUSSIAN: public LSource
{
private:
  double tp;
  double tb; 
  double power;  
public:
  LSource_GAUSSIAN(double t1,double t2,double p1):
  tp(t1),tb(t2),power(p1){};
  double currentft(double t)
  {
    return power*2*exp(-(t-tp)*(t-tp)/tb)/(tb*sqrt(PI)*erfc(tp/tb));
  }
};



class LSource_SHELL: public LSource
{
private:
  void   *dll;
  double (*Lsource_Shell)(double);
  double scale_t;
public:
  LSource_SHELL(void * dp, void * fp, double s_t)
  {
     dll = dp;
     Lsource_Shell = (double (*)(double)) fp;
     scale_t = s_t;
  }
  double currentft(double t)
  {
     return Lsource_Shell(t/scale_t);
  }
  ~LSource_SHELL()  {dlclose(dll);}
};

#endif
