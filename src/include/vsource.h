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
//sevel types of vsource are defined here
//reference: spice vsource model

#ifndef _vsource_h_
#define _vsource_h_
#include <math.h>
#include <dlfcn.h>
#include "typedef.h"

enum VSOURCE_CARD_ERROR
{
  VSOURCE_NO_ERROR,
  VSOURCE_NO_ID,
  VSOURCE_UNKNOW_PARAMETER
};

//virtual base class
class VSource
{
public:
  char label[32];
  virtual double vapp(double t)=0;
  virtual ~VSource(){}
};


class VDC: public VSource
{
private:
  double td;
  double Vdc;
public:
  VDC(double t1,double v1): td(t1),Vdc(v1){};
  double vapp(double t)
{ return t>=td? Vdc:0;};
};



class VSIN: public VSource
{
private:
  double td;
  double V0;
  double Vamp;
  double fre;
  double alpha;
public:
  VSIN(double t1,double v0,double v1,double f1,double a1): td(t1),V0(v0),Vamp(v1),fre(f1),alpha(a1){};
  double vapp(double t)
  { return t>=td? V0+Vamp*exp(-alpha*(t-td))*sin(2*3.14159265359*fre*(t-td)):V0;};
};


class VPULSE: public VSource
{
private:
  double td;
  double tr;
  double tf;
  double pw;
  double pr;
  double Vlo,Vhi;
public:
  VPULSE(double t1,double v1,double v2,double t2,double t3,double t4, double t5):
  td(t1),tr(t2),tf(t3),pw(t4),pr(t5),Vlo(v1),Vhi(v2){};
  double vapp(double t)
  {
    if(t<td)
      return Vlo;
    else
    {
      t-=td;
      while(t>pr) t-=pr;
      if(t<tr)
        return Vlo+t*(Vhi-Vlo)/tr;
      else if(t<tr+pw)
        return Vhi;
      else if(t<tr+pw+tf)
        return Vhi-(t-tr-pw)*(Vhi-Vlo)/tf;
      else    return Vlo;
    }

  }
};


class VEXP: public VSource
{
private:
  double td;
  double trc;
  double tfd;
  double tfc;
  double Vlo,Vhi;
public:
  VEXP(double t1,double v1,double v2,double t2,double t3,double t4):
  td(t1),trc(t2),tfd(t3),tfc(t4),Vlo(v1),Vhi(v2){};
  double vapp(double t)
  {
    if(t<=td)
      return Vlo;
    else if(t<=tfd)
      return Vlo+(Vhi-Vlo)*(1-exp(-(t-td)/trc));
    else
      return Vlo+(Vhi-Vlo)*(1-exp(-(t-td)/trc))+(Vlo-Vhi)*(1-exp(-(t-tfd)/tfc));

  }
};

class VSHELL: public VSource
{
private:
  void   *dll;
  double (*Vapp_Shell)(double);
  double scale_t;
  double scale_V;
public:
  VSHELL(void * dp, void * fp, double s_t,double s_V)
  {
     dll = dp;
     Vapp_Shell = (double (*)(double)) fp;
     scale_t = s_t;
     scale_V = s_V;
  }
  double vapp(double t)
  {
     return scale_V*Vapp_Shell(t/scale_t);
  }
  ~VSHELL()  {dlclose(dll);}
};

#endif
