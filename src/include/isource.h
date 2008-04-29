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
/*  Last update: Dec 21, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
//sevel types of isource are defined here
//reference: spice vsource model

#ifndef _isource_h_
#define _isource_h_
#include <math.h>
#include <dlfcn.h>
#include "typedef.h"

enum ISOURCE_CARD_ERROR
{
  ISOURCE_NO_ERROR,
  ISOURCE_NO_ID,
  ISOURCE_UNKNOW_PARAMETER
};

//virtual base class
class ISource
{
public:
  char label[32];
  virtual double iapp(double t)=0;
  virtual ~ISource() {}
};


class IDC: public ISource
{
private:
  double td;
  double Idc;
public:
  IDC(double t1,double I1): td(t1),Idc(I1){};
  double iapp(double t)
{ return t>=td? Idc:0;};
};



class ISIN: public ISource
{
private:
  double td;
  double Iamp;
  double fre;
public:
  ISIN(double t1,double I1,double f1): td(t1),Iamp(I1),fre(f1){};
  double iapp(double t)
  { return t>=td? Iamp*sin(2*3.14159265359*fre*(t-td)):0;};
};


class IPULSE: public ISource
{
private:
  double td;
  double tr;
  double tf;
  double pw;
  double pr;
  double Ilo,Ihi;
public:
  IPULSE(double t1,double I1,double I2,double t2,double t3,double t4, double t5):
  td(t1),tr(t2),tf(t3),pw(t4),pr(t5),Ilo(I1),Ihi(I2){};
  double iapp(double t)
  {
    if(t<td)
      return Ilo;
    else
    {
      t-=td;
      while(t>pr) t-=pr;
      if(t<tr)
        return Ilo+t*(Ihi-Ilo)/tr;
      else if(t<tr+pw)
        return Ihi;
      else if(t<tr+pw+tf)
        return Ihi-(t-tr-pw)*(Ihi-Ilo)/tf;
      else    return Ilo;
    }

  }
};


class IEXP: public ISource
{
private:
  double td;
  double trc;
  double tfd;
  double tfc;
  double Ilo,Ihi;
public:
  IEXP(double t1,double I1,double I2,double t2,double t3,double t4):
  td(t1),trc(t2),tfd(t3),tfc(t4),Ilo(I1),Ihi(I2){};
  double iapp(double t)
  {
    if(t<=td)
      return Ilo;
    else if(t<=tfd)
      return Ilo+(Ihi-Ilo)*(1-exp(-(t-td)/trc));
    else
      return Ilo+(Ihi-Ilo)*(1-exp(-(t-td)/trc))+(Ilo-Ihi)*(1-exp(-(t-tfd)/tfc));

  }
  ;
};

class ISHELL: public ISource
{
private:
  void   *dll;
  double (*Iapp_Shell)(double);
  double scale_t;
  double scale_mA;
public:
  ISHELL(void * dp, void * fp, double s_t,double s_mA)
  {
     dll = dp;
     Iapp_Shell = (double (*)(double)) fp;
     scale_t = s_t;
     scale_mA = s_mA;
  }
  double iapp(double t)
  {
     return scale_mA*Iapp_Shell(t/scale_t);
  }
  ~ISHELL()  {dlclose(dll);}
};
#endif
