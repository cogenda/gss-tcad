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
/*  Last update: March 29, 2006                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _doping_h_
#define _doping_h_

#include "mathfunc.h"
#include "typedef.h"
const int  Uniform  = 1011;
const int  Gauss    = 1012;
const int  Erf      = 1013;

const double DONOR    = 1.0;
const double ACCEPTOR = -1.0;

enum PROFILE_CARD_ERROR
{
  PROFILE_NO_ERROR,
  PROFILE_UNKNOW_PARAMETER,
  PROFILE_ION_TYPE,
  PROFILE_XLOCATION,
  PROFILE_YLOCATION,
  PROFILE_NEGTIVE_DOSE,
  PROFILE_NEGTIVE_CONCENTR,
  PROFILE_DOPING_CHAR_LENGTH
};


class DopingFunc
{
public:
  double ion;     // N-ion or P-ion
  double xmin;    // bound box
  double xmax;
  double ymin;
  double ymax;
public:
  virtual double profile(double x,double y)=0;
};

class UniformDopingFunc : public DopingFunc
{
public:
  double PEAK;
  UniformDopingFunc(double xleft, double xright, double ytop, double ybottom, double doping_ion, double N)
  {
    PEAK = N;
    xmin = xleft;
    xmax = xright;
    ymin = ybottom;
    ymax = ytop;
    ion = doping_ion;
  }
  double profile(double x,double y)
  {
    if(x>=xmin-1e-6 && x<=xmax+1e-6 && y>=ymin-1e-6 &&y<=ymax+1e-6)
      return ion*PEAK;
    else
      return 0;
  }
};

class GaussDopingFunc : public DopingFunc
{
public:
  double PEAK;
  double XCHAR;
  double YCHAR;
  GaussDopingFunc(double xleft, double xright, double ytop, double ybottom,
                 double doping_ion, double N, double ax, double ay)
  {
    PEAK = N;
    XCHAR = ax;
    YCHAR = ay;
    xmin = xleft;
    xmax = xright;
    ymin = ybottom;
    ymax = ytop;
    ion = doping_ion;
  }
  double profile(double x,double y)
  {
    double dx,dy;
    if(x<xmin)
      dx = exp(-(x-xmin)*(x-xmin)/(XCHAR*XCHAR));
    else if(x>=xmin&&x<=xmax)
      dx = 1.0;
    else
      dx = exp(-(x-xmax)*(x-xmax)/(XCHAR*XCHAR));

    if(y<ymin)
      dy = exp(-(y-ymin)*(y-ymin)/(YCHAR*YCHAR));
    else if(y>=ymin&&y<=ymax)
      dy = 1.0;
    else
      dy = exp(-(y-ymax)*(y-ymin)/(YCHAR*YCHAR));

    return ion*PEAK*dx*dy;
  }
};

class ErfDopingFunc : public DopingFunc
{
public:
  double PEAK;
  double XCHAR;
  double YCHAR;
  ErfDopingFunc(double xleft, double xright, double ytop, double ybottom,
                 double doping_ion, double N, double ax, double ay)
  {
    PEAK = N;
    XCHAR = ax;
    YCHAR = ay;
    xmin = xleft;
    xmax = xright;
    ymin = ybottom;
    ymax = ytop;
    ion = doping_ion;
  }
  double profile(double x,double y)
  {
    double dx,dy;
    dx = (erfc((x-xmax)/XCHAR)-erfc((x-xmin)/XCHAR))/2;
    if(y<ymin)
      dy = exp(-(y-ymin)*(y-ymin)/(YCHAR*YCHAR));
    else if(y>=ymin&&y<=ymax)
      dy = 1.0;
    else
      dy = exp(-(y-ymax)*(y-ymin)/(YCHAR*YCHAR));
    return ion*PEAK*dx*dy;
  }
};
#endif
