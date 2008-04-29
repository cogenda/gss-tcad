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
/*  Last update: Oct 17 2005                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#ifndef _phy_scale_h_
#define _phy_scale_h_

#include "typedef.h"

class PhysicalUnitScale
{
public:
  // use as fundamental physical unit
  double   s_meter;	 // the length unit
  double   s_second;	 // the time unit
  double   s_volt;	 // potential unit
  double   s_coulomb;	 // the charge unit
  double   s_kelvin;	 // the temperature unit
public:
  // set as induced physical unit
  double   s_micron;
  double   s_centimeter;
  double   s_kg;		      // the mass unit
  double   s_joule,s_eV;	      // energy unit
  double   s_W;                       // power unit 
  double   s_ps;
  double   s_A;
  double   s_mA;
public:
  PhysicalUnitScale();
  void SetPhysicalUnitScale(double);
  void SetPhysicalUnitScale(double,double);
  void SetPhysicalUnitScale(double,double,double,double,double);
};

#endif
