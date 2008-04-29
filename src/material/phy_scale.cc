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
/*  Last update: Oct 17 2005                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#include "phy_scale.h"


// default scale value
PhysicalUnitScale::PhysicalUnitScale()
{
  s_centimeter = 1e6;
  s_second = 1e12;
  s_volt = 1.0;
  s_coulomb = 1.0/1.602176462e-19;
  s_kelvin = 1.0/300;

  s_meter = s_centimeter*1e2;
  s_micron = s_centimeter/1e4;
  s_joule = s_coulomb*s_volt;
  s_W  = s_joule/s_second;
  s_kg = s_joule/(s_meter*s_meter)*s_second*s_second;
  s_eV = s_joule*1.602176462e-19;
  s_ps = s_second/1e12;
  s_A  = s_coulomb/s_second;
  s_mA  = 1e-3*s_coulomb/s_second;
}

//user defined scale value : length only
void PhysicalUnitScale::SetPhysicalUnitScale(double s_length)
{

  s_centimeter = s_length;
  s_second = 1e12;
  s_volt = 1.0;
  s_coulomb = 1.0/1.602176462e-19;
  s_kelvin = 1.0/300;
  
  s_meter = s_centimeter*1e2;
  s_micron = s_centimeter/1e4;
  s_joule = s_coulomb*s_volt;
  s_W  = s_joule/s_second;
  s_kg = s_joule/(s_meter*s_meter)*s_second*s_second;
  s_eV = s_joule*1.602176462e-19;
  s_ps = s_second/1e12;
  s_A  = s_coulomb/s_second;
  s_mA  = 1e-3*s_coulomb/s_second;
}

//user defined scale value : length and potential
void PhysicalUnitScale::SetPhysicalUnitScale(double s_length,double s_potential)
{
  s_centimeter = s_length;
  s_second = 1e12;
  s_volt = s_potential;
  s_coulomb = 1.0/1.602176462e-19;
  s_kelvin = 1.0/300;

  s_meter = s_centimeter*1e2;
  s_micron = s_centimeter/1e4;
  s_joule = s_coulomb*s_volt;
  s_W  = s_joule/s_second;
  s_kg = s_joule/(s_meter*s_meter)*s_second*s_second;
  s_eV = s_joule*1.602176462e-19;
  s_ps = s_second/1e12;
  s_A  = s_coulomb/s_second;
  s_mA  = 1e-3*s_coulomb/s_second;
}


//user defined scale value
void PhysicalUnitScale::SetPhysicalUnitScale(double s_length, double s_time,double s_voltage,
                                             double s_charge, double s_temperature)
{
  s_centimeter = s_length;
  s_second = s_time;
  s_volt   = s_voltage;
  s_coulomb= s_charge;
  s_kelvin = s_temperature;

  s_meter = s_centimeter*1e2;
  s_micron = s_centimeter/1e4;
  s_joule = s_coulomb*s_volt;
  s_W  = s_joule/s_second;
  s_kg = s_joule/(s_meter*s_meter)*s_second*s_second;
  s_eV = s_joule*1.602176462e-19;
  s_ps = s_second/1e12;
  s_A  = s_coulomb/s_second;
  s_mA  = 1e-3*s_coulomb/s_second;
}

