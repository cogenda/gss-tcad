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
/*  Last update: Feb 26, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#ifndef _bzonedata_h_
#define _bzonedata_h_
#include "typedef.h"
#include "mesh.h"
#include "material.h"


//----------------------------------------------------------------
class BZoneData
{
public:
  int           zone_index;
  ZONE          *pzone;
  int           material_type;
  int           material;
  int           node_num;
  int           tri_num; 
  virtual   int Init(ZONE*,double,PhysicalUnitScale *)=0;
  BZoneData()   { pzone = 0; }
  virtual  ~BZoneData() {}
};

class SMCZone;
class ISZone;
class ElZone;
class VacuumZone;//```
class PMLZone;


#endif
