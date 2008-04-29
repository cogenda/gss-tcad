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
/*  Last update: Apr 13, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _mesh_h_
#define _mesh_h_

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include "element.h"
#include "interface.h"
#include "zone.h"
#include "cmdbuf.h"

using namespace std;

class MESH
{
protected:

  char        basename[32];               //cgns base name
  double      mesh_scale;                 //this value will be mul to coordinate
  int         bc_num;                     //total bc num  summed by each zone num
  vector<GNode>    gnode;                 //these are global mesh structure
  vector<Tri>      gtri;
  vector<GSegment> gsegment;              //used for mesh refine
public:
  int         zone_num;
  ZONE        *zone;                    
  ZoneInterface zinterface;

protected:
  //variables and functions for drawing mesh
  double      xMin;                     // minimum value of x coordinate
  double      xMax;                     // maximum value of x coordinate
  double      yMin;                     // minimum value of y coordinate
  double      yMax;                     // maximum value of y coordinate
protected:
  ZONE*       Get_zone(const char *);
  ZONE*       Get_zone(const int);
  int         Get_zone_index(const char *);
  Segment*    Get_segment(const int);
  int*        reorder_zone(int);              //needed only by solver.
public:
  int         import_cgns(const char *cgnsfile);
  void        export_mesh(const char *cgnsfile);
  void        build_zone();
  int         zone_to_field();
  int         field_to_zone();
  void        reorder_mesh();
  void        clear_mesh();
public:

  MESH        ():mesh_scale(1e6),bc_num(0),zone_num(0),zone(0)
  {xMin = +1.0e100;   yMin = +1.0e100;   xMax = -1.0e100;   yMax = -1.0e100;}
  ~MESH       () {delete [] zone;}
  void        set_mesh_scale(double scale) {mesh_scale = scale;}

};

#endif
