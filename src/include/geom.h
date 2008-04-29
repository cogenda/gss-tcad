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

#ifndef _geom_h_
#define _geom_h_
#include "typedef.h"

extern double CompAngle(double , double ,double , double ,double , double );
extern int    CircleCenter(double , double , double , double ,double , double ,double &, double &);
extern void   MirrorPoint(double ,double ,double ,double ,double ,double ,double &,double &);
extern double DistancePointLine(double , double ,double , double ,double , double );
extern double Distance(double , double ,double ,double );
extern double TriArea(double, double, double, double, double, double );
extern double MaxTriAngle(double, double, double, double, double, double );
enum IntersectResult { PARALLEL, COINCIDENT, NOT_INTERESECTING, INTERESECTING };
IntersectResult RaySegmentIntersectTest(double, double, double, double, double,double,double,double);
#endif
