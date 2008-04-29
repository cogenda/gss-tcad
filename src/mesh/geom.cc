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
/*  Last update: May 29, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "geom.h"


/* ----------------------------------------------------------------------------
 * CompAngle:  This function computes the angle between three points given
 * by coordinates (x1,y1), (x2,y2) and (x3,y3).  The angle is measured from
 * vector 12 to 32.clock order
 */
double CompAngle(double x1, double y1,double x2, double y2,double x3, double y3)
{
  double xl,yl;                         /* coordinates of vector 12 */
  double xr,yr;                         /* coordinates of vector 32 */
  double r;                             /* temporary number */
  double angle;                         /* angle between vectors 12 & 32 */
  double PI=3.14159265358979323846;
  xl = x1 - x2;  yl = y1 - y2;
  xr = x3 - x2;  yr = y3 - y2;

  r = (xl*xr+yl*yr)/(sqrt(xl*xl+yl*yl)*sqrt(xr*xr+yr*yr));
  if (r > 1.0)
    r =  1.0;
  else if (r < -1.0) r = -1.0;
  angle = acos(r);
  if (xl*yr-xr*yl > 0)
    angle = 2*PI-angle;
  return(angle);

} /* CompAngle */

/* ----------------------------------------------------------------------------
 * CircleCenter:  This function determines the center of the circle that
 * passes through points (x1,y1), (x2,y2) and (x3,y3).  If the points are
 * coplaner then 1 is returned otherwise 0 is
 * returned.
 */
int CircleCenter(double x1, double y1, double x2, double y2,double x3, double y3,
                  double &xc, double &yc)
{
  double xl,yl;                         /* coordinates of vector 12 */
  double xr,yr;                         /* coordinates of vector 32 */
  double rl,rr;                         /* ??? */
  double det;                           /* determinent of coordinate matrix */

  xl = x1 - x2;  yl = y1 - y2;
  xr = x3 - x2;  yr = y3 - y2;

  det = xr * yl - xl * yr;
  if (det == 0.0)
  {
    xc = yc = 0.0;
    return 1;
  }

  rl = (xl * (x1 + x2) + yl * (y1 + y2)) / det;
  rr = (xr * (x3 + x2) + yr * (y3 + y2)) / det;

  xc = -0.5 * (rl * yr - rr * yl);
  yc = -0.5 * (rr * xl - rl * xr);
  return 0;
} /* CircleCenter */


/* ----------------------------------------------------------------------------
 * MirrorPoint:  This function computes the MirrorPoint(xg,yg) of (x,y) with the
 * line determined by (x1,y1) and (x2,y2).
 */
void MirrorPoint(double x,double y,double x1,double y1,double x2,double y2,double &xg,double &yg)
{
  double r,ang1,ang2,angle;


  double r1 = sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y));
  double r2 = sqrt((x2-x)*(x2-x)+(y2-y)*(y2-y));
  if(r1>r2)
  {
    r = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    ang1 = CompAngle(x2,y2,x1,y1,x,y);
    ang2 = CompAngle(x2,y2,x1,y1,x1+r,y1);
    angle = ang1 + ang2;
    xg = x1 + r1*cos(angle);
    yg = y1 + r1*sin(angle);
  }
  else
  {
    r = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    ang1 = CompAngle(x1,y1,x2,y2,x,y);
    ang2 = CompAngle(x1,y1,x2,y2,x2+r,y2);
    angle = ang1 + ang2;
    xg = x2 + r2*cos(angle);
    yg = y2 + r2*sin(angle);

  }
}/* CompAngle */


/* ----------------------------------------------------------------------------
 * DistancePointLine:  This function computes the distance between point (x,y) and
 * line determined by (x1,y1) and (x2,y2).
 */
double DistancePointLine(double x1, double y1,double x2, double y2,double x, double y)
{
  double LineMag = Distance(x1,y1,x2,y2);
  double U = (((x-x1)*(x2-x1))+((y-y1)*(y2-y1))) / (LineMag*LineMag);

  double IntersectionX = x1 + U*(x2-x1);
  double IntersectionY = y1 + U*(y2-y1);

  return  Distance( x, y, IntersectionX, IntersectionY);
}/* DistancePointLine */


/* ----------------------------------------------------------------------------
 * Distance:  This function computes the distance between (x1,y1) and (x2,y2)
 * line determined by (x1,y1) and (x2,y2).
 */
double Distance(double x1, double y1,double x2,double y2)
{
  if(x1==x2 && y1==y2) return 0;
  return sqrt((x1-x2)*(x1-x2) +(y1-y2)*(y1-y2));
}/* Distance */


/* ----------------------------------------------------------------------------
 * TriArea:  This function computes the area of a triangle determined by (x1,y1)
 * (x2,y2) and (x3,y3).
 */
double TriArea(double x1, double y1,double x2, double y2,double x3, double y3)
{
  double d1,d2,d3,d,r,xc,yc;
  d1=Distance(x1,y1,x2,y2);
  d2=Distance(x1,y1,x3,y3);
  d3=Distance(x2,y2,x3,y3);
  d=0.5*(d1+d2+d3);
  return sqrt(d*(d-d1)*(d-d2)*(d-d3));
}/* TriArea */


/* ----------------------------------------------------------------------------
 * TriAngle:  This function compute the max angle of a triangle
 *
 */
double MaxTriAngle(double x1, double y1, double x2, double y2,double x3, double y3)
{
  double a,b,c;
  a=Distance(x1,y1,x2,y2);
  b=Distance(x1,y1,x3,y3);
  c=Distance(x2,y2,x3,y3);
  if(a>b&&a>c) return acos((b*b+c*c-a*a)/(2*b*c));
  if(b>a&&b>c) return acos((a*a+c*c-b*b)/(2*a*c));
  if(c>b&&c>a) return acos((a*a+b*b-c*c)/(2*a*b));
  return 0; //prevent warning of complier
} /* TriAngle */


/* ----------------------------------------------------------------------------
 * RaySegmentIntersectTest:  This function do 2D Ray-Segment Intersection test
 *
 */
IntersectResult RaySegmentIntersectTest(double begin_x, double begin_y, double ex, double ey, double x1,double y1,double x2,double y2)
{
  double end_x =  begin_x + ex/sqrt(ex*ex+ey*ey+1e-20);
  double end_y =  begin_y + ey/sqrt(ex*ex+ey*ey+1e-20);
  double denom  = ((y2 - y1)*(end_x - begin_x)) -((x2 - x1)*(end_y - begin_y));
  double nume_a = ((x2 - x1)*(begin_y - y1)) - ((y2 - y1)*(begin_x - x1));
  double nume_b = ((end_x - begin_x)*(begin_y - y1)) - ((end_y - begin_y)*(begin_x - x1));
  if(denom == 0.0)
  {
    if(nume_a == 0.0 && nume_b == 0.0)
    {
      return COINCIDENT;
    }
    return PARALLEL;
  }

  double ua = nume_a / denom;
  double ub = nume_b / denom;

  if(ub >= 0.0 && ub <= 1.0)
  {
    // Get the intersection point.
    double intersection_x = begin_x + ua*(end_x - begin_x);
    double intersection_y = begin_y + ua*(end_y - begin_y);
    return INTERESECTING;
  }

  return NOT_INTERESECTING;
}

