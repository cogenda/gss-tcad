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
/*  Last update: April 19, 2006                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/


#ifndef _skeleton_h_
#define _skeleton_h_
#include <vector>
#include <map>
#include <math.h>
#define ANSI_DECLARATORS

extern "C"
{
#include "triangle.h"
}

using namespace std;
const int INNER  = 10000;
const int TOP    = 10001;
const int BOTTOM = 10002;
const int LEFT   = 10003;
const int RIGHT  = 10004;


const int XDIR   = 10107;
const int YDIR   = 10108;

const int LabeledSegment = 10201;
const int Interface      =-10202;
const int NoLabel        = 10203;

const int Triangle  = 3;
const int Rectangle = 4;
const int Hexagon   = 6;
const int Ellipse   = 0;


class skeleton_point
{
public:
  double x,y;    //location
  int    index;
  int    eliminated;
  void   set_location(double a,double b) {x=a;y=b;}
  skeleton_point():x(0),y(0),eliminated(0) {}
  skeleton_point(double a,double b):x(a),y(b),eliminated(0) {}
};

class aux_point
{
public:
 double x,y;
 int    index;
 aux_point():x(0),y(0) {}
 aux_point(double a,double b):x(a),y(b) {}
};

class aux_edge
{
public:
 int    segment_mark;
 int    p1,p2; 
};

class skeleton_edge
{
public:
  int    IX,IY;
  int    p1[2];    //p1[index_x][index_y]
  int    p2[2];    //p2[index_x][index_y]
};


struct lt_edge
{
  bool operator()(skeleton_edge e1, skeleton_edge e2) const
  {
    return (e1.p1[1]*e1.IX*e1.IY+e1.p1[0]+e1.p2[1]*e1.IX*e1.IY+e1.p2[0] <
            e2.p1[1]*e2.IX*e2.IY+e2.p1[0]+e2.p2[1]*e2.IX*e2.IY+e2.p2[0]);
  }
};


class skeleton_segment
{
public:
  vector<skeleton_edge>  edge_list; //input
  int    segment_mark;
  char   segment_label[32];
};


enum   mole_grad{MOLE_GRAD_X,MOLE_GRAD_Y};


class skeleton_region
{
public:
  char   label[32];
  char   material[32];
  int    shape;
  vector<double> px;             //half point location within the region
  vector<double> py;             //specify region material  
  
  //for regtangle 
  int    ixmin,ixmax,iymin,iymax;  //bound box
  double xmin,xmax,ymin,ymax;      //bound box
  //for ellipse
  double centrex,centrey;          //the centre of the ellipse
  double major_radii,minor_radii;  //major and minor radii 
  double theta;                    //the rotary angle
  int    division;                 //
  
  double mole_x1;              // for compound materials
  double mole_x1_slope;
  mole_grad mole_x1_grad;
  
  int    node_num;             //output
  int    tri_num;              //output
  vector<int> boundary;        //output
};


class skeleton_line
{
public:
  vector<skeleton_point>  point_list;
  void insert(skeleton_point &a)       { point_list.push_back(a); }
};

class output_edge
{
public:
  int index;
  int p1,p2;
  int interface;
  int r1,r2;
  int bc_type;
  int mark;
};


class mesh_constructor
{
public:
  int  point_num;
  int  aux_point_num;
  skeleton_point **point_array2d;
  vector<aux_point>        aux_point_array1d;
  map<skeleton_edge, int, lt_edge> edge_table;
  vector<aux_edge>         aux_edge_array1d; 
  vector<skeleton_segment> segment_array1d;
  vector<skeleton_region>  region_array1d;
  int  IX,IY;
  double xmin,xmax,ymin,ymax;
  double width,depth;
  int  find_skeleton_line_x(double x);
  int  find_skeleton_line_y(double y);
  int  find_move_skeleton_line_x(double x);
  int  find_move_skeleton_line_y(double y);
  int  x_eliminate(int ixmin,int ixmax, int iymin, int iymax);
  int  y_eliminate(int iymin,int iymax, int ixmin, int ixmax);
  void make_node_index();
  int  set_region_rectangle();
  int  set_segment_rectangle(int,int,int,int,char *);
  int  set_segment_rectangle(int,int,int,int,int,char *);
  int  set_spread_rectangle(int,double,int,int,double,double,double,double);
  int  set_region_ellipse();
public:
  vector<output_edge> out_edge;
  struct triangulateio in,out;
  int  do_mesh(char *tri_arg);
  int  to_cgns(const char *filename);
  mesh_constructor(skeleton_line &lx,skeleton_line &ly);
  ~mesh_constructor();
};


#endif
