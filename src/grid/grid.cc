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
/*  Last update: April 14, 2006                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bsolver.h"
#include "skeleton.h"
#include "log.h"
#include "typedef.h"

/* ----------------------------------------------------------------------------
 * set_xmesh:  This function check and do XMESH card
 * which init mesh line in x direction.
 */
int set_xmesh(list<Cmd>::iterator pcmd, skeleton_line & skeleton_line_x)
{
  // check parameters.
  if(!pcmd->allowed_args(9,"x.min","x.left","x.max","x.right",
                           "width","ratio","n.spaces","h1","h2"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " XMESH: unrecognized parameter(s)!\n";
        gss_log.record();
        return 1;
  }
  // get parameter value from command structure
  double xmax = pcmd->get_number("x.max",1,"x.right",0.0);
  double xmin = pcmd->get_number("x.min",1,"x.left",0.0);
  double width = pcmd->get_number("width",0,0.0);
  double x,dx;
  if(!skeleton_line_x.point_list.size())  // if no previous X.MESH card, must insert the first point.
  {
      skeleton_point  p(xmin,0);
      skeleton_line_x.insert(p);
  }
  else  // xmin is determined by previous X.MESH card.
    xmin = skeleton_line_x.point_list[skeleton_line_x.point_list.size()-1].x;

  if(width>0)   xmax = xmin+width;
  else          width = xmax-xmin;

  double ratio = pcmd->get_number("ratio",0,1.0) ;
  int  nspaces = pcmd->get_integer("n.spaces",0,1);
  
  if(pcmd->is_arg_exist("h1") && pcmd->is_arg_exist("h2"))
  {
        double h1 = pcmd->get_number("h1",0,0.1);
        double h2 = pcmd->get_number("h2",0,h1);
        if(h1==h2)
        {
                ratio = 1.0;
                nspaces = int(width/h1+0.5);
        }
        else
        {
                ratio = (width-h1)/(width-h2);
                nspaces = (int)(log(h2/h1)/log(ratio)+ 1 +0.5);
        }
  }

  else if(pcmd->is_arg_exist("h1"))
  {
        double h1 = pcmd->get_number("h1",0,0.1);
        ratio = pcmd->get_number("ratio",0,1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(width/h1+0.5);
        }
        else
        {
		nspaces = (int)(log(1-(1-ratio)*width/h1)/log(ratio)+0.5);
        }
  }
  
  else if(pcmd->is_arg_exist("h2"))
  {
        double h2 = pcmd->get_number("h2",0,0.1);
        ratio = 1.0/pcmd->get_number("ratio",0,1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(width/h2+0.5);
        }
        else
        {
		double h1 = width*(1-ratio)+h2*ratio;
		nspaces = int(log(h2/h1)/log(ratio)+1.0+0.5);  
        }
  }
  
  x = xmin;

  if(nspaces==1)
  {
      skeleton_point  p(xmax,0);
      skeleton_line_x.insert(p);
      return 0;
  }

  if(ratio==1.0)
    dx = width/nspaces;
  else
    dx = width*(ratio-1)/(std::pow(ratio,nspaces)-1);

  for(int i=1;i<=nspaces;i++)
  {
      x += dx;
      skeleton_point  p(x,0);
      skeleton_line_x.insert(p);
      dx*=ratio;
  }

  return 0;
}

/* ----------------------------------------------------------------------------
 * set_ymesh:  This function check and do YMESH card
 * which init mesh line in y direction.
 */
int set_ymesh(list<Cmd>::iterator pcmd, skeleton_line & skeleton_line_y)
{
  // check parameters.
  if(!pcmd->allowed_args(9,"y.min","y.bottom","y.max","y.top",
                           "depth","ratio","n.spaces","h1","h2"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " YMESH: unrecognized parameter(s)!\n";
        gss_log.record();
        return 1;
  }
  // get parameter value from command structure
  double ymax = pcmd->get_number("y.max",1,"y.top",0.0);
  double ymin = pcmd->get_number("y.min",1,"y.bottom",0.0);
  double depth = pcmd->get_number("depth",0,0.0);
  double y,dy;

  if(!skeleton_line_y.point_list.size()) // if no previous Y.MESH card, must insert the first point.
  {
      skeleton_point  p(0,ymax);
      skeleton_line_y.insert(p);
  }
  else  // ymax is determined by previous Y.MESH card.
    ymax = skeleton_line_y.point_list[skeleton_line_y.point_list.size()-1].y;

  if(depth>0)  ymin = ymax-depth;
  else         depth = ymax- ymin;

  double ratio = pcmd->get_number("ratio",0,1.0) ;
  int  nspaces = pcmd->get_integer("n.spaces",0,1);

  if(pcmd->is_arg_exist("h1") && pcmd->is_arg_exist("h2"))
  {
        double h1 = pcmd->get_number("h1",0,0.1);
        double h2 = pcmd->get_number("h2",0,h1);
        if(h1==h2)
        {
                ratio = 1.0;
                nspaces = int(depth/h1+0.5);
        }
        else
        {
                ratio = (depth-h1)/(depth-h2);
                nspaces = (int)(log(h2/h1)/log(ratio)+ 1 +0.5);
        }
  }
  
  else if(pcmd->is_arg_exist("h1") )
  {
        double h1 = pcmd->get_number("h1",0,0.1);
        ratio = pcmd->get_number("ratio",0,1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(depth/h1+0.5);
        }
        else
        {
		nspaces = (int)(log(1-(1-ratio)*depth/h1)/log(ratio)+0.5);
        }
  }
  
  else if(pcmd->is_arg_exist("h2") )
  {
        double h2 = pcmd->get_number("h2",0,0.1);
        ratio = 1.0/pcmd->get_number("ratio",0,1.0) ;
        if(ratio==1.0)
        {
                nspaces = int(depth/h2+0.5);
        }
        else
        {
		double h1 = depth*(1-ratio)+h2*ratio;
		nspaces = int(log(h2/h1)/log(ratio)+1.0+0.5);  
        }
  }

  y = ymax;

  if(nspaces==1)
  {
      skeleton_point  p(0,ymin);
      skeleton_line_y.insert(p);
      return 0;
  }

  if(ratio==1.0)
    dy = depth/nspaces;
  else
    dy = depth*(ratio-1)/(std::pow(ratio,nspaces)-1);

  for(int i=1;i<=nspaces;i++)
  {
      y -= dy;
      skeleton_point  p(0,y);
      skeleton_line_y.insert(p);
      dy*=ratio;
  }
  return 0;
}

/* ----------------------------------------------------------------------------
 * set_eliminate:  This function check and do ELIMINATE card
 * eliminate unnecessary lines in x or y direction.
 */
int set_eliminate(list<Cmd>::iterator pcmd, mesh_constructor & orig_mesh)
{

  // check parameters.
  if(!pcmd->allowed_args(17,"x.min","x.left","x.max","x.right",
                            "ix.min","ix.left","ix.max","ix.right",
                            "y.min","y.bottom","y.max","y.top",
                            "iy.min","iy.bottom","iy.max","iy.top",
                            "direction"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " ELIMINATE: unrecognized parameter(s)!\n";
        gss_log.record();
        return 1;
  }
  // get parameter value from command structure
  int ixmin = pcmd->get_integer("ix.min",1,"ix.left",0);
  int ixmax = pcmd->get_integer("ix.max",1,"ix.right",orig_mesh.IX-1);
  int iymin = pcmd->get_integer("iy.min",1,"iy.top",0);
  int iymax = pcmd->get_integer("iy.max",1,"iy.bottom",orig_mesh.IY-1);
  char *dir = pcmd->get_string("direction",0,"");
  if(pcmd->is_arg_exist("x.max") || pcmd->is_arg_exist("x.right"))
  {
        double xmax = pcmd->get_number("x.max",1,"x.right",orig_mesh.xmax);
        ixmax = orig_mesh.find_skeleton_line_x(xmax);
  }
  if(pcmd->is_arg_exist("x.min") || pcmd->is_arg_exist("x.left"))
  {
        double xmin = pcmd->get_number("x.min",1,"x.left",orig_mesh.xmin);
        ixmin = orig_mesh.find_skeleton_line_x(xmin);
  }
  if(pcmd->is_arg_exist("y.max") || pcmd->is_arg_exist("y.top"))
  {
        double ymax = pcmd->get_number("y.max",1,"y.top",orig_mesh.ymax);
        iymin = orig_mesh.find_skeleton_line_y(ymax);
  }
  if(pcmd->is_arg_exist("y.min") || pcmd->is_arg_exist("y.bottom"))
  {
        double ymin = pcmd->get_number("y.min",1,"y.bottom",orig_mesh.ymin);
        iymax = orig_mesh.find_skeleton_line_y(ymin);
  }
  if(!strcmp(dir,"COLUMNS"))
        orig_mesh.y_eliminate(iymin,iymax,ixmin,ixmax);
  else if(!strcmp(dir,"ROWS"))
        orig_mesh.x_eliminate(ixmin,ixmax,iymin,iymax);
  else
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " ELIMINATE: wrong eliminate direction!\n";
        gss_log.record();
        return 1;
  }
  return 0;
}




/* ----------------------------------------------------------------------------
 * set_spread:  This function check and do SPREAD card
 * I copyed some code from PISCES. please forgive me...
 */
int set_spread(list<Cmd>::iterator pcmd, mesh_constructor & orig_mesh)
{
  // check parameters.
  if(!pcmd->allowed_args(14,"location","width","upper","lower","encroach",
                            "y.lower","fix.lower","thickness","vol.rat",
                            "grading","middle","y.middle","gr1","gr2"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " SPREAD: unrecognized parameter(s)!\n";
        gss_log.record();
        return 1;
  }
  // get parameter value from command structure

  // get x location of distorted region.
  int location;
  int xpoint;
  if(!pcmd->is_arg_exist("location"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " SPREAD: you must give location of spread region!\n";
        gss_log.record();
        return 1;
  }
  if(pcmd->is_arg_value("location","Left"))
    {location = LEFT;xpoint=0;}
  else if(pcmd->is_arg_value("location","Right"))
    {location = RIGHT;xpoint=orig_mesh.IX-1;}
  else
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " SPREAD: you must give location of spread region!\n";
        gss_log.record();
        return 1;
  }
  double width = pcmd->get_number("width",0,0.0);

  // the upper and lower grid line of distorted region
  int  upperline = pcmd->get_integer("upper",0,0);
  int  lowerline = pcmd->get_integer("lower",0,0);
  // current thickness
  double cthick = orig_mesh.point_array2d[upperline][xpoint].y-orig_mesh.point_array2d[lowerline][xpoint].y;

  // get the new location of upper and lower y grid line
  double yuploc,yloloc,thick;
  if(pcmd->is_arg_exist("y.lower"))
  {
     yuploc = orig_mesh.point_array2d[upperline][xpoint].y;
     yloloc = pcmd->get_number("y.lower",0,0.0);
     thick = yuploc - yloloc;
  }
  else if(pcmd->is_arg_exist("thickness"))
  {
    thick = pcmd->get_number("thickness",0,1.0);
    double vol_rat = pcmd->get_number("vol.rat",0,0.44);
    double dthick = vol_rat*(thick-cthick);
    double uthick = thick-dthick-cthick;
    yuploc=orig_mesh.point_array2d[upperline][xpoint].y + uthick;
    yloloc=orig_mesh.point_array2d[lowerline][xpoint].y - dthick;
  }

  // this parameter give the transition between distorted and undistorted grid. in x direction.
  // from PISCES:
  //   scale so that 80% of thickness change occurs in a
  //   distance equal to the thickness change (for an
  //   encroachment factor of 1)
  double encroach =  pcmd->get_number("encroach",0,1.0);
  if(encroach<0.1) encroach = 0.1;
  const double erfc80 = 1.812386;
  encroach=fabs(erfc80/((thick-cthick)*encroach));

  // get grading parameter(s)
  double grading =  pcmd->get_number("grading",0,1.0);

  orig_mesh.set_spread_rectangle(location,width,upperline,lowerline,yuploc,yloloc,encroach,grading);
  return 0;
}


int set_region_ellipse(list<Cmd>::iterator pcmd, mesh_constructor & orig_mesh)
{
  // check parameters.
  if(!pcmd->allowed_args(9,"shape","centrex","centrey","majorradii","minorradii","theta","division","label","material"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " REGION: unrecognized parameter(s)!\n";
        gss_log.record();
        return 1;
  }
  // get parameter value from command structure
  skeleton_region region;
  region.shape = Ellipse;
  region.division = pcmd->get_integer("division",0,12);
  region.centrex = pcmd->get_number("centrex",0,0.0);
  region.centrey = pcmd->get_number("centrey",0,0.0);
  region.major_radii  = pcmd->get_number("majorradii",0,1.0);
  region.minor_radii  = pcmd->get_number("minorradii",0,region.major_radii);
  region.theta        = pcmd->get_number("theta",0,0.0)/360*(2*3.14159265359);
  if(region.major_radii<=0 || region.minor_radii<=0)
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " REGION: Circle radius should be positive!\n";
        gss_log.record();
        return 1;
  }
  strcpy( region.label, pcmd->get_string("label",0,"")) ;
  strcpy( region.material, pcmd->get_string("material",0,"")) ;
  orig_mesh.region_array1d.push_back(region);
  return 0;
}

 
int set_region_rectangle(list<Cmd>::iterator pcmd, mesh_constructor & orig_mesh)
{
  // check parameters.
  if(!pcmd->allowed_args(24,"x.min","x.left","x.max","x.right",
                            "ix.min","ix.left","ix.max","ix.right",
                            "y.min","y.bottom","y.max","y.top",
                            "iy.min","iy.bottom","iy.max","iy.top",
                            "label","material","x.mole","mole.begin",
                            "mole.slope","mole.end","mole.grad","shape"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " REGION: unrecognized parameter(s)!\n";
        gss_log.record();
        return 1;
  }
  // get parameter value from command structure
  skeleton_region region;
  region.shape = Rectangle;
  region.ixmin = pcmd->get_integer("ix.min",1,"ix.left",0);
  region.ixmax = pcmd->get_integer("ix.max",1,"ix.right",orig_mesh.IX-1);
  region.iymin = pcmd->get_integer("iy.min",1,"iy.top",0);
  region.iymax = pcmd->get_integer("iy.max",1,"iy.bottom",orig_mesh.IY-1);
  strcpy( region.label, pcmd->get_string("label",0,"")) ;
  strcpy( region.material, pcmd->get_string("material",0,"")) ;

  if(pcmd->is_arg_exist("x.max") || pcmd->is_arg_exist("x.right"))
  {
        double xmax = pcmd->get_number("x.max",1,"x.right",orig_mesh.xmax);
        region.ixmax = orig_mesh.find_move_skeleton_line_x(xmax);
  }
  if(pcmd->is_arg_exist("x.min") || pcmd->is_arg_exist("x.left"))
  {
        double xmin = pcmd->get_number("x.min",1,"x.left",orig_mesh.xmin);
        region.ixmin = orig_mesh.find_move_skeleton_line_x(xmin);
  }
  if(pcmd->is_arg_exist("y.max") || pcmd->is_arg_exist("y.top"))
  {
        double ymax = pcmd->get_number("y.max",1,"y.top",orig_mesh.ymax);
        region.iymin = orig_mesh.find_move_skeleton_line_y(ymax);
  }
  if(pcmd->is_arg_exist("y.min") || pcmd->is_arg_exist("y.bottom"))
  {
        double ymin = pcmd->get_number("y.min",1,"y.bottom",orig_mesh.ymin);
        region.iymax = orig_mesh.find_move_skeleton_line_y(ymin);
  }

  if(region.ixmin>=region.ixmax)
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " REGION: Can't locate left/right boundary of region!\n";
        gss_log.record();
        return 1;
  }
  if(region.iymin>=region.iymax)
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " REGION: Can't locate top/bottom boundary of region!\n";
        gss_log.record();
        return 1;
  }
  region.xmin = dmin( orig_mesh.point_array2d[region.iymin][region.ixmin].x,
                      orig_mesh.point_array2d[region.iymax][region.ixmin].x);
  region.ymax = dmax( orig_mesh.point_array2d[region.iymin][region.ixmin].y,
                      orig_mesh.point_array2d[region.iymin][region.ixmax].y);
  region.xmax = dmax( orig_mesh.point_array2d[region.iymin][region.ixmax].x,
                      orig_mesh.point_array2d[region.iymax][region.ixmax].x);
  region.ymin = dmin( orig_mesh.point_array2d[region.iymax][region.ixmin].y,
                      orig_mesh.point_array2d[region.iymax][region.ixmax].y);

  if(IsSingleCompSemiconductor(region.material))
  {
        region.mole_x1 = pcmd->get_number("x.mole",0,0.0);
        region.mole_x1_slope = pcmd->get_number("mole.slope",0,0.0);
        region.mole_x1_grad = MOLE_GRAD_Y;
        if(pcmd->is_arg_value("mole.grad","X.Linear"))
          region.mole_x1_grad = MOLE_GRAD_X;
        else if(pcmd->is_arg_value("mole.grad","Y.Linear"))
          region.mole_x1_grad = MOLE_GRAD_Y;
        if(pcmd->is_arg_exist("mole.end"))
        {
          double mole_x1_end = pcmd->get_number("mole.end",0,0.0);
          double length;
          if(region.mole_x1_grad == MOLE_GRAD_X) length = region.xmax - region.xmin;
          if(region.mole_x1_grad == MOLE_GRAD_Y) length = region.ymax - region.ymin;
          region.mole_x1_slope = (mole_x1_end - region.mole_x1)/length;
        }
  }

  orig_mesh.region_array1d.push_back(region);
  orig_mesh.point_array2d[region.iymax][region.ixmin].eliminated = 0;
  orig_mesh.point_array2d[region.iymax][region.ixmax].eliminated = 0;
  orig_mesh.point_array2d[region.iymin][region.ixmin].eliminated = 0;
  orig_mesh.point_array2d[region.iymin][region.ixmax].eliminated = 0;
  return 0;

}

/* ----------------------------------------------------------------------------
 * set_region:  This function check and do REGION card
 * specify the material of a region
 */
int set_region(list<Cmd>::iterator pcmd, mesh_constructor & orig_mesh) 
{
  int flag;
  string shape = pcmd->get_string("shape",0,"Rectangle");
  if(shape=="Rectangle") flag=set_region_rectangle(pcmd,orig_mesh);
  if(shape=="Ellipse")   flag=set_region_ellipse(pcmd,orig_mesh);
  return flag;
}


/* ----------------------------------------------------------------------------
 * set_segment:  This function check and do SEGMENT card
 * specify the lable and location of a segment.
 * these segments may be specified as electrode later by BOUNDARY card.
 */
int set_segment(list<Cmd>::iterator pcmd, mesh_constructor & orig_mesh)
{
  int ixmin,ixmax,iymin,iymax;
  int ix,iy;
  // check parameters.
  if(!pcmd->allowed_args(23,"x","x.min","x.left","x.max","x.right",
                            "ix","ix.min","ix.left","ix.max","ix.right",
                            "y","y.min","y.bottom","y.max","y.top",
                            "iy","iy.min","iy.bottom","iy.max","iy.top",
                            "label","location","direction"))
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " SEGMENT: unrecognized parameter(s)!\n";
        gss_log.record();
        return 1;
  }
  // get parameter value from command structure
  ix    = pcmd->get_integer("ix",0,0);
  ixmin = pcmd->get_integer("ix.min",1,"ix.left",0);
  ixmax = pcmd->get_integer("ix.max",1,"ix.right",orig_mesh.IX-1);
  iy    = pcmd->get_integer("iy",0,0);
  iymin = pcmd->get_integer("iy.min",1,"iy.top",0);
  iymax = pcmd->get_integer("iy.max",1,"iy.bottom",orig_mesh.IY-1);
  char *label = pcmd->get_string("label",0,"");
  if(pcmd->is_arg_exist("x"))
  {
        double x = pcmd->get_number("x",0,0.0);
        ix = orig_mesh.find_skeleton_line_x(x);
  }
  if(pcmd->is_arg_exist("x.max") || pcmd->is_arg_exist("x.right"))
  {
        double xmax = pcmd->get_number("x.max",1,"x.right",orig_mesh.xmax);
        ixmax = orig_mesh.find_skeleton_line_x(xmax);
  }
  if(pcmd->is_arg_exist("x.min") || pcmd->is_arg_exist("x.left"))
  {
        double xmin = pcmd->get_number("x.min",1,"x.left",orig_mesh.xmin);
        ixmin = orig_mesh.find_skeleton_line_x(xmin);
  }
  if(pcmd->is_arg_exist("y"))
  {
        double y = pcmd->get_number("y",0,0.0);
        iy = orig_mesh.find_skeleton_line_y(y);
  }
  if(pcmd->is_arg_exist("y.max") || pcmd->is_arg_exist("y.top"))
  {
        double ymax = pcmd->get_number("y.max",1,"y.top",orig_mesh.ymax);
        iymin = orig_mesh.find_skeleton_line_y(ymax);
  }
  if(pcmd->is_arg_exist("y.min") || pcmd->is_arg_exist("y.bottom"))
  {
        double ymin = pcmd->get_number("y.min",1,"y.bottom",orig_mesh.ymin);
        iymax = orig_mesh.find_skeleton_line_y(ymin);
  }
  if(pcmd->is_arg_exist("location"))
  {
      if(pcmd->is_arg_value("location","TOP"))
        iymin=iymax=0;
      else if(pcmd->is_arg_value("location","BOTTOM"))
        iymin=iymax=orig_mesh.IY-1;
      else if(pcmd->is_arg_value("location","LEFT"))
        ixmin=ixmax=0;
      else if(pcmd->is_arg_value("location","RIGHT"))
        ixmin=ixmax=orig_mesh.IX-1;
  }

  if(pcmd->is_arg_exist("direction"))
  {
      if(pcmd->is_arg_value("direction","Horizontal"))
        iymin=iymax=iy;
      else if(pcmd->is_arg_value("direction","Vertical"))
        ixmin=ixmax=ix;
  }
  if(ixmin!=ixmax && iymin!=iymax)
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " SEGMENT: segment must be horizontal or vertical!\n";
        gss_log.record();
        return 1;
  }
  if(ixmin==ixmax && iymin==iymax)
  {
        gss_log.string_buf()<<"line " <<pcmd->get_current_lineno()<< " SEGMENT: it seems you give a point instead of segment.\n";
        gss_log.record();
        return 1;
  }
  orig_mesh.set_segment_rectangle(ixmin,ixmax,iymin,iymax,label);

  return 0;

}



/* ----------------------------------------------------------------------------
 * BSolver::init_grid:  driven function of MESH
 */
int BSolver::init_grid(list<Cmd> &cmdlist)
{
  int  flag = 0;
  skeleton_line skeleton_line_x;
  skeleton_line skeleton_line_y;
  char *tri_cmd;
  gss_log.string_buf()<<"\nConstruct Initial Grid...\n";
  gss_log.record();
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
    if(pcmdbuf->is_current_cmd("MESH"))   // It's a MESH card
    {
        flag = 1;
        list<Cmd>::iterator    pcmd = pcmdbuf->get_current_cmd();
        strcpy(ModelFile,pcmd->get_string("modelfile",0,"GSS_Model.cgns"));
        tri_cmd = pcmd->get_string("triangle",0,"pzADq30Q");
    }
  if(!flag) return 0;

  // build x-y grid skeleton
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
  {
    if(pcmdbuf->is_current_cmd("XMESH"))   // It's a X.MESH card
      if(set_xmesh(pcmdbuf->get_current_cmd(),skeleton_line_x)) return 1;
    if(pcmdbuf->is_current_cmd("YMESH"))   // It's a Y.MESH card
      if(set_ymesh(pcmdbuf->get_current_cmd(),skeleton_line_y)) return 1;
  }

  //after the initializtion of x and y mesh lines, we can build rectangle mesh now.
  mesh_constructor orig_mesh(skeleton_line_x,skeleton_line_y);

  //use ELIMINATE and SPREAD card to do necessary modification to the rectangle mesh.
  //after that, use REGION card to set material region.
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
  {
    if(pcmdbuf->is_current_cmd("ELIMINATE"))   // do eliminate
      if(set_eliminate(pcmdbuf->get_current_cmd(),orig_mesh)) return 1;
    if(pcmdbuf->is_current_cmd("SPREAD"))      // do spread
      if(set_spread(pcmdbuf->get_current_cmd(),orig_mesh)) return 1;
    if(pcmdbuf->is_current_cmd("REGION"))      // set region here.
      if(set_region(pcmdbuf->get_current_cmd(),orig_mesh)) return 1;
  }

  orig_mesh.set_region_rectangle();
  orig_mesh.set_region_ellipse();
  orig_mesh.make_node_index();

  // set labeled segment(most of them are electrode) here.
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
  {
    if(pcmdbuf->is_current_cmd("SEGMENT"))
      if(set_segment(pcmdbuf->get_current_cmd(),orig_mesh)) return 1;
  }


  //call Triangle to build triangle mesh
  orig_mesh.do_mesh(tri_cmd);

  //write mesh to cgns file. It will be readed later.
  if(orig_mesh.to_cgns(ModelFile)) 
  {
    gss_log.string_buf()<<"Save mesh error. I Can't open cgns file "<<ModelFile<<"."<<endl;
    gss_log.record();
    return 1;
  }
  // ok read it into main mesh structure.
  // I'm too lazy to fill main mesh structure directly.
  if(import_cgns(ModelFile)) return 1;
  build_zone();
  if(setup_bc())   return 1;
  if(build_zonedata())    return 1;
  setup_doping();
  import_mole_from_cgns(ModelFile);
  if(setup_init_data()) return 1;
  reorder();
  zone_to_field();
  build_least_squares();

  return 0;
}


