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
/*  Last update: May 14, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "typedef.h"
#include "skeleton.h"
#include "material.h"
#include "mathfunc.h"
#include <stdio.h>
#include <cgnslib.h>

mesh_constructor::mesh_constructor(skeleton_line &lx,skeleton_line &ly)
{
  // init Triangle io structure.
  in.pointlist = (double *) NULL;
  in.pointattributelist = (double *) NULL;
  in.pointmarkerlist = (int *) NULL;
  in.segmentlist = (int *) NULL;
  in.segmentmarkerlist = (int *) NULL;
  in.regionlist = (double *)NULL;
  out.pointlist = (double *) NULL;
  out.pointattributelist = (double *) NULL;
  out.pointmarkerlist = (int *) NULL;
  out.trianglelist = (int *) NULL;
  out.triangleattributelist = (double *) NULL;
  out.segmentlist = (int *) NULL;
  out.segmentmarkerlist = (int *) NULL;

  IX = lx.point_list.size();
  IY = ly.point_list.size();
  point_num = IX*IY;
  aux_point_num = 0;
  // set up the 2d array for gird points
  point_array2d = new skeleton_point *[IY];
  for(int i=0;i<IY;i++)
    point_array2d[i]= new skeleton_point [IX];
  for(int i=0;i<IX;i++)
    for(int j=0;j<IY;j++)
    {
      point_array2d[j][i].set_location(lx.point_list[i].x,ly.point_list[j].y);
    }

  xmin  =  point_array2d[0][0].x;
  xmax  =  point_array2d[0][IX-1].x;
  ymax  =  point_array2d[0][0].y;
  ymin  =  point_array2d[IY-1][0].y;

  width = xmax-xmin;
  depth = ymax-ymin;
}

mesh_constructor::~mesh_constructor()
{
  for(int i=0;i<IY;i++)
    delete [] point_array2d[i];
  delete [] point_array2d;
  free(in.pointlist);
  free(in.pointmarkerlist);
  free(in.pointattributelist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);
  free(in.regionlist);

  free(out.pointlist);
  free(out.pointmarkerlist);
  free(out.pointattributelist);
  free(out.trianglelist);
  free(out.triangleattributelist);
  free(out.segmentlist);
  free(out.segmentmarkerlist);
}

int mesh_constructor::find_skeleton_line_x(double x)
{
  int ix=0;
  double dx=1e10;
  for(int i=0;i<IX;i++)
    if(fabs(x-point_array2d[0][i].x)<dx)
    {
      dx=fabs(x-point_array2d[0][i].x);
      ix=i;
    }
  return ix;
}

int mesh_constructor::find_skeleton_line_y(double y)
{
  int iy=0;
  double dy=1e10;
  for(int i=0;i<IY;i++)
    if(fabs(y-point_array2d[i][0].y)<dy)
    {
      dy=fabs(y-point_array2d[i][0].y);
      iy=i;
    }
  return iy;
}


int mesh_constructor::find_move_skeleton_line_x(double x)
{
  int ix=find_skeleton_line_x(x);
  for(int j=0;j<IY;j++)
    point_array2d[j][ix].x=x;
  return ix;
}


int mesh_constructor::find_move_skeleton_line_y(double y)
{
  int iy=find_skeleton_line_y(y);
  for(int i=0;i<IX;i++)
    point_array2d[iy][i].y=y;
  return iy;
}


int mesh_constructor::x_eliminate(int ixmin,int ixmax, int iymin, int iymax)
{
  for(int i=ixmin;i<=ixmax;i++)
  {
    int eliminate_flag = 1;
    for(int j=max(1,iymin+1);j<min(iymax+1,IY-1);j++)
      if(!point_array2d[j][i].eliminated)
      {
        if(eliminate_flag==1)
        {
          point_array2d[j][i].eliminated=1;
          eliminate_flag=0;
        }
        else
          eliminate_flag=1;
      }
  }
  return 0;
}

int mesh_constructor::y_eliminate(int iymin,int iymax, int ixmin, int ixmax)
{
  for(int j=iymin;j<=iymax;j++)
  {
    int eliminate_flag = 1;
    for(int i=max(1,ixmin+1);i<min(ixmax+1,IX-1);i++)
      if(!point_array2d[j][i].eliminated)
      {
        if(eliminate_flag==1)
        {
          point_array2d[j][i].eliminated=1;
          eliminate_flag=0;
        }
        else
          eliminate_flag=1;
      }
  }
  return 0;
}


int mesh_constructor::set_spread_rectangle(int location,double width,int upperline,int lowerline,
                          double yuploc,double yloloc,double encroach,double grading)
{
  double ybtanc = point_array2d[IY-1][0].y;
  double yupanc = point_array2d[upperline][0].y;
  double yloanc = point_array2d[lowerline][0].y;

  double xloc,yupold,yloold,ybot;
  if(location == LEFT)
    xloc = point_array2d[0][0].x + width;
  else
    xloc = point_array2d[0][IX-1].x - width;

  double  erfarg,erfar2,erfval,erfvl2;

  //proc colume by colume
  for(int i=0;i<IX;i++)
  {
      double xco=point_array2d[0][i].x;

      //evaluate error function for this coord.
      erfarg=(xco-xloc)*encroach;
      if (location==RIGHT)
         erfar2=1.5*(erfarg+0.6);
      else
         erfar2=1.5*(erfarg-0.6);
      if (location==LEFT)
      {
        erfval=erfc(erfarg);
        erfvl2=erfc(erfar2);
      }
      else
      {
        erfval=erfc(-erfarg);
        erfvl2=erfc(-erfar2);
      }
      erfval=erfval*0.5;
      erfvl2=erfvl2*0.5;

      //
      // compute new node locations on this column

      //  get upper, lower, bottom current loc.

      yupanc = yupold;
      yupold = point_array2d[upperline][i].y;
      yloanc = yloold;
      yloold = point_array2d[lowerline][i].y;
      ybtanc = ybot;
      ybot = point_array2d[IY-1][i].y;

      // compute upward shift and downward
      double deltup=erfval*(yupold-yuploc);
      double deltlo=erfval*(yloloc-yloold);

      //  compute new locations
      double yupnew=yupold-deltup;
      double ylonew=yloold+deltlo;

      // compute old and new spreads of middle, bottom
      double spmdol=yloold-yupold;
      double spmdnw=spmdol+deltup+deltlo;
      double spbtol=ybot-yloold;
      double spbtnw=spbtol-deltlo;

      // grading ratio
      double y,dy;
      if(grading!=1.0)
      {
        y  =  yupnew;
        dy = (ylonew-yupnew)*(grading-1)/(pow(grading,lowerline-upperline)-1);
      }

      // scan y nodes and move as necessary
      for(int j=0;j<IY;j++)
      {
        // get node number and y-coord.
        double yco=point_array2d[j][i].y;

        // consider top, middle, and bottom regions

        //  top region (j<upperline) shift all by deltup
        if(j<=upperline)
        {
          yco=yco-deltup;
        }
        // bottom region, spread proportionally
        else if (j>=lowerline)
        {
          double rat=(yco-yloold)/spbtol;
          yco=ylonew+rat*spbtnw;
        }
        // middle region spread proportionally unless new grading
        // is requested.
        else
        {
          if(grading==1.0)
          {
            double rat=(yco-yupold)/spmdol;
            yco=yupnew+rat*spmdnw;
          }
          else // new grading requested
          {
                double ycordg = y+dy;
                y += dy;
                dy*= grading;
                // vary from new grading to proportional grading
                // based on erfvl2 (2/3 the spread of erfval and
                // centered at 30% point of erfval instead of
                // 50% point.
                yco=yco+erfvl2*(ycordg-yco);
          }
        }
        point_array2d[j][i].y = yco;
      }
  }
  return 0;
}



void mesh_constructor::make_node_index()
{
  point_num=0;
  for(int i=0;i<IX;i++)
    for(int j=0;j<IY;j++)
    {
      if(!point_array2d[j][i].eliminated)
        point_array2d[j][i].index=point_num++;
      else
        point_array2d[j][i].index=-1;
    }
  for(int i=0;i<aux_point_array1d.size();i++)
     aux_point_array1d[i].index = point_num++;
}


int  mesh_constructor::set_region_rectangle()
{
  // set inner point for each region;
  for(int i=0;i<IX-1;i++)
    for(int j=0;j<IY-1;j++)
    {
      for(int k=region_array1d.size()-1;k>=0;k--)
	if(region_array1d[k].shape==Rectangle &&
	    i>=region_array1d[k].ixmin && i+1<=region_array1d[k].ixmax &&
            j>=region_array1d[k].iymin && j+1<=region_array1d[k].iymax)
        {
          region_array1d[k].px.push_back((point_array2d[j][i].x+point_array2d[j][i+1].x)/2);
          region_array1d[k].py.push_back((point_array2d[j][i].y+point_array2d[j+1][i].y)/2);
          break;
        }
    }
  // set segment bounding
  map<skeleton_edge, int, lt_edge>::iterator pt;
  for(int r=0;r<region_array1d.size();r++)
  {
    if(region_array1d[r].shape!=Rectangle) continue;
    int ixmin = region_array1d[r].ixmin;
    int ixmax = region_array1d[r].ixmax;
    int iymin = region_array1d[r].iymin;
    int iymax = region_array1d[r].iymax;
    skeleton_edge edge;
    edge.IX = IX;
    edge.IY = IY;
    skeleton_segment segment;
    sprintf(segment.segment_label,"%s_Neumann",region_array1d[r].label);
    segment.segment_mark=r+1;
    segment_array1d.push_back(segment);
    //process bottom line
    for(int i=ixmin;i<ixmax;)
    {
      edge.p1[0]=i++;
      edge.p1[1]=iymax;
      while(point_array2d[iymax][i].eliminated) i++;
      edge.p2[0]=i;
      edge.p2[1]=iymax;
      edge_table[edge]=r+1;
      //edge_table.insert(pair<skeleton_edge const,int>(edge,r+1));
    }
    //process top line
    for(int i=ixmin;i<ixmax;)
    {
      edge.p1[0]=i++;
      edge.p1[1]=iymin;
      while(point_array2d[iymin][i].eliminated) i++;
      edge.p2[0]=i;
      edge.p2[1]=iymin;
      edge_table[edge]=r+1;
      //edge_table.insert(pair<skeleton_edge const,int>(edge,r+1));
    }
    //process left line
    for(int i=iymin;i<iymax;)
    {
      edge.p1[0]=ixmin;
      edge.p1[1]=i++;
      while(point_array2d[i][ixmin].eliminated) i++;
      edge.p2[0]=ixmin;
      edge.p2[1]=i;
      edge_table[edge]=r+1;
      //edge_table.insert(pair<skeleton_edge const,int>(edge,r+1));
    }
    //process right line
    for(int i=iymin;i<iymax;)
    {
      edge.p1[0]=ixmax;
      edge.p1[1]=i++;
      while(point_array2d[i][ixmax].eliminated) i++;
      edge.p2[0]=ixmax;
      edge.p2[1]=i;
      edge_table[edge]=r+1;
      //edge_table.insert(pair<skeleton_edge const,int>(edge,r+1));
    }
  }
  return 0;
}


int  mesh_constructor::set_segment_rectangle(int location,int ixmin,int ixmax,int iymin,int iymax,char *label)
{
  int segment_mark = region_array1d.size()+segment_array1d.size()+1;
  skeleton_segment segment;
  strcpy(segment.segment_label,label);
  segment.segment_mark=segment_mark;
  segment_array1d.push_back(segment);
  skeleton_edge edge;
  edge.IX = IX;
  edge.IY = IY;
  if(location==TOP)
  {
    //process top line
    for(int i=ixmin;i<ixmax;)
    {
      while(point_array2d[0][i].eliminated) i++;
      edge.p1[0]=i++;
      edge.p1[1]=0;
      while(point_array2d[0][i].eliminated) i++;
      edge.p2[0]=i;
      edge.p2[1]=0;
      edge_table[edge]=segment_mark;
    }
  }

  if(location==BOTTOM)
  {
    //process bottom line
    for(int i=ixmin;i<ixmax;)
    {
      while(point_array2d[IY-1][i].eliminated) i++;
      edge.p1[0]=i++;
      edge.p1[1]=IY-1;
      while(point_array2d[IY-1][i].eliminated) i++;
      edge.p2[0]=i;
      edge.p2[1]=IY-1;
      edge_table[edge]=segment_mark;
    }
  }
  if(location==LEFT)
  {
    //process left line
    for(int i=iymin;i<iymax;)
    {
      while(point_array2d[i][0].eliminated) i++;
      edge.p1[0]=0;
      edge.p1[1]=i++;
      while(point_array2d[i][0].eliminated) i++;
      edge.p2[0]=0;
      edge.p2[1]=i;
      edge_table[edge]=segment_mark;
    }
  }
  if(location==RIGHT)
  {
    //process right line
    for(int i=iymin;i<iymax;)
    {
      while(point_array2d[i][IX-1].eliminated) i++;
      edge.p1[0]=IX-1;
      edge.p1[1]=i++;
      while(point_array2d[i][IX-1].eliminated) i++;
      edge.p2[0]=IX-1;
      edge.p2[1]=i;
      edge_table[edge]=segment_mark;
    }
  }

  return 0;
}


int  mesh_constructor::set_segment_rectangle(int ixmin,int ixmax,int iymin,int iymax,char *label)
{
  int segment_mark = region_array1d.size()+segment_array1d.size()+1;
  skeleton_segment segment;
  strcpy(segment.segment_label,label);
  segment.segment_mark=segment_mark;
  segment_array1d.push_back(segment);
  skeleton_edge edge;
  edge.IX = IX;
  edge.IY = IY;
  if(iymin==iymax) //horizontal
  {
    //process horizontal line
    for(int i=ixmin;i<ixmax;)
    {
      while(point_array2d[iymin][i].eliminated) i++;
      edge.p1[0]=i++;
      edge.p1[1]=iymin;
      while(point_array2d[iymin][i].eliminated) i++;
      edge.p2[0]=i;
      edge.p2[1]=iymin;
      edge_table[edge]=segment_mark;
    }
  }
  else if(ixmin==ixmax) //vertical
  {
    //process vertical line
    for(int i=iymin;i<iymax;)
    {
      while(point_array2d[i][ixmin].eliminated) i++;
      edge.p1[0]=ixmin;
      edge.p1[1]=i++;
      while(point_array2d[i][ixmin].eliminated) i++;
      edge.p2[0]=ixmin;
      edge.p2[1]=i;
      edge_table[edge]=segment_mark;
    }
  }

  return 0;
}

int mesh_constructor::set_region_ellipse()
{
  
  for(int r=0; r<region_array1d.size();r++)
   if(region_array1d[r].shape==Ellipse)
   {
     //set aux point 
     double theta  = region_array1d[r].theta;
     double theta0 = region_array1d[r].theta;
     double delta_angle = 2*3.14159265359/region_array1d[r].division;
     aux_point tmp_point;
     aux_edge  tmp_edge;
     int   offset=aux_point_array1d.size();
     for(int i=0;i<region_array1d[r].division;i++)
     {
        tmp_point.x = region_array1d[r].centrex + region_array1d[r].major_radii*cos(theta);
	tmp_point.y = region_array1d[r].centrey + region_array1d[r].minor_radii*sin(theta);
	theta+=delta_angle;
	tmp_edge.p1=aux_point_array1d.size();
	aux_point_array1d.push_back(tmp_point);
	tmp_edge.p2=aux_point_array1d.size();
	if(tmp_edge.p2==offset+region_array1d[r].division) tmp_edge.p2=offset;
	// set segment bounding
        tmp_edge.segment_mark=r+1;
	aux_edge_array1d.push_back(tmp_edge);
     }
     //do necessary elimination of neibour points to improve mesh quality.
     for(int i=1;i<IX-1;i++)
      for(int j=1;j<IY-1;j++)
     {
       double x=point_array2d[j][i].x;
       double y=point_array2d[j][i].y;
       double xc=region_array1d[r].centrex;
       double yc=region_array1d[r].centrey;
       double char_length =  (fabs(point_array2d[j][i+1].x-point_array2d[j][i-1].x)+
                              fabs(point_array2d[j+1][i].y-point_array2d[j-1][i].y))/4.0;
       double rsmin=pow(region_array1d[r].minor_radii-char_length,2);
       double rsmax=pow(region_array1d[r].minor_radii+char_length,2);
       double rlmin=pow(region_array1d[r].major_radii-char_length,2);
       double rlmax=pow(region_array1d[r].major_radii+char_length,2);
       if( pow(x-xc,2)/rlmax + pow(y-yc,2)/rsmax < 1 &&
           pow(x-xc,2)/rlmin + pow(y-yc,2)/rsmin > 1  ) 
	 point_array2d[j][i].eliminated=1;	      
     }
     
     region_array1d[r].px.push_back(region_array1d[r].centrex);
     region_array1d[r].py.push_back(region_array1d[r].centrey);
   }  

   return 0;

}


int mesh_constructor::do_mesh(char *tri_arg)
{
  //set point
  in.numberofpoints = point_num;
  in.numberofpointattributes = 0;
  in.pointattributelist = (double *)NULL;
  in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
  in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
  double *ppointlist=in.pointlist;
  int *ppointmarkerlist=in.pointmarkerlist;
  //the points belongs to rectangle region 
  for(int i=0;i<IX;i++)
    for(int j=0;j<IY;j++)
      if(!point_array2d[j][i].eliminated)
      {
        *ppointlist++ = point_array2d[j][i].x;
        *ppointlist++ = point_array2d[j][i].y;
        *ppointmarkerlist++ = 0;
      }
  //the points belongs to the boundary of ellipse region     
  for(int i=0;i<aux_point_array1d.size();i++)
  {
        *ppointlist++ = aux_point_array1d[i].x;
        *ppointlist++ = aux_point_array1d[i].y;
        *ppointmarkerlist++ = 0;
  }
  
  //do necessarily prepare for call triangulate  
  in.numberoftriangles = 0;
  in.numberofcorners = 3;
  in.numberoftriangleattributes = 0;
  in.trianglelist =  (int *) NULL;
  in.trianglearealist = (double *) NULL;
  in.triangleattributelist = NULL;

  in.numberofsegments = edge_table.size()+aux_edge_array1d.size();
  in.segmentlist =  (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
  int *psegmentlist =  in.segmentlist;
  int *psegmentmarkerlist = in.segmentmarkerlist;
  map<skeleton_edge, int, lt_edge>::iterator pt = edge_table.begin();
  for(int i=0;i<edge_table.size();i++)
  {
    *psegmentlist++ = point_array2d[pt->first.p1[1]][pt->first.p1[0]].index;
    *psegmentlist++ = point_array2d[pt->first.p2[1]][pt->first.p2[0]].index;
    *psegmentmarkerlist++ = pt->second;
    pt++;
  }
  for(int i=0;i<aux_edge_array1d.size();i++)
  {
    *psegmentlist++ = aux_point_array1d[aux_edge_array1d[i].p1].index;
    *psegmentlist++ = aux_point_array1d[aux_edge_array1d[i].p2].index;
    *psegmentmarkerlist++ = aux_edge_array1d[i].segment_mark;
    pt++;
  }
  
  in.numberofholes = 0;
  in.numberofregions = region_array1d.size();
  in.regionlist = (double *) malloc(in.numberofregions * 4 * sizeof(double));
  double *pregionlist =  in.regionlist;
  for(int i=0;i<in.numberofregions;i++)
  {
    *pregionlist++ = region_array1d[i].px[0];
    *pregionlist++ = region_array1d[i].py[0];
    *pregionlist++ = double(i);
    *pregionlist++ = 0;
  }

  // Refine the triangulation according to the attached
  // triangle area constraints.

  triangulate(tri_arg, &in, &out, (struct triangulateio *) NULL);

  for(int i=0;i<out.numberofsegments;i++)
  {
    output_edge edge;
    edge.index = i;
    edge.p1 = out.segmentlist[2*i+0];
    edge.p2 = out.segmentlist[2*i+1];
    edge.interface = -1;
    edge.r1 = -1;
    edge.r2 = -1;
    edge.bc_type = -1;
    edge.mark = out.segmentmarkerlist[i];
    out_edge.push_back(edge);
  }

  return 0;
}

typedef struct
{
      int r1;
      int r2;
      vector<int> node_r1;
      vector<int> node_r2;
      int bc_mark;
}Zone_conn;


int mesh_constructor::to_cgns(const char *filename)
{
  int fn,B,Z,C,S,BC,I,SOL,F;
  int size[3];
  double *x,*y;
  int *elem;
  char bcname[32];
  char conn_name[32];
  //this is the main part of the work. the nodes will be classified by
  //its region and insert to each zone

  int  region_num = region_array1d.size();       //total region;
  int **local_index = new int*[region_num];      //reg_index = local_index[reg][global_index]
  int **local_mark = new int*[region_num];       //set local bc mark
  int **conn = new int*[region_num];             //for zone to zone connectivity information

  for(int r=0; r<region_num; r ++)
  {
    local_index[r] = new int[out.numberofpoints];
    local_mark[r] = new int[out.numberofpoints];
    conn[r] = new int[out.numberofpoints];
  }

  int *af = new int[out.numberofpoints];

  for(int r=0;r<region_num;r++)
  {
    int *pf = af;
    for(int j = 0; j < out.numberofpoints; j++)
    {
       local_index[r][j] = -1; //-1 means node not in this region
       local_mark[r][j] = -1; // -1 means no bc_mark
       *pf++ = 0;        //visited flag array
    }
    region_array1d[r].node_num = 0;
    region_array1d[r].tri_num = 0;

    for(int j = 0; j < out.numberoftriangles; j++)     //search for all the nodes
      if (int(out.triangleattributelist[j]+0.5) == r)
      {
        region_array1d[r].tri_num++;
        int nodeA = out.trianglelist[3*j+0];
        int nodeB = out.trianglelist[3*j+1];
        int nodeC = out.trianglelist[3*j+2];
        if (!af[nodeA]) { af[nodeA] = 1;local_index[r][nodeA] = region_array1d[r].node_num++; }
        if (!af[nodeB]) { af[nodeB] = 1;local_index[r][nodeB] = region_array1d[r].node_num++; }
        if (!af[nodeC]) { af[nodeC] = 1;local_index[r][nodeC] = region_array1d[r].node_num++; }
      }
  }
  delete [] af;

  for(int i = 0; i < out_edge.size(); i++)     //search for all the edges
    for(int r=0;r<region_num;r++)
    {
      int p1 = out_edge[i].p1;
      int p2 = out_edge[i].p2;
      if(local_index[r][p1]!=-1 && local_index[r][p2]!=-1)
      {
        if(out_edge[i].r1==-1) out_edge[i].r1 = r;
      else {out_edge[i].interface = 1; out_edge[i].r2 = r;}
        region_array1d[r].boundary.push_back(i);
      }
    }

  // remove old file if exist
  remove(filename);

  // open CGNS file for write
  if(cg_open(filename,MODE_WRITE,&fn))
  {
    return 1;
  }
  // create base (can give any name)
  cg_base_write(fn,"GSS_Mesh",2,2,&B); /*two dimension*/

  // create zone
  for(int r=0;r<region_num;r++)
  {

    size[0] = region_array1d[r].node_num;
    size[1] = region_array1d[r].tri_num;
    size[2] = 0;

    x = new double[size[0]];
    y = new double[size[0]];
    elem= new int[4*size[1]];
    //zone name set as region name
    cg_zone_write(fn,B,region_array1d[r].label,size,Unstructured,&Z);

    cg_goto(fn,B,"Zone_t",Z,"end");
    cg_descriptor_write("RegionType",region_array1d[r].material);
    // write grid coordinates
    for(int j=0;j<out.numberofpoints;j++)
      if(local_index[r][j]>-1)  //process node in this region
      {
        //convert the unit to cm
        x[local_index[r][j]] = out.pointlist[2*j+0]*1e-4;
        y[local_index[r][j]] = out.pointlist[2*j+1]*1e-4;
      }

    cg_coord_write(fn,B,Z,RealDouble,"CoordinateX",x,&C);
    cg_coord_write(fn,B,Z,RealDouble,"CoordinateY",y,&C);

    // set element connectivity here

    int *pelem=elem;
    for(int j = 0; j < out.numberoftriangles; j++)
      if (int(out.triangleattributelist[j]+0.5) == r)  //is the triangle belongs to this region?
      {
        *pelem++ = TRI_3;
        *pelem++ = local_index[r][out.trianglelist[3*j+0]]+1;
        *pelem++ = local_index[r][out.trianglelist[3*j+1]]+1;
        *pelem++ = local_index[r][out.trianglelist[3*j+2]]+1;
      }

    cg_section_write(fn,B,Z,"GridElements",MIXED,1,size[1],0,elem,&S);

    delete [] x;
    delete [] y;
    delete [] elem;


    // write  boundary segment
    int *asegment = new int[2*out.numberofpoints];
    int *psegment = asegment;
    int  segrange[2];
    segrange[0] = size[1]+1;
    int  csegment;
    vector<Zone_conn> zone_conn_array;
    // write labeled boundary
    for(int j=0;j<segment_array1d.size();j++)
    {
      if(segment_array1d[j].segment_mark > region_array1d.size()) //for labeled segment
      {
        psegment = asegment;
        csegment = 0;
        for(int k=0;k<region_array1d[r].boundary.size();k++)
        {
          int segment = region_array1d[r].boundary[k];
          if(out_edge[segment].mark == segment_array1d[j].segment_mark) //the edge belongs to this labeled segment
          {
            out_edge[segment].bc_type = LabeledSegment;
            int nodeA = out_edge[segment].p1;
            int nodeB = out_edge[segment].p2;
            //only record the nodes that belong to this region
            if(local_index[r][nodeA]==-1||local_index[r][nodeB]==-1 )continue;

            *psegment++ = local_index[r][nodeA]+1;
            *psegment++ = local_index[r][nodeB]+1;
            csegment++;
            local_mark[r][nodeA] = j;
            local_mark[r][nodeB] = j;
          }
        }

        if(csegment>0)
        {
          segrange[1] = segrange[0] + csegment - 1;
          cg_section_write(fn,B,Z,segment_array1d[j].segment_label,BAR_2,
                           segrange[0],segrange[1],0,asegment,&S);
          cg_boco_write(fn,B,Z,segment_array1d[j].segment_label,BCTypeUserDefined,
                        ElementRange,2,segrange,&BC);
          segrange[0] = segrange[1]+1;
        }
      }
    }
    //record zone to zone connectivity data if necessary
    if(region_num>1)
    {   //build connectivity data here
      for(int k=0; k<region_num; k++)
      {
          if(k==r) continue;
          Zone_conn zone_conn;
          zone_conn.r1 = r;
          zone_conn.r2 = k;
          for (int s=0; s<segment_array1d.size(); s++)
          {
            for(int j=0; j<out.numberofpoints; j++)
              if(local_index[k][j]!=-1 && local_index[r][j]!=-1 && local_mark[r][j]==s)
            {
              zone_conn.node_r1.push_back(local_index[r][j]+1);
              zone_conn.node_r2.push_back(local_index[k][j]+1);
            }
            zone_conn.bc_mark = s;
            if(zone_conn.node_r1.size()>=2)
              zone_conn_array.push_back(zone_conn);
            zone_conn.node_r1.clear();
            zone_conn.node_r2.clear();
          }
      }
    }

    // write Interface as boundary
    for(int j=0;j<region_num;j++)
    {
      if(j!=r)
      {
        csegment = 0;
        psegment = asegment;
        for(int k=0; k<region_array1d[j].boundary.size(); k++)
        {
          //the edge belongs to region j
          int segment = region_array1d[j].boundary[k];
          int nodeA = out_edge[region_array1d[j].boundary[k]].p1;
          int nodeB = out_edge[region_array1d[j].boundary[k]].p2;
          if(local_index[r][nodeA]==-1||local_index[r][nodeB]==-1 )continue; // does the edge belongs to region r ?
          if(out_edge[segment].bc_type == LabeledSegment) continue;          // skip labeled segment
          out_edge[segment].bc_type = Interface;       // ok , it's an interface edge.
          *psegment++ = local_index[r][nodeA]+1;
          *psegment++ = local_index[r][nodeB]+1;
          local_mark[r][nodeA] = Interface;
          local_mark[r][nodeB] = Interface;
          csegment++;
        }

        if(csegment>0)
        {
          //the name is alpha ordered
          if(strcmp(region_array1d[r].label,region_array1d[j].label)<0)
            snprintf(bcname,31,"IF_%s_to_%s",region_array1d[r].label,region_array1d[j].label);
          else
            snprintf(bcname,31,"IF_%s_to_%s",region_array1d[j].label,region_array1d[r].label);

          segrange[1] = segrange[0] + csegment - 1;
          cg_section_write(fn,B,Z,bcname,BAR_2,
                           segrange[0],segrange[1],0,asegment,&S);
          cg_boco_write(fn,B,Z,bcname,BCTypeUserDefined,ElementRange,2,segrange,&BC);
          segrange[0] = segrange[1]+1;
        }
      }
    }

    //record zone to zone connectivity data if necessary
    if(region_num>1)
    {   //build connectivity data here
      for(int k=0; k<region_num; k++)
      {
          if(k==r) continue;
          Zone_conn zone_conn;
          zone_conn.r1 = r;
          zone_conn.r2 = k;
          for(int j=0; j<out.numberofpoints; j++)
            if(local_index[k][j]!=-1 && local_index[r][j]!=-1 && local_mark[r][j]==Interface)
            {
              zone_conn.node_r1.push_back(local_index[r][j]+1);
              zone_conn.node_r2.push_back(local_index[k][j]+1);
            }
          zone_conn.bc_mark = Interface;
          if(zone_conn.node_r1.size()>=2)
            zone_conn_array.push_back(zone_conn);
          zone_conn.node_r1.clear();
          zone_conn.node_r2.clear();
      }
    }

    // write other boundary segments as Neumann boundary
    psegment = asegment;
    csegment = 0;

    for(int j=0;j<region_array1d[r].boundary.size();j++)
    {
      int segment = region_array1d[r].boundary[j];
      int nodeA = out_edge[segment].p1;
      int nodeB = out_edge[segment].p2;
      //only record the nodes that belong to this region
      if(local_index[r][nodeA]==-1||local_index[r][nodeB]==-1 )continue;
      //skip labeled segment and interface
      if(out_edge[segment].bc_type == LabeledSegment) continue;
      if(out_edge[segment].bc_type == Interface) continue;
      *psegment++ = local_index[r][nodeA]+1;
      *psegment++ = local_index[r][nodeB]+1;
      csegment++;
    }
    if(csegment>0)
    {
      snprintf(bcname,31,"%s_Neumann",region_array1d[r].label);
      segrange[1] = segrange[0] + csegment - 1;
      cg_section_write(fn,B,Z,bcname,BAR_2,
                       segrange[0],segrange[1],0,asegment,&S);
      cg_boco_write(fn,B,Z,bcname,BCTypeUserDefined,ElementRange,2,segrange,&BC);
      segrange[0] = segrange[1]+1;
    }
    delete [] asegment;

    
    for(int j=0; j<zone_conn_array.size(); j++)
    {

        int conn_node = zone_conn_array[j].node_r1.size();
        int r2 = zone_conn_array[j].r2;
        for(int k=0;k<conn_node;k++)
        {
           conn[r][k] = zone_conn_array[j].node_r1[k];
           conn[r2][k] = zone_conn_array[j].node_r2[k];
        }
        if(zone_conn_array[j].bc_mark==Interface)
        {
          //the conn_name is alpha ordered
          if(strcmp(region_array1d[r].label,region_array1d[r2].label)<0)
            snprintf(conn_name,31,"IF_%s_to_%s",region_array1d[r].label,region_array1d[r2].label);
          else
            snprintf(conn_name,31,"IF_%s_to_%s",region_array1d[r2].label,region_array1d[r].label);
          cg_conn_write(fn,B,Z,conn_name,Vertex, Abutting1to1,
                      PointList, conn_node, conn[r],region_array1d[r2].label,Unstructured,
                      PointListDonor,Integer,conn_node, conn[r2],&I);
        }
        else
        {
          cg_conn_write(fn,B,Z,segment_array1d[zone_conn_array[j].bc_mark].segment_label,Vertex, Abutting1to1,
                      PointList, conn_node, conn[r],region_array1d[r2].label,Unstructured,
                      PointListDonor,Integer,conn_node, conn[r2],&I);
        }
    }

    if(IsSingleCompSemiconductor(region_array1d[r].material))
    {
        double *mole_x;
        mole_x = new double[region_array1d[r].node_num];
        if(region_array1d[r].mole_x1_grad == MOLE_GRAD_X)
        {
          for(int j=0;j<out.numberofpoints;j++)
            if(local_index[r][j]>-1)  //process node in this region
          {
              double x = out.pointlist[2*j+0];
              double mole = region_array1d[r].mole_x1
                          + region_array1d[r].mole_x1_slope*(x-region_array1d[r].xmin);
              mole_x[local_index[r][j]] = dmax(mole,0);
          }

        }
        else if(region_array1d[r].mole_x1_grad == MOLE_GRAD_Y)
        {
          for(int j=0;j<out.numberofpoints;j++)
            if(local_index[r][j]>-1)  //process node in this region
          {
              double y = out.pointlist[2*j+1];
              double mole = region_array1d[r].mole_x1
                          - region_array1d[r].mole_x1_slope*(y-region_array1d[r].ymax);
              mole_x[local_index[r][j]] = dmax(mole,0);
          }
        }

        cg_sol_write(fn,B,Z,"Mole",Vertex,&SOL);
        cg_field_write(fn,B,Z,SOL,RealDouble,"mole_x",mole_x,&F);
        delete [] mole_x;
    }

  }

  //all things are done
  for(int r =0; r<region_num; r ++)
  {
    delete [] local_index[r];
    delete [] local_mark[r];
    delete [] conn[r];
  }
  delete  [] local_index;
  delete  [] local_mark;
  delete  [] conn;
  // close CGNS file
  cg_close(fn);
  return 0;
}
