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


#include "element.h"
#include <algorithm>
using namespace std;

void   Segment::set_node_array()
{
  node_array.clear();
  if(!edge_array.size()) return;
  vector<int> edge_array_flag(edge_array.size(),0);
  int p1,p2;
  for(int i=0; i<edge_array.size();i++)
  {
    if(!edge_array_flag[0])
    { 
      //insert the first edge
      p1=edge_array[0].p1;
      p2=edge_array[0].p2;
      node_array.insert(node_array.begin(),p1);
      node_array.push_back(p2);
      edge_array_flag[0]=1;
    }
    //search the edge which has the point p1, inset another point of the edge in front of p1
    for(int j=0; j<edge_array.size();j++)
    {
        if(edge_array_flag[j]) continue;
	if(edge_array[j].p1==p1)  
	{  
	   node_array.insert(node_array.begin(),edge_array[j].p2);
	   edge_array_flag[j]=1;
	   p1=edge_array[j].p2;
	   break;
	} 
	if(edge_array[j].p2==p1)  
	{ 
	  node_array.insert(node_array.begin(),edge_array[j].p1);
	  edge_array_flag[j]=1;
	  p1=edge_array[j].p1;
	  break;
	}
    }
    //search the edge which has the point p2, inset another point of the edge at the end of p2
    for(int j=0; j<edge_array.size();j++)
    {
        if(edge_array_flag[j]) continue;	
	if(edge_array[j].p1==p2)  
	{
	  node_array.push_back(edge_array[j].p2);
	  edge_array_flag[j]=1;
	  p2=edge_array[j].p2;
	  break;
	} 
	if(edge_array[j].p2==p2)  
	{
	  node_array.push_back(edge_array[j].p1);
	  edge_array_flag[j]=1;
	  p2=edge_array[j].p1;
	  break;
	}  
    }
  }
}

int   Segment::get_edge_bc_index(int p1, int p2)
{
  for(int i=0; i<edge_array.size(); i++)
  {
        if(edge_array[i].p1 == p1 && edge_array[i].p2 == p2)
                return bc_index;
        if(edge_array[i].p2 == p1 && edge_array[i].p1 == p2)
                return bc_index;
  }
  return 0;
}

VoronoiCell::VoronoiCell():  nb_num(0) //constructor
{ nb_array  = (int*) NULL;
  inb_array = (int*) NULL;
  elen      = (double*) NULL;
  ilen      = (double*) NULL;
  angle     = (double*) NULL;
  celledge  = (int*) NULL;
}

VoronoiCell::VoronoiCell(const VoronoiCell& org) //copy constructor
{
  x =  org.x;     y =  org.y;
  nb_num =  org.nb_num;
  area   =  org.area;
  bc_index = org.bc_index;
  sa = org.sa;
  sb = org.sb;
  sc = org.sc;
  nb_array  = (int*) NULL;
  inb_array = (int*) NULL;
  elen      = (double*) NULL;
  ilen      = (double*) NULL;
  angle     = (double*) NULL;
  celledge  = (int*) NULL;
  if(nb_num)
  {
    nb_array = new int[nb_num];
    inb_array= new int[nb_num];
    elen = new double[nb_num];
    ilen = new double[nb_num];
    angle= new double[nb_num];
    celledge = new int[nb_num];
    for(int i=0;i<nb_num; i++)
    {
      nb_array[i] = org.nb_array[i];
      inb_array[i] = org.inb_array[i];
      elen[i] = org.elen[i];
      ilen[i] = org.ilen[i];
      angle[i] = org.angle[i];
      celledge[i] = org.celledge[i];
    }
  }
}

VoronoiCell::~VoronoiCell()
{
  delete [] nb_array;
  delete [] inb_array;
  delete [] elen;
  delete [] ilen;
  delete [] angle;
  delete [] celledge;
}
