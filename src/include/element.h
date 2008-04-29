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

#ifndef _element_h_
#define _element_h_

#include "typedef.h"
#include <vector>
using namespace std;

class   Node                            // node
{
public:
  double          x,y;                  // x and y coordinate
  int             bc_index;             // the index to boundary of the node
  int             zone_index;           // the zone index;
  int             g_index;              // the index of this node in gnode structure
};

class   GNode   :public Node            // global node
{
public:
  int             local_index;          // the index in the zone
};


class Edge                              // boundary edge
{
public:
  int p1,p2;
  int bc_index;
  int tri_index;                        // this segment edge belongs to which triangle
};


class Segment                           // segment(the group of boundary edges)
{
public:
  char      label[32];                  // the label of this segment
  int       bc_index;                   // boundary index
  int       interface;                  // -1 means a boundary segment,>=0 is the interface index
  int       zone_index;                 // this segment belongs to which zone
  int       g_index;                    // the index of this segment in GSegment structure
  int       edge_num;
  vector<Edge> edge_array;              // segment edges.
  vector<int>  node_array;              // segment nodes.
  void      set_node_array();           // convert edge_array to node_array here
  int       get_edge_bc_index(int, int);
};

class GSegment   :public  Segment       // global segment
{
public:
  int       local_index;                // the index in the zone
};



class Tri                               // triangle
{
public:
  int     node[3];                      // three nodes, local index
  int     g_node[3];                    // the global index of 3 nodes
  double  edge_len[3];                  // the length of 3 edges: a,b and c
  double  angle[3];                     // the degree of 3 angles: A, B abd C
  double  xc,yc;                        // the location of circle center
  double  d[3];                         // the distance from circle center to each edge: da, db and dc
  double  s[3];                         // partial area of each region seperated by da, db and dc
  int     bc[3];                        // the bc index of 3 edge
  double  area;                         // the area of triangle
  int     zone_index;                   // the zone index
  int     local_index;                  // the tri index in the zone
};



class VEdge                     // Voronoi edge
{
public:
  int        cell1;             //2 cell
  int        cell2;
  int        c1_nb;             //cell2 index in cell 1's nb_array
  int        c2_nb;             //cell1 index in cell 2's nb_array
  double     x1,y1;
  double     x2,y2;
  double     length;
  int        flag;
};

class VoronoiCell               //VoronoiCell
{
public:
  double        x,y;              //cell centre
  int           nb_num;
  int           *nb_array;        //it's neighbor nodes
  int           *inb_array;       //inverse index of neighbor nodes
  double        *elen;
  double        *ilen;
  double        *angle;
  double        area;
  int           *celledge;
  int           bc_index;         //the index to boundary of the node
  double        sa,sb,sc;         //for get gradian by least-squares method
  VoronoiCell();
  VoronoiCell(const VoronoiCell&);
  ~VoronoiCell();

};



#endif
