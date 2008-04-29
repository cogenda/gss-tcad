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

#ifndef _zone_h_
#define _zone_h_

#include "element.h"
#include <string>
using namespace std;

class DAVCELL
{
private:
  int cell_num;             // total number of elements
  VoronoiCell *cell_array;  // dynamic array
public:
  int size()                            { return cell_num;    }
  VoronoiCell *GetPointer (int i)       { return cell_array + i; }
  const VoronoiCell & operator[]  (int i) { return cell_array[i];  }
public:
  void Init(int c)                      { cell_num = c;  cell_array = new VoronoiCell[c]; }
  void clear()                          { delete [] cell_array; cell_num = 0; cell_array = 0;}
  DAVCELL ()                            { cell_num = 0; cell_array = 0; }
  ~DAVCELL ()                           { if(cell_array) delete [] cell_array; }
};



class ZONE
{
public:
  int                zone_index;                 //zone array index
  char               zonename[32];               //cgns zone name
  char               zonelabel[32];              //zone material label
  int                offset;                     //the offset of the node index in global index
  int                zone_conn_num;              //the zone connect number
  vector<string>     donor_zone_name;            //the array of dornor zone name
  vector<string>     connect_name;
  vector<int>        donor_zone_array;           //the array of dornor zone index

  double               field_area_vcell;           //field area summed  by vcell !may have problem
  double               field_area_tcell;           //field area summed  by tri cell
  vector<Node>       danode;
  vector<Tri>        datri;
  vector<Segment>    dasegment;                  //used for mesh refine
  vector<VEdge>      davedge;
  DAVCELL            davcell;

public:
  int       build_tricell_data();
  int       setup_voronoicell();
  int       build_voronoicell_data();
};

#endif
