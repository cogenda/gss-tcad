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
/*  Last update: Jan 01, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _interface_h_
#define _interface_h_
#include "zone.h"
#include <vector>
using namespace std;

class Interface         // a collection of interface edges
{
public:
  char conn_name[32];
  char zone1_name[32];
  char zone2_name[32];
  int  zone1;
  int  zone2;
  ZONE *pzone1;
  ZONE *pzone2;
  int  flag;
  int node_num;
  vector<int> index_array1;
  vector<int> index_array2;
  vector<int> gindex_array;
  int   Find_neighbor_zone_index(int z);
  int   Find_neighbor_node_index(int z,int n);
  ZONE *Get_neighbor_zone_pointer(int z);
  void  clear();
};


class ZoneInterface
{
public:
  vector<Interface> interface;
public:
  int  Find(const char* conn);
  Interface & operator[](int i)     { return interface[i];}
  int   size()                      { return interface.size(); }
  void  Clear_flag_all()            { for(int i=0;i<interface.size();i++) interface[i].flag = 0;}
  void  Set_flag_all()              { for(int i=0;i<interface.size();i++) interface[i].flag = 1;}
  void  Clear_flag(const char *conn){ interface[Find(conn)].flag = 0;}
  void  Set_flag(const char *conn)  { interface[Find(conn)].flag = 1;}
  int   Flag(const char *conn)      { return interface[Find(conn)].flag; }
  void  Clear_flag(int i)           { interface[i].flag = 0;}
  void  Set_flag(int i)             { interface[i].flag = 1;}
  int   Flag(int i)                 { return interface[i].flag; }
  int   IsBelong(int inf,int z);
  vector<int> * Get_node_array(int inf,int z);

};

#endif
