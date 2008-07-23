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
/*  Last update: Jan 01, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include <string.h>
#include "interface.h"


int Interface::Find_neighbor_zone_index(int z)
{
  if (z==zone1) return zone2;
  if (z==zone2) return zone1;
  return -1;
}

ZONE * Interface::Get_neighbor_zone_pointer(int z)
{
  if (z==zone1) return pzone2;
  if (z==zone2) return pzone1;
  return (ZONE *)0;
}

int Interface::Find_neighbor_node_index(int z, int n)
{
  if(z==zone1)
  {
    for(int i=0;i<index_array1.size();i++)
      if(index_array1[i]==n) return index_array2[i];
    return -1;
  }
  else
  {
    for(int i=0;i<index_array2.size();i++)
      if(index_array2[i]==n) return index_array1[i];
    return -1;
  }
}

void  Interface::clear()
{
  node_num = 0;
  index_array1.clear();
  index_array2.clear();
  gindex_array.clear();
}
//------------------------------------------------------------------

int  ZoneInterface::Find(const char *conn)
{
  if(interface.size()==0) return -1;
  for(int i=0;i<interface.size();i++)
  {
    if(!strcmp(interface[i].conn_name,conn))
      return i;
  }
  return -1;
}


int   ZoneInterface::IsBelong(int inf,int z)
{
  if(interface[inf].zone1 == z) return 1;
  else if(interface[inf].zone2 == z) return 1;
  return 0;
}

vector<int> * ZoneInterface::Get_node_array(int inf,int z)
{
  if(interface[inf].zone1 == z)
    return &interface[inf].index_array1;
  else if(interface[inf].zone2 == z)
    return &interface[inf].index_array2;
  return   (vector<int> *)0; //prevent warning of complier
}
