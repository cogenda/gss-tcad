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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <cgnslib.h>
#include "geom.h"
#include "mesh.h"
#include "log.h"

//simple error hanfle function
inline int cgns_error_handle(const char *s,int index_file)
{
  gss_log.string_buf()<<s<<ends;
  gss_log.record();
  cg_close(index_file);
  return 1;
}

/* ----------------------------------------------------------------------------
 * MESH::import_cgns:  This function read node location, cell node, bounary
 * definition and solution data from cgns file
 */
int MESH::import_cgns(const char *cgnsfile)
{
  int index_file,index_base,index_zone,index_sect,index_I;
  int nbase,ndescriptor,nsection,nbndry,nconns,npnts,ndata_donor,ElementDataSize,BCNum;
  char sectionname[32],descriptorname[32],boconame[32],solname[32],fieldname[32],connectname[32],donorname[32];

  int  cell_dim,physical_dim,isize[3],iparent_flag,iparentdata;
  int  normalindex,normallistflag,ndataset,nsol,normallist;
  double *x,*y;
  int *ielem;
  int  err;
  ZoneType_t     zonetype,donor_zonetype;
  ElementType_t  elemtype;
  BCType_t       bocotype;
  PointSetType_t ptset,donor_ptset_type;
  DataType_t     normaldatatype,donor_datatype;
  GridLocation_t location;
  GridConnectivityType_t connect_type;

  gss_log.string_buf()<<"\n--------------------------------------------------\n";  
  gss_log.record();
  
  //--------------------------------------------------------------------------------------
  gss_log.string_buf()<<"\nLoading cgns...\n";  
  gss_log.record();
  //open the cgns file
  if(cg_open(cgnsfile,MODE_READ,&index_file))
  {
    gss_log.string_buf()<<"error, can't open cgns file!\n";
    gss_log.record();
    return 1;
  }
  gss_log.string_buf()<<"cgns file     : "<<cgnsfile<<"\n";  
  gss_log.record();

  //--------------------------------------------------------------------------------------
  //cgns file should have only one base
  assert(!cg_nbases(index_file,&nbase));
  if(nbase!=1)
    return cgns_error_handle("multibases found,not a supported cgns file",index_file);

  //--------------------------------------------------------------------------------------
  //read base information
  index_base=1; //because only have one base
  assert(!cg_base_read(index_file,index_base,basename,&cell_dim,&physical_dim));
  if(physical_dim!=2)
    return cgns_error_handle("wrong physical_dim found,not a supported cgns file",index_file);
  gss_log.string_buf()<<"base name     : " << basename<<"\n";  
  gss_log.record();

  //--------------------------------------------------------------------------------------
  //cgns file may have several zones
  assert(!cg_nzones(index_file,1,&zone_num));
  gss_log.string_buf()<<"This cgns file have "<<zone_num<<" zones\n\n";  
  gss_log.record();
  zone = new ZONE[zone_num];
  int node_offset = 0;
  //read zone information
  for(index_zone=1;index_zone<=zone_num;index_zone++,node_offset+=isize[0])
  {
    ZONE *pzone = &zone[index_zone-1];
    pzone->zone_index = index_zone-1;

    assert(!cg_zone_type(index_file,index_base,index_zone,&zonetype));
    if(zonetype!=Unstructured)
      return cgns_error_handle("Unstructured zone needed,not a supported cgns file",index_file);
    assert(!cg_zone_read(index_file,index_base,index_zone,pzone->zonename,isize));
    gss_log.string_buf()<<"read zone "<<index_zone<<"\n";             gss_log.record();
    gss_log.string_buf()<<"zone name     : "<<pzone->zonename<<"\n";  gss_log.record();
    gss_log.string_buf()<<"points        : "<<isize[0]<<"\n";         gss_log.record();
    gss_log.string_buf()<<"elements      : "<<isize[1]<<"\n";         gss_log.record();

    // this block has problems
    char *label;
    cg_goto(index_file,index_base,"Zone_t",index_zone,"end");
    cg_ndescriptors(&ndescriptor);
    if(ndescriptor)
    {
      cg_descriptor_read(ndescriptor,descriptorname,&label);
      strcpy(pzone->zonelabel,label);
    }
    else
      strcpy(pzone->zonelabel,"Unspecified");
    free(label);
    // end of this block

    pzone->offset   = node_offset;
    //alloc double array for x and y coordinate
    x = new double[isize[0]]; assert(x);
    y = new double[isize[0]]; assert(y);

    //read the coordinate
    int index = 1;
    assert(!cg_coord_read(index_file,index_base,index_zone,"CoordinateX",RealDouble,&index,&isize[0],x));
    assert(!cg_coord_read(index_file,index_base,index_zone,"CoordinateY",RealDouble,&index,&isize[0],y));


    //assign coordinates to danode
    Node tmp_node;
    for(int i=0;i<isize[0];i++)
    {   //scale length to my own standerd
      tmp_node.x = x[i]*mesh_scale;
      tmp_node.y = y[i]*mesh_scale;
      tmp_node.bc_index = 0;
      tmp_node.zone_index = index_zone-1;
      tmp_node.g_index = -1;
      if (tmp_node.x > xMax) xMax = tmp_node.x;
      if (tmp_node.x < xMin) xMin = tmp_node.x;
      if (tmp_node.y > yMax) yMax = tmp_node.y;
      if (tmp_node.y < yMin) yMin = tmp_node.y;
      pzone->danode.push_back(tmp_node);
    }
    delete [] x;
    delete [] y;

    //--------------------------------------------------------------------------------------
    // find out how many sections
    assert(!cg_nsections(index_file,index_base,index_zone,&nsection));

    //read section information
    for(index_sect = 1; index_sect<=nsection; index_sect++)
    {
      assert(!cg_section_read(index_file,index_base,index_zone,index_sect,
                              sectionname,&elemtype,&index,&isize[1],&nbndry,&iparent_flag));
      if(elemtype==MIXED) break;
    }
    if(index_sect>nsection)
      return cgns_error_handle("error elem type,not a supported cgns file",index_file);
    //read elements
    assert(!cg_ElementDataSize(index_file,index_base,index_zone,index_sect,&ElementDataSize));

    ielem= new int[ElementDataSize]; assert(ielem);
    assert(!cg_elements_read(index_file,index_base,index_zone,index_sect,ielem,&iparentdata));
    //alloc memory for datri

    Tri tmp_Tri;
    int *pelem = ielem;
    for(int i=0;i<=isize[1]-index;i++)
    {
      if(*pelem == TRI_3)
      {
        pelem++;
        tmp_Tri.node[0] = *(pelem++)-1; //cgns node index begin with 1. sub it
        tmp_Tri.node[1] = *(pelem++)-1;
        tmp_Tri.node[2] = *(pelem++)-1;
        tmp_Tri.bc[0] = 0;
        tmp_Tri.bc[1] = 0;
        tmp_Tri.bc[2] = 0;
        tmp_Tri.zone_index = index_zone-1;
        tmp_Tri.local_index = i;
      }
      else
      {
        delete [] ielem;
        return cgns_error_handle("Unsupported mesh strucure, GSS reads triangular mesh only!",index_file);
      }

      pzone->datri.push_back(tmp_Tri);
    }

    //read zone connect information

    cg_nconns(index_file,index_base,index_zone,&nconns);
    pzone->zone_conn_num = nconns;
    pzone->donor_zone_array.clear();
    //pzone->donor_zone_name.resize(nconns);

    if(nconns>0)
    {
      Interface pitf;
      for(index_I=1;index_I<=nconns;index_I++)
      {
        cg_conn_info(index_file,index_base,index_zone,index_I,connectname ,
                     &location ,&connect_type ,
                     &ptset , &npnts ,donorname ,
                     &donor_zonetype , &donor_ptset_type ,
                     &donor_datatype , &ndata_donor );
        pzone->donor_zone_name.push_back(donorname);
        pzone->connect_name.push_back(connectname);
        if(zinterface.Find(connectname)==-1) //the interface is not set
        {
          strcpy(pitf.conn_name,connectname);
          strcpy(pitf.zone1_name,pzone->zonename);
          strcpy(pitf.zone2_name,donorname);
          pitf.flag = 0;
          pitf.node_num = npnts;
          int *index_array1 = new int[npnts];
          int *index_array2 = new int[npnts];

          cg_conn_read(index_file,index_base,index_zone,index_I,index_array1,
                       donor_datatype,index_array2);
          for(int j=0;j<pitf.node_num;j++)
          {
            pitf.index_array1.push_back(index_array1[j]-1);
            pitf.index_array2.push_back(index_array2[j]-1);
          }
          delete [] index_array1;
          delete [] index_array2;
          zinterface.interface.push_back(pitf);
          pitf.index_array1.clear();
          pitf.index_array2.clear();
        }
        gss_log.string_buf()<<"zone Connectivity "<<connectname<<" read "<<npnts<<" nodes\n";
        gss_log.record();
      }

    }

    //---------------------------------------------------------------------------------------
    //now read ZoneBC
    //the boundary marker will be assigned to each boundary node

    //how many BC in the ZoneBC
    assert(!cg_nbocos(index_file,index_base,index_zone,&BCNum));
    gss_log.string_buf()<<"BC number     : "<<BCNum<<"\n";
    gss_log.record();
    Segment tmp_segment;
    for(int BCIndex=1;BCIndex<=BCNum;BCIndex++)
    {
      assert(!cg_boco_info(index_file,index_base,index_zone,BCIndex,boconame,
                           &bocotype,&ptset,&npnts,&normalindex,&normallistflag,&normaldatatype,&ndataset));
      int start,end;

      //read section information about boundary segment
      for(index_sect = 1; index_sect<=nsection; index_sect++)
      {
        assert(!cg_section_read(index_file,index_base,index_zone,index_sect,
                                sectionname,&elemtype,&start,&end,&nbndry,&iparent_flag));
        if(!strcmp(boconame,sectionname)&&elemtype==BAR_2) break;
      }
      if(index_sect>nsection)
        return cgns_error_handle("boundary elem error,not a supported cgns file",index_file);
      //read elements
      assert(!cg_ElementDataSize(index_file,index_base,index_zone,index_sect,&ElementDataSize));
      assert(!cg_elements_read(index_file,index_base,index_zone,index_sect,ielem,&iparentdata));


      int bc_elem_range[2];
      cg_boco_read(index_file,index_base,index_zone,BCIndex,bc_elem_range,&normallist);
      if(bc_elem_range[0]!=start||bc_elem_range[1]!=end)
        return cgns_error_handle("boundary elem error,not a supported cgns file",index_file);

      tmp_segment.bc_index = BCIndex+bc_num;
      strcpy(tmp_segment.label,boconame);
      tmp_segment.edge_num = end-start+1;
      tmp_segment.interface=-1;
      tmp_segment.zone_index=index_zone-1;
      tmp_segment.g_index = -1;
      Edge tmp_edge;
      for(int i=0;i<end-start+1;i++)
      {
        tmp_edge.p1 = ielem[2*i]-1;
        tmp_edge.p2 = ielem[2*i+1]-1;
        tmp_edge.bc_index = BCIndex+bc_num;
        tmp_segment.edge_array.push_back(tmp_edge);
        pzone->danode[ielem[2*i]-1].bc_index   = BCIndex + bc_num;
        pzone->danode[ielem[2*i+1]-1].bc_index = BCIndex + bc_num;
      }
      pzone->dasegment.push_back(tmp_segment);
      tmp_segment.edge_array.clear();
      gss_log.string_buf()<<"   BC "<<BCIndex<<" "<<boconame<<"... "<<end-start+1<<" edges read\n";
      gss_log.record();

    }
    bc_num += BCNum;  //the zone bc_index must be different in global
    delete [] ielem;

    gss_log.string_buf()<<"zone "<<index_zone<<" read finished\n\n";
    gss_log.record();
  }

  //post process
  for(int i=0;i<zinterface.interface.size();i++)
  {
    zinterface[i].zone1 = Get_zone(zinterface[i].zone1_name)->zone_index;
    zinterface[i].zone2 = Get_zone(zinterface[i].zone2_name)->zone_index;
    zinterface[i].pzone1 = &zone[zinterface[i].zone1];
    zinterface[i].pzone2 = &zone[zinterface[i].zone2];
  }
  for(int i=0;i<zone_num;i++)
    for(int j=0;j<zone[i].zone_conn_num;j++)
      zone[i].donor_zone_array.push_back(Get_zone(zone[i].donor_zone_name[j].c_str())->zone_index);

  for(int z=0;z<zone_num;z++)
    for(int j=0;j<zone[z].dasegment.size();j++)
      for(int k=0;k<zinterface.interface.size();k++)
      {
        if(zinterface.IsBelong(k,z))
        {
          if(!strcmp(zone[z].dasegment[j].label,zinterface.interface[k].conn_name))
              zone[z].dasegment[j].interface=k;
        }
      }

  for(int i=0;i<zone_num;i++)
    for(int j=0;j<zone[i].dasegment.size();j++)
    {
      zone[i].dasegment[j].set_node_array();
    }

  gss_log.string_buf()<<"mesh read finished\n\n";
  gss_log.record();
  cg_close(index_file);
  return 0;

}



/* ----------------------------------------------------------------------------
 * MESH::build_mesh:  an interface function to call setup_voronoicell()
 * build_voronoicell_data() and build_ghost_cell() for each zone
 */
void  MESH::build_zone()
{
  for(int i=0;i<zone_num;i++)
  {
    zone[i].build_tricell_data();
    zone[i].setup_voronoicell();
    zone[i].build_voronoicell_data();
  }
}



/* ----------------------------------------------------------------------------
 * MESH::reorder_zone:  This function reorder the nodes in the specified zone .
 * for reduce the bandwidth of matrix.
 * danodes,datri,dasegment and correlative interface may be modified.
 * it calls setup_voronoicell() first to get adjacent information .
 * when it was done, voronoicell should be rebuilt.
 * it returns a pointer to an array which contain convertor index
 * used as new_index = old_to_new[old_index],this array must be deleted by hand
 */
int * MESH::reorder_zone(int zone_index)
{

  ZONE *pzone = &zone[zone_index];
  pzone->setup_voronoicell();
  //find the zone node which has the minimal x and y coordinates.
  int start = 0;
  for(int i=0;i<pzone->danode.size();i++)
    if( pzone->danode[i].x<=pzone->danode[start].x && pzone->danode[i].y<=pzone->danode[start].y)
      start = i;
    else if(pzone->danode[i].x<pzone->danode[start].x)
      start = i;

  //ok ,begin at this node. do Breadth-First Search
  int *af = new int[pzone->danode.size()];
  int *new_order = new int[pzone->danode.size()];
  int *old_to_new = new int[pzone->danode.size()];
  int index = 0;
  for(int i=0;i<pzone->danode.size();i++) af[i] = 0;
  queue<int> Q;
  af[start] =1; Q.push(start);
  while(!Q.empty())
  {
    int current = Q.front();   Q.pop();
    new_order[index] = current; //record the order
    old_to_new[current] = index;
    index++;
    for(int j=0;j<pzone->davcell[current].nb_num;j++)
      if(!af[pzone->davcell[current].nb_array[j]])
      {
        Q.push(pzone->davcell[current].nb_array[j]);
        af[pzone->davcell[current].nb_array[j]] = 1;
      }
  }
  //we have got the new order. do post process
  vector<Node>   tmpdanode = pzone->danode;
  for(int i=0;i<pzone->danode.size();i++)
  {
    pzone->danode[i] = tmpdanode[new_order[i]];
  }

  //the reorder process will not change the clock wise order.
  for(int i=0;i<pzone->datri.size();i++)
  {
    pzone->datri[i].node[0] = old_to_new[pzone->datri[i].node[0]];
    pzone->datri[i].node[1] = old_to_new[pzone->datri[i].node[1]];
    pzone->datri[i].node[2] = old_to_new[pzone->datri[i].node[2]];
  }

  for(int i=0;i<pzone->dasegment.size();i++)
  {
    for(int j=0;j<pzone->dasegment[i].edge_array.size();j++)
    {
      pzone->dasegment[i].edge_array[j].p1 = old_to_new[pzone->dasegment[i].edge_array[j].p1];
      pzone->dasegment[i].edge_array[j].p2 = old_to_new[pzone->dasegment[i].edge_array[j].p2];
    }
    pzone->dasegment[i].set_node_array();

    if(pzone->dasegment[i].interface>=0)
    {
      int inf = pzone->dasegment[i].interface;
      vector<int> *pinf_node = zinterface.Get_node_array(inf,zone_index);
      vector<int>:: iterator pint = pinf_node->begin();
      for(int j=0; j< pinf_node->size(); j++,pint++)
        *pint = old_to_new[*pint];
    }
  }

  delete [] af;
  delete [] new_order;
  return old_to_new;

}


/* ----------------------------------------------------------------------------
 * MESH::reorder_mesh:  interface function to call reorder_zone
 *
 */
void   MESH::reorder_mesh()
{

  for(int i=0;i<zone_num;i++)
    delete reorder_zone(i);
}



/* ----------------------------------------------------------------------------
 * MESH::Get_zone:  This function get zone pointer by
 * zone name
 */
ZONE*  MESH::Get_zone(const char *name)
{
  for(int i=0;i<zone_num;i++)
    if(!strcmp(zone[i].zonename,name)) return &zone[i];
  return (ZONE*)NULL;
}

/* ----------------------------------------------------------------------------
 * MESH::Get_zone_index:  This function get zone index by
 * zone name
 */
int  MESH::Get_zone_index(const char *name)
{
  for(int i=0;i<zone_num;i++)
    if(!strcmp(zone[i].zonename,name)) return i;
  return -1;
}

/* ----------------------------------------------------------------------------
 * MESH::Get_zone:  This function get zone pointer by
 * zone index
 */
ZONE*  MESH::Get_zone(const int z)
{
  for(int i=0;i<zone_num;i++)
    if(zone[i].zone_index == z) return &zone[i];
  return (ZONE*)NULL;
}

/* ----------------------------------------------------------------------------
 * MESH::Get_segment:  This function get zone segment pointer by
 * global bc_index(from 1 to bc_num)
 */
Segment* MESH::Get_segment(const int bc)
{
  for(int i=0;i<zone_num;i++)
    for(int j=0;j<zone[i].dasegment.size();j++)
      if(bc == zone[i].dasegment[j].bc_index)
        return &zone[i].dasegment[j];

  return        (Segment*)NULL;
}

/* ----------------------------------------------------------------------------
 * MESH::zone_to_field:  This function gather local zone mesh into a global
 * mesh structrue. for adaptive mesh refinement
 */

int  MESH::zone_to_field()
{
  gnode.clear();
  gtri.clear();
  gsegment.clear();
  int gnode_num = 0;
  GNode node;
  int total_node_num = 0;
  for(int i=0;i<zone_num;i++)
    total_node_num+=zone[i].danode.size();
  int *sign = new int [total_node_num];
  int *reorder = new int [total_node_num];
  int *p = sign;
  zinterface.Clear_flag_all();
  //eat duplicate interface node.
  for(int i=0;i<zone_num;i++)
  {
    for(int j=0;j<zone[i].danode.size();j++)
    {
      sign[zone[i].offset + j] = zone[i].offset + j;
    }
    for(int j=0;j<zone[i].connect_name.size();j++)
    {
      int z2 = zone[i].donor_zone_array[j];
      const char *connect_name = zone[i].connect_name[j].c_str();       
      if(zinterface.Flag(connect_name)==0)
        zinterface.Set_flag(connect_name);
      else
      {
        int inf =  zinterface.Find(connect_name);
        vector<int>* pinf_node1 = zinterface.Get_node_array(inf,i);
        vector<int>* pinf_node2 = zinterface.Get_node_array(inf,z2);

        for(int k=0;k<zinterface[inf].node_num;k++)
          sign[zone[z2].offset + pinf_node2->at(k)] = sign[zone[i].offset + pinf_node1->at(k)];
      }
    }
  }

  for(int i=0;i<zone_num;i++)
    for(int j=0;j<zone[i].danode.size();j++)
      if(sign[zone[i].offset + j]==zone[i].offset + j)
      {
        node.x = zone[i].danode[j].x;
        node.y = zone[i].danode[j].y;
        node.bc_index = zone[i].danode[j].bc_index;
        node.local_index = j;
        node.zone_index = i;
        node.g_index = gnode.size();
        zone[i].danode[j].g_index = gnode.size();
        gnode.push_back(node);
        reorder[zone[i].offset + j] = gnode_num++;
      }

  //set g_index of duplicate node.
  zinterface.Clear_flag_all();
  for(int i=0;i<zone_num;i++)
  {
    for(int j=0;j<zone[i].connect_name.size();j++)
    {
      int z2 = zone[i].donor_zone_array[j];
      const char *connect_name = zone[i].connect_name[j].c_str();       
      if(zinterface.Flag(connect_name)==0)
        zinterface.Set_flag(connect_name);
      else
      {
        int inf =  zinterface.Find(connect_name);
        vector<int>* pinf_node1 = zinterface.Get_node_array(inf,i);
        vector<int>* pinf_node2 = zinterface.Get_node_array(inf,z2);

        for(int k=0;k<zinterface[inf].node_num;k++)
          zone[z2].danode[pinf_node2->at(k)].g_index = zone[i].danode[pinf_node1->at(k)].g_index;
      }
    }
  }

  // set global triangles
  for(int i=0;i<zone_num;i++)
    for(int j=0;j<zone[i].datri.size();j++)
    {
      Tri tmp_tri = zone[i].datri[j];
      tmp_tri.g_node[0] = reorder[sign[zone[i].offset + tmp_tri.node[0]]];
      tmp_tri.g_node[1] = reorder[sign[zone[i].offset + tmp_tri.node[1]]];
      tmp_tri.g_node[2] = reorder[sign[zone[i].offset + tmp_tri.node[2]]];

      gtri.push_back(tmp_tri);
    }

  //I droped duplicated zone interface segment. so the bc_index is not continuous.
  //but it dosen't matter for getting enough information to build a Planar Straight Line Graph (PSLG) .
  GSegment tmp_segment;
  for(int i=0;i<zone_num;i++)
  {
    for(int j=0;j<zone[i].dasegment.size();j++)
    {
      int flag = -1;
      Segment zonesegment = zone[i].dasegment[j];
      //if two zone segmants have same name, skip one
      for(int k=0;k<gsegment.size();k++)
        if(!strcmp(gsegment[k].label,zonesegment.label))  flag = k;
      if(flag!=-1)
      {
        zone[i].dasegment[j].g_index = flag;
        continue;
      }
      //else copy it to global level
      tmp_segment.bc_index = zonesegment.bc_index;
      tmp_segment.interface = zonesegment.interface;
      strcpy(tmp_segment.label,zonesegment.label);
      tmp_segment.edge_num = zonesegment.edge_num;
      tmp_segment.g_index = gsegment.size();
      zone[i].dasegment[j].g_index = gsegment.size();
      tmp_segment.zone_index = i;
      tmp_segment.local_index = j;
      for(int k=0;k<zonesegment.edge_num;k++)
      {
        Edge tmp_edge = zonesegment.edge_array[k];
        tmp_edge.p1 = reorder[sign[zone[i].offset+zonesegment.edge_array[k].p1]];
        tmp_edge.p2 = reorder[sign[zone[i].offset+zonesegment.edge_array[k].p2]];
        tmp_edge.bc_index = zonesegment.edge_array[k].bc_index;
        tmp_segment.edge_array.push_back(tmp_edge);
      }
      gsegment.push_back(tmp_segment);
      tmp_segment.edge_array.clear();
    }

  }

  for(int i=0;i<gsegment.size();i++)
  {
    gsegment[i].set_node_array();

    for(int h=0; h<gsegment[i].edge_array.size();h++)
    {
            Edge tmp_edge = gsegment[i].edge_array[h];
            tmp_edge.p1 = gsegment[i].edge_array[h].p1;
            tmp_edge.p2 = gsegment[i].edge_array[h].p2;
    }
  }
  delete [] sign;
  delete [] reorder;
  return 0;
}


/* ----------------------------------------------------------------------------
 * MESH::field_to_zone:  This function build local zone mesh from  global
 * mesh structrue.
 */

int  MESH::field_to_zone()
{

  //the nodes will be classified by its zone index and insert to each zone
  int **local_index = new int*[zone_num];       //zone_index = local_index[reg][global_index]
  int *zonenodes = new int[zone_num];

  for(int i=0; i<zone_num; i++)
  {
    local_index[i] = new int[gnode.size()];
  }

  int *af = new int[gnode.size()];

  for(int i=0;i<zone_num;i++)
  {
    int *pf = af;
    for(int j = 0; j < gnode.size(); j++)
    {
      local_index[i][j] = -1; //-1 means node not in this zone
      *pf++ = 0;        //visited flag array
    }
    zonenodes[i] = 0;

    for(int j = 0; j < gtri.size(); j++)     //search for all the nodes
      if (gtri[j].zone_index == i)
      {
        if (!af[gtri[j].g_node[0]]) { af[gtri[j].g_node[0]] = 1; local_index[i][gtri[j].g_node[0]] = zonenodes[i]++; }
        if (!af[gtri[j].g_node[1]]) { af[gtri[j].g_node[1]] = 1; local_index[i][gtri[j].g_node[1]] = zonenodes[i]++; }
        if (!af[gtri[j].g_node[2]]) { af[gtri[j].g_node[2]] = 1; local_index[i][gtri[j].g_node[2]] = zonenodes[i]++; }
      }
  }
  delete [] af;

  int node_offset = 0;
  zinterface.Clear_flag_all();
  for(int i=0; i<zone_num; node_offset+=zonenodes[i],i++)
  {
    zone[i].danode.resize(zonenodes[i]);
    for(int j=0;j<gnode.size();j++)
      if(local_index[i][j]>-1)  //process node in this zone
      {
        gnode[j].zone_index = i;
        gnode[j].local_index = local_index[i][j];
        zone[i].danode[local_index[i][j]] = gnode[j];
        zone[i].danode[local_index[i][j]].g_index = j;
      }

    zone[i].offset = node_offset;

    zone[i].datri.clear();
    for(int j=0;j<gtri.size();j++)
      if (gtri[j].zone_index == i)  //is the triangle belongs to this zone?
      {
        Tri tmp_tri = gtri[j];
        tmp_tri.node[0] = local_index[i][gtri[j].g_node[0]];
        tmp_tri.node[1] = local_index[i][gtri[j].g_node[1]];
        tmp_tri.node[2] = local_index[i][gtri[j].g_node[2]];
        tmp_tri.local_index = zone[i].datri.size();
        zone[i].datri.push_back(tmp_tri);
      }


    //process boundary segment here
    for(int j=0;j<zone[i].dasegment.size();j++)
    {
      zone[i].dasegment[j].edge_num = 0;
      zone[i].dasegment[j].edge_array.clear();
    }

    for(int j=0;j<zone[i].dasegment.size();j++)
      for(int k=0;k<gsegment.size();k++)
        if(!strcmp(zone[i].dasegment[j].label,gsegment[k].label))
        {
          zone[i].dasegment[j].edge_num = gsegment[k].edge_num;
          for(int h=0; h<gsegment[k].edge_array.size();h++)
          {
            Edge tmp_edge = gsegment[k].edge_array[h];
            tmp_edge.p1 = local_index[i][gsegment[k].edge_array[h].p1];
            tmp_edge.p2 = local_index[i][gsegment[k].edge_array[h].p2];
            zone[i].dasegment[j].edge_array.push_back(tmp_edge);
          }
        }

    for(int j=0;j<zone[i].dasegment.size();j++)
    {
      vector<int> *pinfi,*pinfj;
      zone[i].dasegment[j].set_node_array();
      if(zone[i].dasegment[j].interface!=-1)
      {
        int inf = zone[i].dasegment[j].interface;
        zinterface[inf].clear();
        pinfi =  zinterface.Get_node_array(inf,i);
        int zn = zinterface[inf].Find_neighbor_zone_index(i);
        pinfj =  zinterface.Get_node_array(inf,zn);
        zinterface[inf].node_num = zone[i].dasegment[j].node_array.size();
        for(int h=0; h<zone[i].dasegment[j].node_array.size();h++)
        {
          int g_index = zone[i].danode[zone[i].dasegment[j].node_array[h]].g_index;
          pinfi->push_back(local_index[i][g_index]);
          pinfj->push_back(local_index[zn][g_index]);
        }
      }
    }

    //reset the bc_index for duplicate interface node
    for(int j=0;j< zone[i].dasegment.size();j++)
      for(int k=0;k<zone[i].dasegment[j].node_array.size();k++)
      {
        zone[i].danode[zone[i].dasegment[j].node_array[k]].bc_index = zone[i].dasegment[j].bc_index;
      }

  }

  // all things are done
  for(int i =0; i<zone_num; i++)
  {
    delete [] local_index[i];
  }
  delete  [] local_index;
   delete  [] zonenodes ;
  return 0;

}



/* ----------------------------------------------------------------------------
 * MESH::export_cgns:  This function write mesh to cgns file.
 */
void MESH::export_mesh(const char *cgnsfile)
{
  int fn,B,Z,C,S,BC,I,SOL,F;
  int size[3];
  double *x,*y;
  int *elem;
  char bcname[32];

  /* remove old file if exist*/
  remove(cgnsfile);

  /* open CGNS file for write */
  cg_open(cgnsfile,MODE_WRITE,&fn);

  /* create base (can give any name) */
  cg_base_write(fn,"GSS_CoreFile",2,2,&B); /*two dimension*/

  /* create zone */
  for(int i=0;i<zone_num;i++)
  {

    size[0] = zone[i].danode.size();
    size[1] = zone[i].datri.size();
    size[2] = 0;

    x = new double[size[0]];
    y = new double[size[0]];
    elem= new int[5*size[1]];
    //set zone
    cg_zone_write(fn,B,zone[i].zonename,size,Unstructured,&Z);
    cg_goto(fn,B,"Zone_t",Z,"end");
    cg_descriptor_write("RegionType",zone[i].zonelabel);

    // write grid coordinates
    for(int j=0;j<zone[i].danode.size();j++)
    {
      //convert the unit
      x[j] = zone[i].danode[j].x/mesh_scale;
      y[j] = zone[i].danode[j].y/mesh_scale;
    }

    cg_coord_write(fn,B,Z,RealDouble,"CoordinateX",x,&C);
    cg_coord_write(fn,B,Z,RealDouble,"CoordinateY",y,&C);

    // set element connectivity here
    int *pelem=elem;


    for(int j = 0; j < zone[i].datri.size(); j++)
    {
        *pelem++ = TRI_3;
        *pelem++ = zone[i].datri[j].node[0]+1;
        *pelem++ = zone[i].datri[j].node[1]+1;
        *pelem++ = zone[i].datri[j].node[2]+1;
    }

    cg_section_write(fn,B,Z,"GridElements",MIXED,1,size[1],0,elem,&S);

    delete [] x;
    delete [] y;
    delete [] elem;

    // write boundary segment

    int *asegment = new int[2*zone[i].danode.size()];
    int  segrange[2];
    segrange[0] = size[1]+1;
    int  csegment;
    for(int j=0;j<zone[i].dasegment.size();j++)
    {
      csegment = 0;
      int *psegment = asegment;

      for(int k=0; k<zone[i].dasegment[j].edge_num; k++)
      {
        *psegment++ = zone[i].dasegment[j].edge_array[k].p1+1;
        *psegment++ = zone[i].dasegment[j].edge_array[k].p2+1;
        csegment++;
      }

      if(csegment>0)
      {
        segrange[1] = segrange[0] + csegment - 1;
        cg_section_write(fn,B,Z,zone[i].dasegment[j].label,BAR_2,segrange[0],segrange[1],0,asegment,&S);
        cg_boco_write(fn,B,Z,zone[i].dasegment[j].label,BCTypeUserDefined,ElementRange,2,segrange,&BC);
        segrange[0] = segrange[1]+1;
      }

    }
    delete [] asegment;

    //write zone to zone connectivity data if necessary
    if(zone_num>1)
    {   //build connectivity data here
      int inf;
      for(int j=0;j<zone[i].connect_name.size();j++)
        if( (inf = zinterface.Find(zone[i].connect_name[j].c_str())) > -1)
        {
          int z2 = zone[i].donor_zone_array[j];
          vector<int> *p1 = zinterface.Get_node_array(inf,i);
          vector<int> *p2 = zinterface.Get_node_array(inf,z2);
          int *conn_local = new int[p1->size()];
          int *conn_donor = new int[p2->size()];
          for(int k=0; k<p1->size();k++)
          {
            conn_local[k]=p1->at(k)+1;
            conn_donor[k]=p2->at(k)+1;
          }
          cg_conn_write(fn,B,Z,zone[i].connect_name[j].c_str(),Vertex, Abutting1to1,
                        PointList, p1->size(), conn_local,zone[z2].zonename,Unstructured,
                        PointListDonor,Integer,p2->size(),conn_donor ,&I);
          delete        [] conn_local;
          delete        [] conn_donor;
        }

    }

  }
  // close CGNS file
  cg_close(fn);
}

/* ----------------------------------------------------------------------------
 * MESH::clear_mesh:  This function delete mesh structure.
 */
void    MESH::clear_mesh()
{
  gnode.clear();
  gtri.clear();
  gsegment.clear();
  delete [] zone;
  zone = 0;
  zinterface.interface.clear();
  bc_num = 0;
  zone_num =0;

}

