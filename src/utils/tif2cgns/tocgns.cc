#include <math.h>
#include <stdio.h>
#include "cgnslib.h"
#include "tifdata.h"
#include "matdefine.h"

const int Electrode     = 1001;
const int Interface     = 2001;

void to_cgns(const char *filename)
{
  int fn,B,Z,C,S,BC,I,SOL,F;
  int size[3];
  double *x,*y;
  int *elem;
  char mat_name[32];
  char bcname[32];
  //this is the main part of the work. the nodes will be classified by
  //its region and insert to each zone

  size_t  region_num = region_array.size();	 //total region;
  int **local_index = new int*[region_num];    //reg_index = local_index[reg][global_index]
  int *conn_nodes = new int[region_num];         //for zone to zone connectivity node num
  int **conn = new int*[region_num];           //for zone to zone connectivity information

  for(size_t i=0; i<region_num; i ++)
  {
    local_index[i] = new int[node_array.size()];
    conn[i] = new int[node_array.size()];
  }

  int *af = new int[node_array.size()];

  for(size_t i=0;i<region_num;i++)
  {
    for(size_t j = 0; j < node_array.size(); j++) local_index[i][j] = -1; //-1 means node not in this region
    int *pf = af;
    for(size_t j = 0; j < node_array.size(); j++) *pf++ = 0;        //visited flag array
    region_array[i].node_num = 0;
    region_array[i].tri_num = 0;
    conn_nodes[i] = 0;

    for(size_t j = 0; j < tri_array.size(); j++)     //search for all the nodes
      if (tri_array[j].region == i)
      {
        region_array[i].tri_num++;
        int nodeA = tri_array[j].c1;
        int nodeB = tri_array[j].c2;
        int nodeC = tri_array[j].c3;
        if (!af[nodeA]) { af[nodeA] = 1;local_index[i][nodeA] = region_array[i].node_num++; }
        if (!af[nodeB]) { af[nodeB] = 1;local_index[i][nodeB] = region_array[i].node_num++; }
        if (!af[nodeC]) { af[nodeC] = 1;local_index[i][nodeC] = region_array[i].node_num++; }
      }
  }
  delete [] af;

  /* remove old file if exist*/
  remove(filename);

  /* open CGNS file for write */
  cg_open(filename,MODE_WRITE,&fn);

  /* create base (can give any name) */
  cg_base_write(fn,"Medici_To_GSS",2,2,&B); /*two dimension*/

  /* create zone */
  for(size_t i=0;i<region_num;i++)
  {

    size[0] = region_array[i].node_num;
    size[1] = region_array[i].tri_num;
    size[2] = 0;

    x = new double[size[0]];
    y = new double[size[0]];
    elem= new int[4*size[1]];
    //zone name set as region name
    cg_zone_write(fn,B,region_array[i].name,size,Unstructured,&Z);

    cg_goto(fn,B,"Zone_t",Z,"end");
    cg_descriptor_write("RegionType",FormatMaterialString(region_array[i].type));
    // write grid coordinates
    for(size_t j=0;j<node_array.size();j++)
      if(local_index[i][j]>-1)  //process node in this region
      {
        //convert the unit to cm
        x[local_index[i][j]] = node_array[j].x*1e-4;
        y[local_index[i][j]] = node_array[j].y*1e-4;
      }

    cg_coord_write(fn,B,Z,RealDouble,"CoordinateX",x,&C);
    cg_coord_write(fn,B,Z,RealDouble,"CoordinateY",y,&C);

    // set element connectivity here

    int *pelem=elem;
    for(size_t j = 0; j < tri_array.size(); j++)
      if (tri_array[j].region == i)  //is the triangle belongs to this region?
      {
        //the triangle in TIF is counter clock wise,
        //but I inversed the y axis befor. i have to reorder it from c3 to c1.
        *pelem++ = TRI_3;
        *pelem++ = local_index[i][tri_array[j].c3]+1;
        *pelem++ = local_index[i][tri_array[j].c2]+1;
        *pelem++ = local_index[i][tri_array[j].c1]+1;
      }

    cg_section_write(fn,B,Z,"GridElements",MIXED,1,size[1],0,elem,&S);



    // write Electrode boundary segment
    // for each interface
    int *asegment = new int[2*node_array.size()];
    int *psegment = asegment;
    int  segrange[2];
    segrange[0] = size[1]+1;
    int  csegment;
    for(size_t j=0;j<interface_array.size();j++)
    {
      if(interface_array[j].region<=0) //an electrode boundary
      {
        psegment = asegment;
        csegment = 0;
        for(int k=0;k<interface_array[j].boundary.size();k++)
        {
          int segment = interface_array[j].boundary[k];
          edge_array[segment].bc_type = Electrode;
          int nodeA = edge_array[segment].point1;
          int nodeB = edge_array[segment].point2;
          //only record the nodes that belong to this region
          if(local_index[i][nodeA]==-1||local_index[i][nodeB]==-1 )continue;

          *psegment++ = local_index[i][nodeA]+1;
          *psegment++ = local_index[i][nodeB]+1;
          csegment++;

        }
        if(csegment>0)
        {
          segrange[1] = segrange[0] + csegment - 1;
          cg_section_write(fn,B,Z,interface_array[j].name,BAR_2,
                           segrange[0],segrange[1],0,asegment,&S);
          cg_boco_write(fn,B,Z,interface_array[j].name,BCTypeUserDefined,
                        ElementRange,2,segrange,&BC);
          segrange[0] = segrange[1]+1;
        }
      }

    }


    //write zone to zone connectivity data if necessary
    if(region_num>1)
    {	//build connectivity data here
      for(size_t k=0; k<region_num; k++)
        for(size_t j=0; j<node_array.size(); j++)
        {
          if(k==i) continue;
          if(local_index[k][j]!=-1 &&local_index[i][j]!=-1)
          {
            conn[k][conn_nodes[k]++]=local_index[k][j]+1;
            conn[i][conn_nodes[i]++]=local_index[i][j]+1;
          }
        }
    }

    char conn_name[32];
    int  *conn_i=conn[i];
    for(size_t j=0; j<region_num; j++)
      if(j!=i&&conn_nodes[j]>0)
      {	//the conn_name is alpha ordered
        if(strcmp(region_array[i].name,region_array[j].name)<0)
          sprintf(conn_name,"IF_%s_to_%s",region_array[i].name,region_array[j].name);
        else
          sprintf(conn_name,"IF_%s_to_%s",region_array[j].name,region_array[i].name);
        cg_conn_write(fn,B,Z,conn_name,Vertex, Abutting1to1,
                      PointList, conn_nodes[j], conn_i,region_array[j].name,Unstructured,
                      PointListDonor,Integer,conn_nodes[j], conn[j],&I);

        conn_i+=conn_nodes[j];
      }

    for(size_t j=0; j<region_num; j++)
      conn_nodes[j]=0;

    // write boundary segment(Interface boundary) this boundary is not needed for
    // later computing, just for compatible with SGFramework

    for(size_t j=0;j<region_num;j++)
    {
      if(j!=i)
      {
        csegment = 0;
        psegment = asegment;

        for(size_t k=0; k<region_array[j].boundary.size(); k++)
        {
          int nodeA = edge_array[region_array[j].boundary[k]].point1;
          int nodeB = edge_array[region_array[j].boundary[k]].point2;
          //only record the nodes that belong to this region
          if(local_index[i][nodeA]==-1||local_index[i][nodeB]==-1 )continue;

          edge_array[region_array[j].boundary[k]].bc_type = Interface;
          *psegment++ = local_index[i][nodeA]+1;
          *psegment++ = local_index[i][nodeB]+1;
          csegment++;
        }

        if(csegment>0)
        {
          //the name is alpha ordered
          if(strcmp(region_array[i].name,region_array[j].name)<0)
            sprintf(bcname,"IF_%s_to_%s",region_array[i].name,region_array[j].name);
          else
            sprintf(bcname,"IF_%s_to_%s",region_array[j].name,region_array[i].name);

          segrange[1] = segrange[0] + csegment - 1;
          cg_section_write(fn,B,Z,bcname,BAR_2,
                           segrange[0],segrange[1],0,asegment,&S);
          cg_boco_write(fn,B,Z,bcname,BCTypeUserDefined,ElementRange,2,segrange,&BC);
          segrange[0] = segrange[1]+1;
	  continue;
        }
	
	//the fucking TSPUREM may not write interface segments for both side of regions.
	//I spent many hours for this problem. here is the fix
	for(size_t k=0; k<region_array[i].boundary.size(); k++)
        {
          int nodeA = edge_array[region_array[i].boundary[k]].point1;
          int nodeB = edge_array[region_array[i].boundary[k]].point2;
          //only record the nodes that belong to this region
          if(local_index[j][nodeA]==-1||local_index[j][nodeB]==-1 )continue;
          edge_array[region_array[i].boundary[k]].bc_type = Interface;
          *psegment++ = local_index[i][nodeA]+1;
          *psegment++ = local_index[i][nodeB]+1;
          csegment++;
        }

        if(csegment>0)
        {
          //the name is alpha ordered
          if(strcmp(region_array[i].name,region_array[j].name)<0)
            sprintf(bcname,"IF_%s_to_%s",region_array[i].name,region_array[j].name);
          else
            sprintf(bcname,"IF_%s_to_%s",region_array[j].name,region_array[i].name);

          segrange[1] = segrange[0] + csegment - 1;
          cg_section_write(fn,B,Z,bcname,BAR_2,
                           segrange[0],segrange[1],0,asegment,&S);
          cg_boco_write(fn,B,Z,bcname,BCTypeUserDefined,ElementRange,2,segrange,&BC);
          segrange[0] = segrange[1]+1;
        }  
      }
    }

    // write boundary segment(Neumann boundary)
    // for each boundary
    psegment = asegment;
    csegment = 0;

    for(size_t j=0;j<region_array[i].boundary.size();j++)
    {
      int segment = region_array[i].boundary[j];
      int nodeA = edge_array[segment].point1;
      int nodeB = edge_array[segment].point2;
      //only record the nodes that belong to this region
      if(local_index[i][nodeA]==-1||local_index[i][nodeB]==-1 )continue;
      //skip Electrode boundary and interface
      if(edge_array[segment].bc_type == Electrode) continue;
      if(edge_array[segment].bc_type == Interface) continue;
      *psegment++ = local_index[i][nodeA]+1;
      *psegment++ = local_index[i][nodeB]+1;
      csegment++;
    }
    if(csegment>0)
    {
      sprintf(bcname,"%s_Neumann",region_array[i].name);
      segrange[1] = segrange[0] + csegment - 1;
      cg_section_write(fn,B,Z,bcname,BAR_2,
                       segrange[0],segrange[1],0,asegment,&S);
      cg_boco_write(fn,B,Z,bcname,BCTypeUserDefined,ElementRange,2,segrange,&BC);
      segrange[0] = segrange[1]+1;
    }
    delete [] asegment;

    // write doping data to cgns in the solution data node
    // the region type must be Si,GaAs,Ge or Semi at present
    // other type not supported yet

    if(IsSemiconductor(region_array[i].type))
    {
      double *Nd=  new double[node_array.size()];
      double *Na=  new double[node_array.size()];
      int Net=-1,Total=-1;
      for(size_t j=0;j<solhead.sol_name_array.size();j++)
      {
        if(!strcmp(solhead.sol_name_array[j].c_str(),"Net"))  Net=j;
        if(!strcmp(solhead.sol_name_array[j].c_str(),"Total"))  Total=j;
      }
      if(Net != -1 && Total != -1)
      {
        for(size_t j = 0; j < soldata.size(); j++)
        {
          int node_index = soldata[j].index;
          if(soldata[j].material!=region_array[i].type) continue;
          if(local_index[i][node_index]>-1)
          {
            Nd[local_index[i][node_index]] = fabs(soldata[j].data_array[Total] + soldata[j].data_array[Net])/2;
            Na[local_index[i][node_index]] = fabs(soldata[j].data_array[Total] - soldata[j].data_array[Net])/2;
          }
        }
        cg_sol_write(fn,B,Z,"Doping",Vertex,&SOL);
        cg_field_write(fn,B,Z,SOL,RealDouble,"Nd",Nd, &F);
        cg_field_write(fn,B,Z,SOL,RealDouble, "Na",Na, &F);
      }
      delete [] Nd;
      delete [] Na;
    }

    // write mole data to cgns in the solution data node
    // medici only support Single Compound Semiconductor
    if(IsSingleCompSemiconductor(region_array[i].type))
    {
      if(component_array.size())
      {
        double *xmole=  new double[region_array[i].node_num];
        for(size_t j = 0; j < region_array[i].node_num; j++)
          for(size_t k = 0; k < component_array.size(); k++)
          {
            if(component_array[k].region == i)
              component_array[k].mole(x[j],y[j],xmole[j]);
          }
        cg_sol_write(fn,B,Z,"Mole",Vertex,&SOL);
        cg_field_write(fn,B,Z,SOL,RealDouble,"mole_x",xmole, &F);
        delete [] xmole;
      }

      else
      {
        for(size_t j=0;j<parameter_array.size();j++)
          if(parameter_array[j].region-1 == i)
          {
            for(size_t k=0;k<parameter_array[j].parameter_name_array.size();k++)
              if(parameter_array[j].parameter_name_array[k] == "XMOLE")
              {
                double *xmole=  new double[ region_array[i].node_num];
                for(int h = 0; h <  region_array[i].node_num; h++)
                  xmole[h] = parameter_array[j].parameter_value_array[k];
                cg_sol_write(fn,B,Z,"Mole",Vertex,&SOL);
                cg_field_write(fn,B,Z,SOL,RealDouble,"mole_x",xmole, &F);
                delete [] xmole;
                break;
              }
          }
      }
    }

    delete [] x;
    delete [] y;
    delete [] elem;

  }

  /* all things are done */
  for(size_t i =0; i<region_num; i ++)
  {
    delete [] local_index[i];
    delete [] conn[i];
  }
  delete  [] local_index;
  delete  [] conn;
  delete  [] conn_nodes;
  /* close CGNS file */
  cg_close(fn);
}
