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

#include <list>
#include <stack>
#include <queue>
#include <algorithm>
#include <math.h>
#include "zone.h"
#include "geom.h"
#include "log.h"
using namespace std;

/* ----------------------------------------------------------------------------
 * ZONE::build_tricell_data:  This function compute circle center, and the distance between
 * the circle center and each edge. if circle center locates out of the triangle, negtive the
 * corresponding distance.
 */
int ZONE::build_tricell_data()
{
  double PI=3.14159265358979323846;
  for(int i=0;i<datri.size();i++)
  {
    double x1,y1,x2,y2,x3,y3;
    double a,b,c;
    x1 = danode[datri[i].node[0]].x;
    y1 = danode[datri[i].node[0]].y;
    x2 = danode[datri[i].node[1]].x;
    y2 = danode[datri[i].node[1]].y;
    x3 = danode[datri[i].node[2]].x;
    y3 = danode[datri[i].node[2]].y;
    datri[i].edge_len[0] = a = Distance(x3,y3,x2,y2);
    datri[i].edge_len[1] = b = Distance(x1,y1,x3,y3);
    datri[i].edge_len[2] = c = Distance(x1,y1,x2,y2);
    datri[i].angle[0] = acos((b*b+c*c-a*a)/(2*b*c));
    datri[i].angle[1] = acos((c*c+a*a-b*b)/(2*c*a));
    datri[i].angle[2] = acos((a*a+b*b-c*c)/(2*a*b));
    datri[i].area = TriArea(x1,y1,x2,y2,x3,y3);
    if(CircleCenter(x1,y1,x2,y2,x3,y3,datri[i].xc,datri[i].yc)) 
    {
      sprintf(log_buf,"warning! CircleCenter: the three nodes of Tri %d of zone %d coplanar.\n",i,zone_index);
      GSS_LOG();
    }
    datri[i].d[0] = DistancePointLine(x2,y2,x3,y3,datri[i].xc,datri[i].yc);
    datri[i].d[1] = DistancePointLine(x1,y1,x3,y3,datri[i].xc,datri[i].yc);
    datri[i].d[2] = DistancePointLine(x1,y1,x2,y2,datri[i].xc,datri[i].yc);
    if(CompAngle(x3,y3,x1,y1,x2,y2)>PI/2) datri[i].d[0]*=-1;
    if(CompAngle(x1,y1,x2,y2,x3,y3)>PI/2) datri[i].d[1]*=-1;
    if(CompAngle(x2,y2,x3,y3,x1,y1)>PI/2) datri[i].d[2]*=-1;
    datri[i].s[0] = 0.5*datri[i].edge_len[0]*datri[i].d[0];
    datri[i].s[1] = 0.5*datri[i].edge_len[1]*datri[i].d[1];
    datri[i].s[2] = 0.5*datri[i].edge_len[2]*datri[i].d[2];
    for(int j=0;j<dasegment.size();j++)
    {
      if(dasegment[j].get_edge_bc_index(datri[i].node[0],datri[i].node[1]))
         datri[i].bc[2] = dasegment[j].bc_index;
      if(dasegment[j].get_edge_bc_index(datri[i].node[1],datri[i].node[2]))
         datri[i].bc[0] = dasegment[j].bc_index;
      if(dasegment[j].get_edge_bc_index(datri[i].node[2],datri[i].node[0]))
         datri[i].bc[1] = dasegment[j].bc_index;
    }
  }
  return 0;
}

/* ----------------------------------------------------------------------------
 * ZONE::setup_voronoicell:  This function setup voronoicell data structure
 * the adjacency infomation will be build here.
 * if the Triangle is count clock wise ordered, this function works well
 */

typedef struct
{
  int node1,node2;
  int flag1,flag2;
}
ATABLE;


int ZONE::setup_voronoicell()
{
  int adjtabsize = 0;

  vector<ATABLE> **zone_atable;
  zone_atable = new vector<ATABLE>*[danode.size()];
  for(int i=0;i<danode.size();i++)
    zone_atable[i] = new vector<ATABLE>;

  ATABLE tmp_atable;
  tmp_atable.flag1=0; tmp_atable.flag2=0;

  davcell.clear();
  davcell.Init(danode.size());   //init davcell, cell num = node num

  //insert triangle into each node's adjacent table
  for(int i=0;i<datri.size();i++)
  {

    //triangle
    tmp_atable.node1 =  datri[i].node[1];
    tmp_atable.node2 =  datri[i].node[2];
    zone_atable[datri[i].node[0]]->push_back(tmp_atable);

    tmp_atable.node1 =  datri[i].node[2];
    tmp_atable.node2 =  datri[i].node[0];
    zone_atable[datri[i].node[1]]->push_back(tmp_atable);

    tmp_atable.node1 =  datri[i].node[0];
    tmp_atable.node2 =  datri[i].node[1];
    zone_atable[datri[i].node[2]]->push_back(tmp_atable);

  }


  stack<int> clockw;
  queue<int> cclockw;
  for(int i=0;i<davcell.size();i++)
  {
    vector<ATABLE>::iterator pt =  zone_atable[i]->begin();
    int v1,v2;

    v1 = pt->node1;
    v2 = pt->node2;
    //search it's neighbor triangles. begin from zone_atable[0]
    //clock wise search.push the nodes into a stack
    clockw.push(v1); pt->flag1=1;
    for(int num=0;num<zone_atable[i]->size();num++)
      for(pt=zone_atable[i]->begin();pt!=zone_atable[i]->end();pt++)
        if(pt->node2==v1&&pt->flag1==0&&pt->flag2==0)
        {

          v1=pt->node1;         pt->flag1 =1;   pt->flag2 =1;
          if(v1!=v2) clockw.push(v1);

        }
    //count clock wise search, push into a queue
    pt =  zone_atable[i]->begin();
    pt->flag2 = 1;
    cclockw.push(v2);

    for(int num=0;num<zone_atable[i]->size();num++)
      for(pt=zone_atable[i]->begin();pt!=zone_atable[i]->end();pt++)
        if(pt->node1==v2&&pt->flag1==0&&pt->flag2==0)
        {
          v2=pt->node2;         pt->flag1 =1;   pt->flag2 =1;
          cclockw.push(v2);
        }
    //alloc memory for each vcell, set pointers,cp bc_index form danode
    VoronoiCell *pCell = davcell.GetPointer(i);
    int clockwsize = clockw.size();
    int cclockwsize = cclockw.size();
    int neighbor = clockw.size() + cclockw.size();

    pCell->nb_num = neighbor;
    pCell->bc_index = danode[i].bc_index;
    pCell->nb_array = new int[neighbor];
    pCell->elen  = new double[neighbor];
    pCell->ilen  = new double[neighbor];
    pCell->angle = new double[neighbor];
    pCell->celledge = new int[neighbor];
    pCell->inb_array = new int[neighbor];

    for(int j=0;j<pCell->nb_num;j++)  pCell->celledge[j] = -1;

    //now read the neighbors of this node in count clock order
    for(int j=0;j<clockwsize;j++)
    {
      pCell->nb_array[j] = clockw.top();
      clockw.pop();
    }
    for(int j=clockwsize;j<cclockwsize+clockwsize;j++)
    {
      pCell->nb_array[j] = cclockw.front();
      cclockw.pop();
    }

    adjtabsize += neighbor;

  }

  //adjtabsize equal to vedge*2
  davedge.clear();
  davedge.resize(adjtabsize/2);


  //compute invert index of neighbour nodes
  //i == davcell[j].nb_array[davcell[i].inb_array[j]]
  for(int i=0;i<davcell.size();i++)
  {
    VoronoiCell *pCell = davcell.GetPointer(i);
    for(int j=0;j<davcell[i].nb_num;j++)
      for(int k=0;k<davcell[davcell[i].nb_array[j]].nb_num;k++)
      {
        if(i==davcell[davcell[i].nb_array[j]].nb_array[k])
          pCell->inb_array[j]=k;
      }
  }

  //not needed any more, free them.
  for(int i=0;i<davcell.size();i++)
    delete zone_atable[i];
  delete [] zone_atable ;
  return 0;

}

/* ----------------------------------------------------------------------------
 * ZONE::setup_voronoicell:  This function setup voronoicell data structure
 * the adjacency infomation will be build here.
 * if the Triangle is not count clock wise ordered, but the zone is rectangular.
 * please use this  function.
 */
/*
int ZONE::setup_voronoicell()
{
        int adjtabsize=0;
        list<int>::iterator pt;
        list<int>**pneighborlist;
        davcell.clear();
        davcell.Init(danode.size());   //init davcell, cell num = node num

        //compute max triangle, it may be useful
        double max_angle = 0;
        for(int i=0;i<datri.size();i++)
        {
                 double angle= MaxTriAngle(danode[datri[i].A].x,danode[datri[i].A].y,
                                           danode[datri[i].B].x,danode[datri[i].B].y,
                                           danode[datri[i].C].x,danode[datri[i].C].y);
                 if(max_angle<angle) max_angle=angle;
                 if(datri[i].D!=-1)
                        angle= MaxTriAngle(danode[datri[i].A].x,danode[datri[i].A].y,
                                           danode[datri[i].C].x,danode[datri[i].C].y,
                                           danode[datri[i].D].x,danode[datri[i].D].y);
                 if(max_angle<angle) max_angle=angle;
        }
        //alloc list. each voronoi cell will have it's own neighborlist.
        //my god ,tens of thounds of lists i should alloc.
        //stupid code

        pneighborlist=new list<int>* [davcell.size()];
        for(int i=0;i<davcell.size();i++)       pneighborlist[i] = new list<int>;
        //first for each cell, insert all neighbor nodes into the each node's neighborlist.
        for(int i=0;i<datri.size();i++)
        {

                //triangle
                pneighborlist[datri[i].A]->push_front(datri[i].C);
                pneighborlist[datri[i].A]->push_front(datri[i].B);
                pneighborlist[datri[i].B]->push_front(datri[i].A);
                pneighborlist[datri[i].B]->push_front(datri[i].C);
                pneighborlist[datri[i].C]->push_front(datri[i].B);
                pneighborlist[datri[i].C]->push_front(datri[i].A);

        }

        //then remove the overlap nodes in pneighborlist
        for(int i=0;i<davcell.size();i++)
          for(pt=pneighborlist[i]->begin();pt!=pneighborlist[i]->end();)
            if(*pt==i||count(pneighborlist[i]->begin(),pneighborlist[i]->end(),*pt)>1 )
                pt = pneighborlist[i]->erase(pt);
            else
                pt++;

        //alloc memory for each vcell, set pointers,cp bc_index form danode
        VoronoiCell *pCell = davcell.GetPointer(0);
        for(int i=0;i<davcell.size();i++,pCell++)
        {
                pCell->nb_num = pneighborlist[i]->size();
                pCell->nb_array = new int[pneighborlist[i]->size()];
                pCell->elen  = new double[pneighborlist[i]->size()];
                pCell->ilen  = new double[pneighborlist[i]->size()];
                pCell->angle = new double[pneighborlist[i]->size()];
                pCell->celledge = new int[pneighborlist[i]->size()];
                pCell->inb_array = new int[pneighborlist[i]->size()];

                adjtabsize += pneighborlist[i]->size();
        }
        //adjtabsize equal to vedge*2
        davedge.clear();
        davedge.resize(adjtabsize/2);

        // use angle to horizontal line to reorder the neighbour nodes as count clocked
        double x1,y1,x2,y2,x3,y3,r,angle;
        map<double,int>  adjmap;
        map<double,int> ::iterator padjmap;
        list<int> NODElist;
        list<int>::iterator pNODElist, pstart;

        for(int i=0;i<davcell.size();i++)
        {       //compute each angle to horizontal line
                //use container map to order them
                for(pt=pneighborlist[i]->begin();pt!=pneighborlist[i]->end();pt++)
                {
                        x2 = danode[i].x;               y2 = danode[i].y;
                        x1 = danode[*pt].x;             y1 = danode[*pt].y;
                        r  = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
                        x3 = danode[i].x+r;             y3 = danode[i].y;
                        angle = CompAngle(x1,y1,x2,y2,x3,y3);
                        adjmap.insert(map<double,int>::value_type(angle,*pt));
                }
                int    *pnb_index = davcell.GetPointer(i)->nb_array;
                double *pangle = davcell.GetPointer(i)->angle;

                for(padjmap = adjmap.begin();padjmap!=adjmap.end();padjmap++)
                {
                        *pnb_index++ = padjmap->second;
                        *pangle++ = padjmap->first;
                }
                adjmap.clear();

                //boundary cell may have problem. the nb_array may be dissevered by bd nodes
                //i have to reorder it to make sure that bd nodes are located at first and last
                //positions of nb_array array
                int do_reorder = 0;
                if( davcell[i].bc_index != 0) //it is a bd cell
                {       for(int j=0;j<davcell[i].nb_num-1;j++)
                           if(fabs(davcell[i].angle[j]-davcell[i].angle[j+1])>max_angle)
                                 {do_reorder=j+1;break;}

                        // the angle between line i,j and line i,j+1 are larger than max_angle
                        // denote that j and j+1 are bd nodes.
                        // j+1 will be set at the beginning of nb_array ,j set as last value
                }
                if(do_reorder)
                {

                        //push all the neighbors into a list
                        for(int j=0;j<davcell[i].nb_num;j++)
                                NODElist.push_back(davcell[i].nb_array[j]);
                        pNODElist=NODElist.begin();
                        while((*pNODElist)!=davcell[i].nb_array[do_reorder]) pNODElist++;
                        //find do_reorder point
                        pstart = pNODElist;
                        pnb_index = davcell.GetPointer(i)->nb_array;
                        //change the order
                        for(;pNODElist!=NODElist.end();pNODElist++)
                                *pnb_index++ = *pNODElist;
                        for(pNODElist=NODElist.begin();pNODElist!=pstart;pNODElist++)
                                *pnb_index++ = *pNODElist;
                        //recoder finished
                        NODElist.clear();
                        do_reorder=0;
                        //recompute angle
                        pangle = davcell.GetPointer(i)->angle;
                        for(int j=0;j<davcell[i].nb_num;j++)
                        {
                                x2 = danode[i].x;               y2 = danode[i].y;
                                Node nb_node = danode[davcell[i].nb_array[j]];
                                x1 = nb_node.x;
                                y1 = nb_node.y;
                                r  = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
                                x3 = x2+r;
                                y3 = y2;
                                *pangle++ = CompAngle(x1,y1,x2,y2,x3,y3);
                        }
                 }


        }
        //compute invert index of neighbour nodes
        //davcell[i] == davcell[davcell[j].nb_array[davcell[i].inb_array[j]]]
        for(int i=0;i<davcell.size();i++)
        {
                VoronoiCell *pCell = davcell.GetPointer(i);
                for(int j=0;j<davcell[i].nb_num;j++)
                  for(int k=0;k<davcell[davcell[i].nb_array[j]].nb_num;k++)
                  {
                        if(i==davcell[davcell[i].nb_array[j]].nb_array[k])
                                pCell->inb_array[j]=k;
                  }
        }

        //not needed any more, free them.
        for(int i=0;i<davcell.size();i++)
                delete pneighborlist[i];
        delete [] pneighborlist ;
        return 0;

}
*/


/* ----------------------------------------------------------------------------
 * ZONE::setup_voronoicell: the remain work is eary to do. ilen,elen and area
 * will be setted here.
 */
int ZONE::build_voronoicell_data()
{
  //copy some data from danode
  for(int i=0;i<davcell.size();i++)
  {
    VoronoiCell *pCell = davcell.GetPointer(i);
    pCell->x = danode[i].x;
    pCell->y = danode[i].y;
    pCell->bc_index = danode[i].bc_index;
  }

  //comput angle: the angle from horizontal line to it's neighbor node
  double x1,y1,x2,y2,x3,y3,r,angle;
  for(int i=0;i<davcell.size();i++)
  {
    double *pangle = davcell.GetPointer(i)->angle;
    for(int j=0;j<davcell[i].nb_num;j++)
    {
      x2 = danode[i].x;         y2 = danode[i].y;
      Node nb_node = danode[davcell[i].nb_array[j]];
      x1 = nb_node.x;
      y1 = nb_node.y;
      r  = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
      x3 = x2+r;
      y3 = y2;
      *pangle++ = CompAngle(x1,y1,x2,y2,x3,y3);
    }
  }

  //comput ilen: length between node i and it's neighbor
  for(int i=0;i<davcell.size();i++)
  {
    VoronoiCell *pCell = davcell.GetPointer(i);
    for(int j=0;j<davcell[i].nb_num;j++)
    {
      Node n = danode[davcell[i].nb_array[j]];
      pCell->ilen[j]=Distance(danode[i].x,danode[i].y,n.x,n.y);
    }
  }

  //next to comput elen , length of Voronoi Cell's Edge
  //at the same time ,davedge will be set
  list<double> listxc,listyc; //two list to store location of CircleCenter
  list<double>::iterator pxc,pyc;
  double xc,yc,xc_next,yc_next;
  vector<VEdge>:: iterator pEdge = davedge.begin();
  int nedge = 0;
  for(int i=0;i<davcell.size();i++)
  {
    //fill the list with CircleCenter
    for(int j=0;j<davcell[i].nb_num-1;j++)
    {
      int error=CircleCenter(danode[i].x,danode[i].y,
                   danode[davcell[i].nb_array[j]].x,
                   danode[davcell[i].nb_array[j]].y,
                   danode[davcell[i].nb_array[j+1]].x,
                   danode[davcell[i].nb_array[j+1]].y,
                   xc,yc);
      if(error) 
      {
        sprintf(log_buf,"warning! CircleCenter: the three nodes of vcell %d of zone %d coplanar.\n",i,zone_index);
        GSS_LOG();
      }
      listxc.push_back(xc);
      listyc.push_back(yc);
    }
    int nb_first= davcell[i].nb_array[0];
    int nb_last = davcell[i].nb_array[davcell[i].nb_num-1];
    //if all the 3 points are bd ndoes
    if (davcell[i].bc_index&&davcell[nb_first].bc_index&&
        davcell[nb_last].bc_index)
    {
      xc = (danode[i].x+danode[nb_first].x)/2.0;
      yc = (danode[i].y+danode[nb_first].y)/2.0;
      listxc.push_front(xc);
      listyc.push_front(yc);
      xc = (danode[i].x+danode[nb_last].x)/2.0;
      yc = (danode[i].y+danode[nb_last].y)/2.0;
      listxc.push_back(xc);
      listyc.push_back(yc);
    }
    else //compute CircleCenter of the last triangle
    {
      int error=CircleCenter(danode[i].x,danode[i].y,
                   danode[nb_first].x,
                   danode[nb_first].y,
                   danode[nb_last].x,
                   danode[nb_last].y,
                   xc,yc);
      if(error) 
      {
        sprintf(log_buf,"warning! CircleCenter: the three nodes of vcell %d of zone %d coplanar.\n",i,zone_index);
        GSS_LOG();
      }
      listxc.push_front(xc);
      listyc.push_front(yc);
      listxc.push_back(xc);
      listyc.push_back(yc);
    }

    pxc=listxc.begin();
    pyc=listyc.begin();
    VoronoiCell *pCell = davcell.GetPointer(i);

    for(int j=0;j<davcell[i].nb_num;j++)
    {
      xc = *pxc++;       yc = *pyc++;
      xc_next = *pxc;    yc_next = *pyc;
      pCell->elen[j] = Distance(xc,yc,xc_next,yc_next);

      // each vedge belongs to 2 vcells

      int i_index = davcell[i].inb_array[j];
      if(davcell[i].celledge[j]==-1 && davcell[davcell[i].nb_array[j]].celledge[i_index]==-1)
      {       //this edge not registered. set data to it
        pEdge->cell1 = i;
        pEdge->cell2 = davcell[i].nb_array[j];
        pEdge->c1_nb = j;
        pEdge->c2_nb = davcell[i].inb_array[j];

        pEdge->length = pCell->elen[j];
        pEdge->x1 = xc;       pEdge->y1 = yc;
        pEdge->x2 = xc_next;  pEdge->y2 = yc_next;
        davcell.GetPointer(i)->celledge[j] = nedge;
        davcell.GetPointer(davcell[i].nb_array[j])->celledge[i_index] = nedge;
        pEdge++; nedge++;
      }
    }

    listxc.clear();
    listyc.clear();
  }

  // compute each cell's volume and zone area
  field_area_vcell = 0;
  for(int i=0;i<davcell.size();i++)
  {
    VoronoiCell *pCell = davcell.GetPointer(i);
    pCell->area = 0;
    for(int j=0;j<davcell[i].nb_num;j++)
      pCell->area += davcell[i].elen[j]*davcell[i].ilen[j]/4.0;
    field_area_vcell += pCell->area;
  }
  return 0;

}
