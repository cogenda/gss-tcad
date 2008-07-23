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
/*  Last update: Dec 19, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cgnslib.h>
#include <math.h>
#include <iostream>
#include "bsolver.h"
#include "geom.h"
#include "typedef.h"
#define ANSI_DECLARATORS

extern "C"
{
#include "triangle.h"
}

inline double FunLinear(PetscScalar value)
{
  return double(value);
}

inline double FunSignedLog(PetscScalar value)
{
  return double(dsign(value)*log(1.0+fabs(value)));
}

/* ----------------------------------------------------------------------------
 * BSolver::set_scale  This function compute the refine scale
 */
double  BSolver::set_refine_scale(Tri & tri, RefineDefine & ref)
{
  double M1=0,M2=0,M3=0,M,Scale=1.1;
  double (*Measure) (PetscScalar);
  int z = tri.zone_index;

  if(zonedata[z]->material_type == Semiconductor)
  {
    if(ref.Measure==Linear)
      Measure = FunLinear;
    if(ref.Measure==SignedLog)
      Measure = FunSignedLog;
    SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);

    int p1 = tri.node[0];
    int p2 = tri.node[1];
    int p3 = tri.node[2];

    if(ref.Variable == Potential)
    {
      M1 = Measure(pzonedata->fs[p1].P);
      M2 = Measure(pzonedata->fs[p2].P);
      M3 = Measure(pzonedata->fs[p3].P);
      M =  dmax(dmax(fabs(M1-M2),fabs(M2-M3)),fabs(M3-M1))+1e-6;
      if(M>ref.Dispersion)
        Scale = ref.Dispersion/M;
      if(Scale<ref.DivisionRatio)        Scale = ref.DivisionRatio;
      return Scale;
    }
    if(ref.Variable == Doping)
    {
      // the mesh size must small than Debye length.
      PetscScalar doping = fabs(pzonedata->aux[p1].Net_doping()
                               +pzonedata->aux[p2].Net_doping()
	                       +pzonedata->aux[p3].Net_doping())/3.0
			       +1.0*pow(scale_unit.s_centimeter,-3);
      PetscScalar  Ld = sqrt(pzonedata->aux[p1].eps*pzonedata->mt->kb*pzonedata->fs[p1].T/(pzonedata->mt->e*pzonedata->mt->e*doping));

      M1 = Measure(pzonedata->aux[p1].Net_doping()*pow(scale_unit.s_centimeter,3));
      M2 = Measure(pzonedata->aux[p2].Net_doping()*pow(scale_unit.s_centimeter,3));
      M3 = Measure(pzonedata->aux[p3].Net_doping()*pow(scale_unit.s_centimeter,3));
      M =  dmax(dmax(fabs(M1-M2),fabs(M2-M3)),fabs(M3-M1))+1e-6;
      if(M>ref.Dispersion)
        Scale = Ld*Ld/(2*tri.area) < ref.Dispersion/M ? Ld*Ld/(2*tri.area):ref.Dispersion/M;
      //if(Scale>Ld/sqrt(2*tri.area))   Scale = Ld*Ld/(2*tri.area);
      if(Scale<ref.DivisionRatio)       Scale = ref.DivisionRatio;
      return Scale;
    }

  }
  else if(zonedata[z]->material_type == Insulator)
  {
    return      1.1;
  }
  else if(zonedata[z]->material_type == Conductor)
  {
    return      1.1;
  }
  else if(zonedata[z]->material_type == Vacuum)
  {
    return      1.1;
  }
  else if(zonedata[z]->material_type == PML)
  {
    return      1.1;
  }
  return 0.0; //prevent warning of complier
}

/* ----------------------------------------------------------------------------
 * BSolver::refine  This function call triangulate to refine the mesh
 *
 */

int BSolver::refine(RefineDefine & ref)
{
  struct triangulateio in, out;

  //do necessarily prepare for call triangulate
  in.numberofpoints = gnode.size();
  in.numberofpointattributes = 10; //P, n, p, Nd, Na, P, As, Sb, B and mole_x

  in.pointattributelist = (double *)malloc(in.numberofpoints * in.numberofpointattributes * sizeof(double));
  in.pointlist = (double *) malloc(in.numberofpoints * 2 * sizeof(double));
  in.pointmarkerlist = (int *) malloc(in.numberofpoints * sizeof(int));
  for(int i=0;i<in.numberofpoints;i++)
  {
    in.pointlist[2*i] = gnode[i].x;
    in.pointlist[2*i+1] = gnode[i].y;
    in.pointmarkerlist[i] = gnode[i].bc_index;
  }
  double *ppointattributelist = in.pointattributelist;
/*
  // set all the interface node in global node array as semiconductor node
  for(int z=0;z<zone_num;z++)
    if(zonedata[z]->material_type == Semiconductor)
    {
      for(int i=0;i<zone[z].danode.size();i++)
      {
        gnode[zone[z].danode[i].g_index].local_index = i;
        gnode[zone[z].danode[i].g_index].zone_index = z;
      }
    }
*/
  for(int i=0;i<in.numberofpoints;i++)
  {
    int z = gnode[i].zone_index;
    if(zonedata[z]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      int n = gnode[i].local_index;
      //*ppointattributelist++ = log(1.0+pzonedata->fs[n].Nd*pow(scale_unit.s_centimeter,3));
      //*ppointattributelist++ = log(1.0+pzonedata->fs[n].Na*pow(scale_unit.s_centimeter,3));
      *ppointattributelist++ = double(pzonedata->fs[n].P);
      *ppointattributelist++ = double(pzonedata->fs[n].n*pow(scale_unit.s_centimeter,3));
      *ppointattributelist++ = double(pzonedata->fs[n].p*pow(scale_unit.s_centimeter,3));
      
      *ppointattributelist++ = double(pzonedata->aux[n].Nd*pow(scale_unit.s_centimeter,3));
      *ppointattributelist++ = double(pzonedata->aux[n].Na*pow(scale_unit.s_centimeter,3));
      *ppointattributelist++ = double(pzonedata->aux[n].P*pow(scale_unit.s_centimeter,3));
      *ppointattributelist++ = double(pzonedata->aux[n].As*pow(scale_unit.s_centimeter,3));
      *ppointattributelist++ = double(pzonedata->aux[n].Sb*pow(scale_unit.s_centimeter,3));
      *ppointattributelist++ = double(pzonedata->aux[n].B*pow(scale_unit.s_centimeter,3));
      
      *ppointattributelist++ = double(pzonedata->aux[n].mole_x);
    }
    else if(zonedata[z]->material_type == Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      int n = gnode[i].local_index;
      *ppointattributelist++ = double(pzonedata->fs[n].P);
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
    }
    else
    {
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
      *ppointattributelist++ = 0;
    }
  }

  in.numberoftriangles = gtri.size();
  in.numberofcorners = 3;
  in.numberoftriangleattributes = 1;
  in.trianglelist =  (int *) malloc(in.numberoftriangles * 3 * sizeof(int));
  in.trianglearealist = (double *) malloc(in.numberoftriangles * sizeof(double));
  in.triangleattributelist = (double *) malloc(in.numberoftriangles * sizeof(double));
  int *ptrianglelist = in.trianglelist;
  double *ptrianglearealist = in.trianglearealist;
  double *ptriangleattributelist = in.triangleattributelist;
  for(int i=0;i<gtri.size();i++)
  {
     *ptrianglelist++ = gtri[i].g_node[0];
     *ptrianglelist++ = gtri[i].g_node[1];
     *ptrianglelist++ = gtri[i].g_node[2];
     *ptriangleattributelist++ = gtri[i].zone_index;
  }
  for(int i=0;i<gtri.size();i++)
  {
    double Scale=0.5;
    Scale = set_refine_scale(gtri[i],ref);
    *ptrianglearealist++ = gtri[i].area*Scale;
  }
  in.numberofsegments =  0;
  for(int i=0;i<gsegment.size();i++)
    in.numberofsegments += gsegment[i].edge_num;
  in.segmentlist =  (int *) malloc(in.numberofsegments * 2 * sizeof(int));
  in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
  int *psegmentlist =  in.segmentlist;
  int *psegmentmarkerlist = in.segmentmarkerlist;
  for(int i=0;i<gsegment.size();i++)
    for(int j=0;j<gsegment[i].edge_num;j++)
    {
      *psegmentlist++ = gsegment[i].edge_array[j].p1;
      *psegmentlist++ = gsegment[i].edge_array[j].p2;
      *psegmentmarkerlist++ = gsegment[i].bc_index;
    }

  in.numberofholes = 0;
  in.numberofregions = 0;


  // Make necessary initializations so that Triangle can return a
  // triangulation in `out'.

  out.pointlist = (double *) NULL;
  out.pointattributelist = (double *) NULL;
  out.pointmarkerlist = (int *) NULL;
  out.trianglelist = (int *) NULL;
  out.triangleattributelist = (double *) NULL;
  out.segmentlist = (int *) NULL;
  out.segmentmarkerlist = (int *) NULL;

  // Refine the triangulation according to the attached
  // triangle area constraints.

  triangulate(ref.tri_cmd, &in, &out, (struct triangulateio *) NULL);

  gnode.resize(out.numberofpoints);
  gtri.resize(out.numberoftriangles);
  for(int i=0; i<gsegment.size(); i++)
  {
    gsegment[i].edge_array.clear();
    gsegment[i].edge_num = 0;
  }

  ppointattributelist = out.pointattributelist;
  for(int i=0;i<out.numberofpoints;i++)
  {
    gnode[i].x = out.pointlist[2*i];
    gnode[i].y = out.pointlist[2*i+1];
    gnode[i].bc_index = out.pointmarkerlist[i];
    gnode[i].g_index = i;
    gnode[i].zone_index = -1;
    gnode[i].local_index = -1;
  }

  ptrianglelist = out.trianglelist;
  for(int i=0;i<out.numberoftriangles;i++)
  {
    gtri[i].zone_index = static_cast<int>( out.triangleattributelist[i]+0.5 );
    gtri[i].g_node[0] = *ptrianglelist++;
    gtri[i].g_node[1] = *ptrianglelist++;
    gtri[i].g_node[2] = *ptrianglelist++;
    gtri[i].bc[0] = 0;
    gtri[i].bc[1] = 0;
    gtri[i].bc[2] = 0;
  }

  for(int j=0;j<gsegment.size();j++)
    for(int i=0;i<out.numberofsegments;i++)
      if(out.segmentmarkerlist[i] == gsegment[j].bc_index)
      {
        Edge tmp_edge;
        tmp_edge.p1 = out.segmentlist[2*i];
        tmp_edge.p2 = out.segmentlist[2*i+1];
        tmp_edge.bc_index = out.segmentmarkerlist[i];
        gsegment[j].edge_array.push_back(tmp_edge);
        gsegment[j].edge_num++;
      }

  field_to_zone();
  build_zone();


  clear_data();
  bc.clear();
  if(setup_bc())   return 1;
  if(build_zonedata())    return 1;

  //re-arrange doping data
  if(doping_func.size()) //if we have analytic function?
        setup_doping();
  else // we have to use 
  {
    ppointattributelist = out.pointattributelist;
    for(int z=0;z<zone_num;z++)
      if(zonedata[z]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].danode.size();i++)
      {
        int index = zone[z].danode[i].g_index;
        pzonedata->aux[i].Nd = fabs((ppointattributelist[10*index+3]))/pow(scale_unit.s_centimeter,3);
        pzonedata->aux[i].Na = fabs((ppointattributelist[10*index+4]))/pow(scale_unit.s_centimeter,3);
	pzonedata->aux[i].P  = fabs((ppointattributelist[10*index+5]))/pow(scale_unit.s_centimeter,3);
        pzonedata->aux[i].As = fabs((ppointattributelist[10*index+6]))/pow(scale_unit.s_centimeter,3);
	pzonedata->aux[i].Sb = fabs((ppointattributelist[10*index+7]))/pow(scale_unit.s_centimeter,3);
        pzonedata->aux[i].B  = fabs((ppointattributelist[10*index+8]))/pow(scale_unit.s_centimeter,3);
      }
    }
  }

  //mole function
  for(int z=0;z<zone_num;z++)
    if(IsSingleCompSemiconductor(zonedata[z]->pzone->zonelabel))
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].danode.size();i++)
      {
        int index = zone[z].danode[i].g_index;
        pzonedata->aux[i].mole_x = fabs(ppointattributelist[10*index+9]);
        if(pzonedata->aux[i].mole_x > 1.0) pzonedata->aux[i].mole_x = 1.0;
      }
    }
  
  setup_init_data();

  // Free all allocated arrays, including those allocated by Triangle.

  free(in.pointlist);
  free(in.pointmarkerlist);
  free(in.pointattributelist);
  free(in.trianglelist);
  free(in.triangleattributelist);
  free(in.trianglearealist);
  free(in.segmentlist);
  free(in.segmentmarkerlist);

  free(out.pointlist);
  free(out.pointmarkerlist);
  free(out.pointattributelist);
  free(out.trianglelist);
  free(out.triangleattributelist);
  free(out.segmentlist);
  free(out.segmentmarkerlist);

  return 0;
}
