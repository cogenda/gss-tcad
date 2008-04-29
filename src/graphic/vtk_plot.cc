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
/*  Last update: Oct 30, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#define HAVE_VTK
#ifdef HAVE_VTK

#include <unistd.h>
#include "bsolver.h"
#include "log.h"

int BSolver::vtk_display_data(PlotDefine & plot_define)
{
    // This is the child process.
    char filename[32];
    char *vtk_wish = new char [strlen(getenv("GSS_DIR"))+32];
    char *vtk_script = new char [strlen(getenv("GSS_DIR"))+32];
    sprintf(vtk_wish,"%s/bin/vtk",getenv("GSS_DIR"));
    sprintf(vtk_script,"%s/bin/vtk_plot.tcl",getenv("GSS_DIR"));

    if ( plot_define.Variable == Potential)
    {
      sprintf(filename,"data_%s.vtk","Potential");
      vtk_output_file(filename);
      if(execl(vtk_wish, "vtk",vtk_script,filename,"Potential",0)==-1)
      { gss_log.string_buf()<<"VTK plot error!";   gss_log.record();}
    }
    else if (plot_define.Variable == DopingNa)
    {
      sprintf(filename,"data_%s.vtk","Na");
      vtk_output_file(filename);
      if(execl(vtk_wish, "vtk",vtk_script,filename,"Na",0)==-1)
      { gss_log.string_buf()<<"VTK plot error!";   gss_log.record();}
    }
    else if (plot_define.Variable == DopingNd)
    {
      sprintf(filename,"data_%s.vtk","Nd");
      vtk_output_file(filename);
      if(execl(vtk_wish, "vtk",vtk_script,filename,"Nd",0)==-1)
      { gss_log.string_buf()<<"VTK plot error!";   gss_log.record();}
    }
    else if (plot_define.Variable == ElecDensity)
    {
      sprintf(filename,"data_%s.vtk","ElecDensity");
      vtk_output_file(filename);
      if(execl(vtk_wish, "vtk",vtk_script,filename,"ElecDensity",0)==-1)
      { gss_log.string_buf()<<"VTK plot error!";   gss_log.record();}
    }
    else if (plot_define.Variable == HoleDensity)
    {
      sprintf(filename,"data_%s.vtk","HoleDensity");
      vtk_output_file(filename);
      if(execl(vtk_wish, "vtk",vtk_script,filename,"HoleDensity",0)==-1)
      { gss_log.string_buf()<<"VTK plot error!";   gss_log.record();}
    }
    else if (plot_define.Variable == Temperature)
    {
      sprintf(filename,"data_%s.vtk","Temperature");
      vtk_output_file(filename);
      if(execl(vtk_wish, "vtk",vtk_script,filename,"Temperature",0)==-1)
      { gss_log.string_buf()<<"VTK plot error!";   gss_log.record();}
    }
    else if (plot_define.Variable == OpticalEz)
    {
      sprintf(filename,"data_%s.vtk","OpticalEz");
      vtk_output_file(filename);
      if(execl(vtk_wish, "vtk",vtk_script,filename,"OpticalEz",0)==-1)
      { gss_log.string_buf()<<"VTK plot error!";   gss_log.record();}
    }
    delete [] vtk_wish;
    delete [] vtk_script;
    return 0;
}



int BSolver::vtk_output_file(char *filename)
{
  FILE *fp;
  fp = fopen(filename,"w");
  if(fp==NULL)
  {
    sprintf(log_buf,"Error, can't open file for vtk data output\n");
    GSS_LOG();
    return 1;
  }

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"GSS data file\n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"DATASET UNSTRUCTURED_GRID\n\n");
  fprintf(fp,"POINTS %d double\n",gnode.size());
  for(int i=0;i<gnode.size();i++)
    fprintf(fp,"%e %e %lf\n",gnode[i].x/scale_unit.s_micron,gnode[i].y/scale_unit.s_micron,0.0);
  fprintf(fp,"\nCELLS %d %d\n",gtri.size(),4*gtri.size());
  for(int i=0;i<gtri.size();i++)
    fprintf(fp,"%d %d %d %d\n",3,gtri[i].g_node[0],gtri[i].g_node[1],gtri[i].g_node[2]);
  fprintf(fp,"CELL_TYPES %d\n",gtri.size());
  for(int i=0;i<gtri.size();i++)
    fprintf(fp,"%d\n",5);
  fprintf(fp,"\n");

  fprintf(fp,"\nPOINT_DATA %d\n",gnode.size());
  
  fprintf(fp,"SCALARS Potential double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].P/scale_unit.s_volt));
    }
    else if(zonedata[zone_index]->material_type == Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].P/scale_unit.s_volt));
    }
    else if(zonedata[zone_index]->material_type == Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].P/scale_unit.s_volt));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }

  fprintf(fp,"SCALARS ElecDensity double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].n*pow(scale_unit.s_centimeter,3)));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }

  fprintf(fp,"SCALARS HoleDensity double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].p*pow(scale_unit.s_centimeter,3)));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }


  fprintf(fp,"SCALARS Nd double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->aux[local_index].Nd*pow(scale_unit.s_centimeter,3)));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }

  fprintf(fp,"SCALARS Na double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->aux[local_index].Na*pow(scale_unit.s_centimeter,3)));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }


  fprintf(fp,"SCALARS Temperature double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].T/scale_unit.s_kelvin));
    }
    else if(zonedata[zone_index]->material_type == Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].T/scale_unit.s_kelvin));
    }
    else if(zonedata[zone_index]->material_type == Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->fs[local_index].T/scale_unit.s_kelvin));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }
  
  fprintf(fp,"SCALARS opticalEz double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpEz/scale_unit.s_volt*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == Insulator)
    {
     ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[zone_index]);
     fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpEz/scale_unit.s_volt*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == Conductor)
    {
     ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[zone_index]);
     fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpEz/scale_unit.s_volt*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == Vacuum)
    {
     VacuumZone *pzonedata = dynamic_cast< VacuumZone * >(zonedata[zone_index]);
     fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpEz/scale_unit.s_volt*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == PML)
    {
      PMLZone *pzonedata = dynamic_cast< PMLZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpEz/scale_unit.s_volt*scale_unit.s_centimeter));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }
  
  fprintf(fp,"SCALARS OptG double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OptG*pow(scale_unit.s_centimeter,3)*scale_unit.s_second));
    }
    else
      fprintf(fp,"%e\n",0.0);
  }
  
  fprintf(fp,"SCALARS opticalHz double 1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpHz/scale_unit.s_A*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == Insulator)
    {
     ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[zone_index]);
     fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpHz/scale_unit.s_A*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == Conductor)
    {
     ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[zone_index]);
     fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpHz/scale_unit.s_A*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == Vacuum)
    {
     VacuumZone *pzonedata = dynamic_cast< VacuumZone * >(zonedata[zone_index]);
     fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpHz/scale_unit.s_A*scale_unit.s_centimeter));
    }
    else if(zonedata[zone_index]->material_type == PML)
    {
      PMLZone *pzonedata = dynamic_cast< PMLZone * >(zonedata[zone_index]);
      fprintf(fp,"%e\n",double(pzonedata->aux[local_index].OpHz/scale_unit.s_A*scale_unit.s_centimeter));
    }
  }
  
  fprintf(fp,"VECTORS EField double \n");
  for(int i=0;i<gnode.size();i++)
  {
    int zone_index = gnode[i].zone_index;
    int local_index = gnode[i].local_index;
    if(zonedata[zone_index]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[zone_index]);
      fprintf(fp,"%e %e 0\n",double(pzonedata->aux[local_index].Ex), double(pzonedata->aux[local_index].Ey));
    }
    else if(zonedata[zone_index]->material_type == Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[zone_index]);
      fprintf(fp,"%e %e 0\n",double(pzonedata->aux[local_index].Ex), double(pzonedata->aux[local_index].Ey));
    }
    else
      fprintf(fp,"0 0 0\n");
  }
  
  
  
  fclose(fp);
  return 0;
}

#endif
