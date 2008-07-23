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
/*  Last update: Nov 29, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#include "bsolver.h"
#include "log.h"

int BSolver::probe_open(int solver_type)
{
  for(int i=0;i<probe_define.size();i++)
  {
    BaseBC * pbc = bc.Get_pointer(probe_define[i].bc_index);
    int node_num=pbc->psegment->node_array.size();
    if(probe_define[i].Append)
      probe_define[i].pFile = fopen(probe_define[i].ProbeFile.c_str(),"a");
    else
      probe_define[i].pFile = fopen(probe_define[i].ProbeFile.c_str(),"w");
    if(!probe_define[i].pFile)
    {sprintf(log_buf,"I can't open file %s for writing!\n",probe_define[i].ProbeFile.c_str());GSS_LOG();return 1;}
    fprintf(probe_define[i].pFile,"#Region: %s   Segment: %s\n",probe_define[i].Region.c_str(), probe_define[i].Segment.c_str());
    fprintf(probe_define[i].pFile,"#%d\tX\tY\t\n",node_num);
    for(int p=0;p<node_num;p++)
    {
      int node = pbc->psegment->node_array[p];
      double x = zone[probe_define[i].zone_index].danode[node].x/scale_unit.s_micron;
      double y = zone[probe_define[i].zone_index].danode[node].y/scale_unit.s_micron;
      fprintf(probe_define[i].pFile,"#%d\t%f\t%f\n",p,x,y);
    }
  }
  return 0;
}



int BSolver::probe(int solver_type, PetscScalar  x)
{
  for(int i=0;i<probe_define.size();i++)
  {
    BaseBC * pbc = bc.Get_pointer(probe_define[i].bc_index);
    int node_num=pbc->psegment->node_array.size();
    //write information
    if(!probe_define[i].flag)
    {
      switch(solver_type)
      {
      case DCSWEEP_VSCAN:
        fprintf(probe_define[i].pFile,"#DCSWEEP_VSCAN(V)\t[%s]\n",probe_define[i].VariableName.c_str());break;
      case DCSWEEP_ISCAN:
        fprintf(probe_define[i].pFile,"#DCSWEEP_ISCAN(mA)\t[%s]\n",probe_define[i].VariableName.c_str());break;
      case TRANSIENT:
        fprintf(probe_define[i].pFile,"#TRANSIENT(ps)\t[%s]\n",probe_define[i].VariableName.c_str());break;
      case EQUILIBRIUM:
        fprintf(probe_define[i].pFile,"#EQUILIBRIUM\t[%s]\n",probe_define[i].VariableName.c_str());break;
      case STEADYSTATE:
        fprintf(probe_define[i].pFile,"#STEADYSTATE\t[%s]\n",probe_define[i].VariableName.c_str());break;
      default:
        fprintf(probe_define[i].pFile,"#[%s]\n",probe_define[i].VariableName.c_str());break;
      }
      probe_define[i].flag=true;
    }
    //write scan information
    switch(solver_type)
    {
    case DCSWEEP_VSCAN:
      fprintf(probe_define[i].pFile,"%e\t",double(x)/scale_unit.s_volt);break;
    case DCSWEEP_ISCAN:
      fprintf(probe_define[i].pFile,"%e\t",double(x)/scale_unit.s_mA);break;
    case TRANSIENT:
      fprintf(probe_define[i].pFile,"%e\t",double(x));break;
    }
    //record data for semiconductor zone only
    if( probe_define[i].Variable==Doping         ||
        probe_define[i].Variable==NetDoping      ||
        probe_define[i].Variable==DopingNd       ||
        probe_define[i].Variable==DopingNa       ||
        probe_define[i].Variable==Phosphorus     ||
        probe_define[i].Variable==Arsenic        ||
        probe_define[i].Variable==Antimony       ||
        probe_define[i].Variable==Boron          ||
	probe_define[i].Variable==ElecDensity    ||
        probe_define[i].Variable==HoleDensity    ||
        probe_define[i].Variable==ElecTemp       ||
        probe_define[i].Variable==HoleTemp       ||
        probe_define[i].Variable==PhiP           ||
        probe_define[i].Variable==PhiN           ||
        probe_define[i].Variable==Phi_Intrinsic  ||
        probe_define[i].Variable==QuantumEc      ||
        probe_define[i].Variable==QuantumEv      ||
        probe_define[i].Variable==OpticalEz      ||
        probe_define[i].Variable==OpticalHz      ||
        probe_define[i].Variable==OpticalG)
    {
      if(zonedata[probe_define[i].zone_index]->material_type == Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[probe_define[i].zone_index]);
        for(int p=0;p<node_num;p++)
        {
          int node = pbc->psegment->node_array[p];
          switch(probe_define[i].Variable)
          {
          case DopingNd:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Total_Nd()*pow(scale_unit.s_centimeter,3)));break;
          case DopingNa:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Total_Na()*pow(scale_unit.s_centimeter,3)));break;
          case Doping:
            fprintf(probe_define[i].pFile,"%e\t",double((pzonedata->aux[node].Total_doping())*pow(scale_unit.s_centimeter,3)));break;
          case NetDoping:
            fprintf(probe_define[i].pFile,"%e\t",double((pzonedata->aux[node].Net_doping())*pow(scale_unit.s_centimeter,3)));break;
	  case Phosphorus:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].P*pow(scale_unit.s_centimeter,3)));break;
          case Arsenic:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].As*pow(scale_unit.s_centimeter,3)));break;
	  case Antimony:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Sb*pow(scale_unit.s_centimeter,3)));break;
          case Boron:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].B*pow(scale_unit.s_centimeter,3)));break;
	  case ElecDensity:
            fprintf(probe_define[i].pFile,"%e\t",double(fabs(pzonedata->fs[node].n)*pow(scale_unit.s_centimeter,3)));break;
          case HoleDensity:
            fprintf(probe_define[i].pFile,"%e\t",double(fabs(pzonedata->fs[node].p)*pow(scale_unit.s_centimeter,3)));break;
          case ElecTemp:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].Tn/scale_unit.s_kelvin));break;
          case HoleTemp:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].Tp/scale_unit.s_kelvin));break;
          case Phi_Intrinsic:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].phi_intrinsic/scale_unit.s_volt));break;
          case PhiP:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].phip/scale_unit.s_volt));break;
          case PhiN:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].phin/scale_unit.s_volt));break;
          case QuantumEc:
            fprintf(probe_define[i].pFile,"%e\t",double((pzonedata->fs[node].Eqc+pzonedata->fs[node].P+pzonedata->aux[node].affinity)/scale_unit.s_volt));break;
          case QuantumEv:
            fprintf(probe_define[i].pFile,"%e\t",double((pzonedata->fs[node].Eqv+pzonedata->fs[node].P+pzonedata->aux[node].affinity+pzonedata->aux[node].Eg)/scale_unit.s_volt));break;
          case OpticalEz:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].OpEz/scale_unit.s_volt*scale_unit.s_centimeter));break;
          case OpticalHz:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].OpHz/scale_unit.s_A*scale_unit.s_centimeter));break;
          case OpticalG:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].OptG*pow(scale_unit.s_centimeter,3)));break;
          }
        }
      }
    }
    //record data for all zones
    if( probe_define[i].Variable==Potential ||
        probe_define[i].Variable==EFieldX   ||
        probe_define[i].Variable==EFieldY   ||
        probe_define[i].Variable==Temperature)
    {
      if(zonedata[probe_define[i].zone_index]->material_type == Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[probe_define[i].zone_index]);
        for(int p=0;p<node_num;p++)
        {
          int node = pbc->psegment->node_array[p];
          switch(probe_define[i].Variable)
          {
          case Potential:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].P/scale_unit.s_volt));break;
          case EFieldX:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Ex/scale_unit.s_volt*scale_unit.s_centimeter));break;
          case EFieldY:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Ey/scale_unit.s_volt*scale_unit.s_centimeter));break;
          case Temperature:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].T/scale_unit.s_kelvin));break;
          }
        }
      }

      if(zonedata[probe_define[i].zone_index]->material_type == Insulator)
      {
        ISZone *pzonedata = dynamic_cast<ISZone *> (zonedata[probe_define[i].zone_index]);
        for(int p=0;p<node_num;p++)
        {
          int node = pbc->psegment->node_array[p];
          switch(probe_define[i].Variable)
          {
          case Potential:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].P/scale_unit.s_volt));break;
          case EFieldX:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Ex/scale_unit.s_volt*scale_unit.s_centimeter));break;
          case EFieldY:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Ey/scale_unit.s_volt*scale_unit.s_centimeter));break;
          case Temperature:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].T/scale_unit.s_kelvin));break;
          }
        }
      }

      if(zonedata[probe_define[i].zone_index]->material_type == Conductor)
      {
        ElZone *pzonedata = dynamic_cast<ElZone *> (zonedata[probe_define[i].zone_index]);
        for(int p=0;p<node_num;p++)
        {
          int node = pbc->psegment->node_array[p];
          switch(probe_define[i].Variable)
          {
          case Potential:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].P/scale_unit.s_volt));break;
          case EFieldX:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Ex/scale_unit.s_volt*scale_unit.s_centimeter));break;
          case EFieldY:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->aux[node].Ey/scale_unit.s_volt*scale_unit.s_centimeter));break;
          case Temperature:
            fprintf(probe_define[i].pFile,"%e\t",double(pzonedata->fs[node].T/scale_unit.s_kelvin));break;
          }
        }
      }

    }
    fprintf(probe_define[i].pFile,"\n");
    fflush(probe_define[i].pFile);
  }

  return 0;
}

int BSolver::probe_close()
{
  for(int i=0;i<probe_define.size();i++)
  {
    if(probe_define[i].pFile)
      fclose(probe_define[i].pFile);
  }
  probe_define.clear();
  return 0;
}
