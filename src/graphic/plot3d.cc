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
/*  Last update: Dec 31, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "config.h"
#if (defined(HAVE_WIN32) || defined(HAVE_X11))

#include "bsolver.h"
#include "xgraph.h"
#include "wgraph.h"
#include "grafix3d.h"
#include "log.h"

static PlotDefine *plotdef;
static TriPlot3D* triplot;
static BSolver  * bsolver;

inline double FunLinear(PetscScalar value)
{
  return value;
}

inline double FunSignedLog10(PetscScalar value)
{
  return dsign(value)*log10(1.0+fabs(value));
}

void TitleEtc ( PlotDefine &plot, int x1, int y1, int ResFactor )
{
  char txt[64];
  int ix;
  int iy = int(CHARPIXELHIGH * ResFactor * 1.2);

  if      (plot.Resolution == RES_640x480)  ix = 56;
  else if (plot.Resolution == RES_1024x768) ix = 64;
  else if (plot.Resolution == RES_800x600)  ix = 96;
  //set the unit for each variable
  string unit;
  switch(plot.Variable)
  {
  case DopingNd:
  case DopingNa:
  case Doping:
  case ElecDensity:
  case HoleDensity:
    unit = "(cm^-3)"; break;
  case ElecTemp:
  case HoleTemp:
  case Temperature:
    unit = "(K)"; break;
  case Potential:
  case Phi_Intrinsic:
  case PhiP:
  case PhiN:
    unit = "(V)"; break;
  case EFieldX:
  case EFieldY:
  case OpticalEx:
  case OpticalEy:
  case OpticalEz:
    unit = "(V/cm)"; break;
  case QuantumEc:
  case QuantumEv:
    unit = "(eV)"; break;
  case OpticalG:
    unit = "(cm^-3 s^-1)"; break;
  case OpticalHx:
  case OpticalHy:
  case OpticalHz:
    unit = "(A/cm)"; break;
  }

  x1 += 4 * ResFactor * CHARPIXELWIDE;
  y1 *= ResFactor;
  sprintf(txt, "GSS 3D graphic output");
  GRText(ix + x1, y1 += iy, txt);
  sprintf(txt, "Variable : %s %s",plot.VariableName.c_str(),unit.c_str());
  GRText(ix + x1, y1 += iy, txt);
  sprintf(txt, "AZ = %e, EL = %e",plot.az, plot.el);
  GRText(ix + x1, y1 += iy, txt);

#ifdef HAVE_X11
  sprintf(txt, "Press ESC to exit");
  GRText(ix + x1, y1 += iy, txt);
#endif
}

void Redraw3D()
{
  triplot->PlotGrid();
  triplot->PlotMesh();
  triplot->PlotData(TRUE, TRUE);
  TitleEtc(*plotdef,0, 5, 1);
}


/*-----------------------------------------------------------------------
 * Function:     BSolver::plot_data
 * Scope:        Public
 * Description:  This function draws 3D data and mesh on the screen
 */
void BSolver::plot_data(PlotDefine &plot)
{
  int    MAXx,MAXy;
  int    x,y; //location of pointer
  int    x_old,y_old; //last location of pointer
  AINODE *aainode=0,*painode=0;
  PNT3D  *apnt3dInput=0;
  int    nodes=0,triangles=0;
  int    daz=10,del=5;
  double (*Measure) (PetscScalar);
  plotdef = &plot;

  if(plot.Measure==Linear)
    Measure = FunLinear;
  if(plot.Measure==SignedLog)
    Measure = FunSignedLog10;

  //data for semiconductor zone only
  if(plot.Variable==Doping||plot.Variable==DopingNd||plot.Variable==DopingNa||
      plot.Variable==ElecDensity||plot.Variable==HoleDensity||
      plot.Variable==ElecTemp||plot.Variable==HoleTemp||
      plot.Variable==PhiP||plot.Variable==PhiN||plot.Variable==Phi_Intrinsic||
      plot.Variable==QuantumEc||plot.Variable==QuantumEv||
      plot.Variable==OpticalG)
  {
    for(int z=0;z<zone_num;z++)
      if(zonedata[z]->material_type == Semiconductor)
      {
        nodes += zone[z].danode.size();
        triangles += zone[z].datri.size();
      }
    aainode     = new AINODE[triangles];
    apnt3dInput = new PNT3D[nodes];
    painode=aainode;
    int offset = 0;
    for(int z=0;z<zone_num;z++)
    {
      if(zonedata[z]->material_type == Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast<SMCZone *>(zonedata[z]);
        for(int i=0;i<zone[z].datri.size();i++)
        {
          (*painode)[0] = zone[z].datri[i].node[0] + offset;
          (*painode)[1] = zone[z].datri[i].node[1] + offset;
          (*painode)[2] = zone[z].datri[i].node[2] + offset;
          (*painode)[3] = -1u;

          painode++;
        }

        for(int j=0;j<zone[z].danode.size();j++)
        {
          apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
          apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
          switch(plot.Variable)
          {
          case DopingNd:
            apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Nd*pow(scale_unit.s_centimeter,3));break;
          case DopingNa:
            apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Na*pow(scale_unit.s_centimeter,3));break;
          case Doping:
            apnt3dInput[j+offset].z = Measure((pzonedata->aux[j].Nd-pzonedata->aux[j].Na)*pow(scale_unit.s_centimeter,3));break;
          case ElecDensity:
            apnt3dInput[j+offset].z = Measure(fabs(pzonedata->fs[j].n)*pow(scale_unit.s_centimeter,3));break;
          case HoleDensity:
            apnt3dInput[j+offset].z = Measure(fabs(pzonedata->fs[j].p)*pow(scale_unit.s_centimeter,3));break;
          case ElecTemp:
            apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].Tn/scale_unit.s_kelvin);break;
          case HoleTemp:
            apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].Tp/scale_unit.s_kelvin);break;
          case Phi_Intrinsic:
            apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].phi_intrinsic/scale_unit.s_volt);break;
          case PhiP:
            apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].phip/scale_unit.s_volt);break;
          case PhiN:
            apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].phin/scale_unit.s_volt);break;
          case QuantumEc:
            apnt3dInput[j+offset].z = Measure((pzonedata->fs[j].Eqc+pzonedata->fs[j].P+pzonedata->aux[j].affinity)/scale_unit.s_volt);break;
          case QuantumEv:
            apnt3dInput[j+offset].z = Measure((pzonedata->fs[j].Eqv+pzonedata->fs[j].P+pzonedata->aux[j].affinity+pzonedata->aux[j].Eg)/scale_unit.s_volt); break;
          case OpticalG:
            apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OptG*pow(scale_unit.s_centimeter,3)*scale_unit.s_second); break;
          }
        }
        offset+=zone[z].danode.size();
      }
    }
  }

  // data for semiconductor zone, insulator zone and electrode zone
  if(plot.Variable==Potential||plot.Variable==EFieldX||plot.Variable==EFieldY||plot.Variable==Temperature)
  {
    for(int z=0;z<zone_num;z++)
    {
      if(zonedata[z]->material_type == Semiconductor ||
          zonedata[z]->material_type == Insulator ||
          zonedata[z]->material_type == Conductor)
      {
        nodes += zone[z].danode.size();
        triangles += zone[z].datri.size();
      }
    }
    aainode     = new AINODE[triangles];
    apnt3dInput = new PNT3D[nodes];
    painode=aainode;
    int offset = 0;
    for(int z=0;z<zone_num;z++)
    {
      if(zonedata[z]->material_type == Semiconductor ||
          zonedata[z]->material_type == Insulator ||
          zonedata[z]->material_type == Conductor)
      {
        for(int i=0;i<zone[z].datri.size();i++)
        {
          (*painode)[0] = zone[z].datri[i].node[0] + offset;
          (*painode)[1] = zone[z].datri[i].node[1] + offset;
          (*painode)[2] = zone[z].datri[i].node[2] + offset;
          (*painode)[3] = -1u;
          painode++;
        }

        if(zonedata[z]->material_type == Semiconductor)
        {
          SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case Potential:
              apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].P/scale_unit.s_volt);break;
            case EFieldX:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Ex/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case EFieldY:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Ey/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case Temperature:
              apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].T/scale_unit.s_kelvin);break;
            }
          }
        }

        if(zonedata[z]->material_type == Insulator)
        {
          ISZone *pzonedata = dynamic_cast<ISZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case Potential:
              apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].P/scale_unit.s_volt);break;
            case EFieldX:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Ex/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case EFieldY:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Ey/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case Temperature:
              apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].T/scale_unit.s_kelvin);break;
            }
          }
        }

        if(zonedata[z]->material_type == Conductor)
        {
          ElZone *pzonedata = dynamic_cast<ElZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case Potential:
              apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].P/scale_unit.s_volt);break;
            case EFieldX:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Ex/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case EFieldY:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].Ey/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case Temperature:
              apnt3dInput[j+offset].z = Measure(pzonedata->fs[j].T/scale_unit.s_kelvin);break;
            }
          }
        }
        offset += zone[z].danode.size();
      }
    }
  }

  // data for all zones
  if(plot.Variable==OpticalEx||plot.Variable==OpticalEy||plot.Variable==OpticalEz||plot.Variable==OpticalHx||
      plot.Variable==OpticalHy||plot.Variable==OpticalHz)
  {
    for(int z=0;z<zone_num;z++)
    {
      if(zonedata[z]->material_type == Semiconductor ||
          zonedata[z]->material_type == Insulator||
          zonedata[z]->material_type == Conductor||
          zonedata[z]->material_type == Vacuum||
          zonedata[z]->material_type == PML)
      {
        nodes += zone[z].danode.size();
        triangles += zone[z].datri.size();
      }
    }
    aainode     = new AINODE[triangles];
    apnt3dInput = new PNT3D[nodes];
    painode=aainode;
    int offset = 0;
    for(int z=0;z<zone_num;z++)
    {
      if(zonedata[z]->material_type == Semiconductor || zonedata[z]->material_type == Insulator||
          zonedata[z]->material_type==Conductor||zonedata[z]->material_type==Vacuum||zonedata[z]->material_type==PML)
      {
        for(int i=0;i<zone[z].datri.size();i++)
        {
          (*painode)[0] = zone[z].datri[i].node[0] + offset;
          (*painode)[1] = zone[z].datri[i].node[1] + offset;
          (*painode)[2] = zone[z].datri[i].node[2] + offset;
          (*painode)[3] = -1u;
          painode++;
        }

        if(zonedata[z]->material_type == Semiconductor)
        {
          SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case OpticalEx://The scale parameter can be found in phy_scale.cc```
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEx/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEy/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEz/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalHx:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHx/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHy/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHz/scale_unit.s_A*scale_unit.s_centimeter);break;
            }
          }
        }

        if(zonedata[z]->material_type == Insulator)
        {
          ISZone *pzonedata = dynamic_cast<ISZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case OpticalEx://The scale parameter can be found in phy_scale.cc```
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEx/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEy/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEz/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalHx:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHx/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHy/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHz/scale_unit.s_A*scale_unit.s_centimeter);break;
            }
          }
        }
        if(zonedata[z]->material_type==Conductor)
        {
          ElZone *pzonedata = dynamic_cast<ElZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case OpticalEx://The scale parameter can be found in phy_scale.cc```
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEx/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEy/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEz/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalHx:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHx/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHy/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHz/scale_unit.s_A*scale_unit.s_centimeter);break;
            }
          }
        }
        if(zonedata[z]->material_type==Vacuum)
        {
          VacuumZone *pzonedata = dynamic_cast<VacuumZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case OpticalEx://The scale parameter can be found in phy_scale.cc```
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEx/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEy/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEz/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalHx:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHx/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHy/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHz/scale_unit.s_A*scale_unit.s_centimeter);break;
            }
          }
        }
        if(zonedata[z]->material_type==PML)
        {
          PMLZone *pzonedata = dynamic_cast<PMLZone *> (zonedata[z]);
          for(int j=0;j<zone[z].danode.size();j++)
          {
            apnt3dInput[j+offset].x = zone[z].danode[j].x*1e4/scale_unit.s_centimeter;
            apnt3dInput[j+offset].y = zone[z].danode[j].y*1e4/scale_unit.s_centimeter;
            switch(plot.Variable)
            {
            case OpticalEx://The scale parameter can be found in phy_scale.cc```
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEx/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEy/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalEz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpEz/scale_unit.s_volt*scale_unit.s_centimeter);break;
            case OpticalHx:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHx/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHy:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHy/scale_unit.s_A*scale_unit.s_centimeter);break;
            case OpticalHz:
              apnt3dInput[j+offset].z = Measure(pzonedata->aux[j].OpHz/scale_unit.s_A*scale_unit.s_centimeter);break;
            }
          }
        }
        offset += zone[z].danode.size();
      }
    }
  }

  if(plot.Resolution == RES_640x480)
  {MAXx = 640; MAXy = 480;}
  else if(plot.Resolution == RES_1024x768)
  {MAXx = 1024; MAXy = 768;}
  else if(plot.Resolution == RES_800x600)
  {MAXx = 800; MAXy = 600;}

  triplot = new TriPlot3D(triangles, nodes,
                          aainode, apnt3dInput,
                          plot.az, plot.el,
                          plot.Resolution, plot.Style,
                          50, -1, 0,
                          plot.Measure==SignedLog);

#ifdef HAVE_X11
  GROpenGraphWin("3D Surface Plot", "Plot", 0, 0, MAXx, MAXy, NULL);
  int fDone = FALSE;
  short int keycode;
  int arg1, arg2;
  while (!fDone)
  {
    switch(getevent(&keycode, &arg1, &arg2,&x,&y))
    {
    case EXPOSE:
      triplot->PlotGrid();
      triplot->PlotMesh();
      triplot->PlotData(TRUE, TRUE);
      TitleEtc(plot,0, 5, 1);
      flushdisplay();
      break;

    case CONFIGURENOTIFY:
      break;
    case BUTTONPRESS:
      x_old=x;
      y_old=y;
    case BUTTONRELEASE:
      break;
    case BUTTONMOTION:
      daz =  isign(x_old-x)*int(sqrt(fabs((x_old-x)/1.0)));
      del =  isign(y_old-y)*int(sqrt(fabs((y_old-y)/1.0)));
      GRClearGraphWin();
      triplot->GetViewAspect(plot.az+=daz, plot.el+=del);
      triplot->Xform3DData();             // Transform the data 3D->3D
      triplot->Xform2DData();             // Transform the data 3D->2D
      triplot->DeterminePlotOrder();      // Sort the order of appearance
      triplot->PlotGrid();
      triplot->PlotMesh();
      triplot->PlotData(TRUE, TRUE);
      TitleEtc(plot,0, 5, 1);
      flushdisplay();
      break;
    case KEYPRESS:
      if(keycode == ESCAPE)   fDone = 1;
      break;
    }
  }
  if(plot.GenerateTIFF)
  {
    gss_log.string_buf()<<"Saving Screen to TIFF File "<<plot.TIFFFileName<<"...";
    gss_log.record();
    if(GRSaveScreen (plot.TIFFFileName.c_str(), MAXx, MAXy))
    {
      gss_log.string_buf()<<"Open TIFF File failed."<<endl;
      gss_log.record();
    }
    else
    {
      gss_log.string_buf()<<"ok"<<endl;
      gss_log.record();
    }
  }
  GRFreeGraphics();
#endif

#ifdef HAVE_WIN32
  GROpenGraphWin("3D Surface Plot", "Plot", 0, 0, MAXx, MAXy, Redraw3D);
#endif


  if (plot.GeneratePS)
  {
    gss_log.string_buf()<<"Saving Screen to PS File "<<plot.PSFileName<<"...";
    gss_log.record();
    triplot->nResFactor = 20;
    triplot->cx = 468 * triplot->nResFactor;
    triplot->cy = 468 * triplot->nResFactor;
    triplot->cxLBorder = 72  * triplot->nResFactor;
    triplot->cyLBorder = 162 * triplot->nResFactor;
    triplot->Xform2DData();
    if(GRPrintOpen(plot.PSFileName.c_str(), triplot->nResFactor))
    {
      gss_log.string_buf()<<"Open PS File failed."<<endl;
      gss_log.record();
    }
    else
    {
      triplot->PlotGrid();
      triplot->PlotMesh();
      triplot->PlotData(TRUE, TRUE);
      TitleEtc(plot, 0, 15, triplot->nResFactor);
      GRPrintClose();
      gss_log.string_buf()<<"ok"<<endl;
      gss_log.record();
    }
  }

  delete triplot;
  delete [] apnt3dInput;
  delete [] aainode;

}

#endif
