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
/*  Last update: Dec 31, 2005                                                */
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
#include "log.h"

static BSolver  * bsolver;
static PlotDefine *plotdef;
///////////////////////////////////////////////////////////////////////////////
/*-----------------------------------------------------------------------
 * Function:     BSolver::FillPoly
 * Scope:        Private
 * Description:  This function draws and fills a polygon.  The polygon is
 *               limited to 32 sides
 */
void BSolver::FillPoly (int c, double *adnode, int nFColor, int nBColor)
{
  // return if the number of sides is greater than 32
  if (c > 32) return;

  // draw and fill the polygon
  GRSetColor(nFColor);
  GRSetFillColor(nBColor);
  int an[64];
  for(int i = 0; i < c; ++i)
  {
    an[2*i]   = xMAP(adnode[2*i]);
    an[2*i+1] = yMAP(adnode[2*i+1]);
  }
  GRFillPoly(c, an);
}



/*-----------------------------------------------------------------------
 * Function:     BSolver::plot_mesh_redraw
 * Scope:        Private
 * Description:  This function redraw the mesh to screen.
 */
void BSolver::plot_mesh_redraw(void *dummy)
{
  PlotDefine *plot = (PlotDefine *)dummy;
  // initialize local variables
  int cte = 0;
  int cre = 0;
  int cObtuse = 0;
  int x=3*wx/4,y=10,d=20;
  cx = wy - 2*cxBorder - 23;
  cy = wy - 2*cyBorder - 23;
  char szGlobalBuf[256];
  double adnode[8];

  //draw triangle from zone level
  for(int z = 0; z < zone_num; z++)
    for(int itri = 0; itri < zone[z].datri.size(); itri++)
    {
      Tri ptri = zone[z].datri[itri];
      adnode[0] = zone[z].danode[ptri.node[0]].x;
      adnode[1] = zone[z].danode[ptri.node[0]].y;
      adnode[2] = zone[z].danode[ptri.node[1]].x;
      adnode[3] = zone[z].danode[ptri.node[1]].y;
      adnode[4] = zone[z].danode[ptri.node[2]].x;
      adnode[5] = zone[z].danode[ptri.node[2]].y;

      int Fcolor=GR_WHITE;
      int Bcolor;
      switch (zonedata[z]->material_type)
      {
      case Semiconductor:
      {
        SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
        double doping = pzonedata->aux[ptri.node[0]].Nd - pzonedata->aux[ptri.node[0]].Na +
                        pzonedata->aux[ptri.node[1]].Nd - pzonedata->aux[ptri.node[1]].Na +
                        pzonedata->aux[ptri.node[2]].Nd - pzonedata->aux[ptri.node[2]].Na;
        if(doping>0) //N-Type
          Bcolor= GR_LIGHTGREEN;
        else         //P-Type
          Bcolor= Rainbow_15;
        break;
      }
      case Insulator:
      {
          if (zonedata[z]->material==SiO2)
            Bcolor= GR_CYAN;//Rainbow_3;
          else if (zonedata[z]->material==Nitride)
            Bcolor= GR_MAGENTA;
	  else   Bcolor= Rainbow_7;
          break;
      }
      case Conductor:
        Bcolor= GR_BLUE;
        break;
      case Vacuum:
        Bcolor= GR_LIGHTCYAN;
        break;
      case PML:
        Bcolor= Rainbow_3;
        break;
      }
      FillPoly(3, adnode, Fcolor, Bcolor);
    }

  //draw Voronoi Cell
  for(int z = 0; z < zone_num; z++)
  {
    GRSetColor(GR_BLACK);
    GRSetFillColor(GR_BLACK);
    GRSetColor(GR_WHITE);
    sprintf(szGlobalBuf, "zone %d lable : %s ",z+1,zone[z].zonename);
    GRText(x+10, y+=d,    szGlobalBuf);
    sprintf(szGlobalBuf, "   %d nodes",int(zone[z].danode.size()));
    GRText(x+10, y+=d,    szGlobalBuf);
    sprintf(szGlobalBuf, "   %d elements",int(zone[z].datri.size()));
    GRText(x+10, y+=d, szGlobalBuf);
  }
#ifdef HAVE_X11
  GRSetColor(GR_WHITE);
  sprintf(szGlobalBuf, "Press ESC to exit");
  GRText(x+10, y+40, szGlobalBuf);
#endif
}


void Redraw2D()
{
  bsolver->plot_mesh_redraw(plotdef);
}

/*-----------------------------------------------------------------------
 * Function:     BSolver::draw_mesh_screen
 * Scope:        public
 * Description:  This function draw the mesh to screen
 */
void BSolver::plot_mesh_screen(PlotDefine &plot)
{

  int fDone = 0;
  int x,y; //location of pointer
  bsolver = this;
  plotdef = &plot;

  if(plot.Resolution == RES_640x480)
  {wx = 640; wy = 480;}
  else if(plot.Resolution == RES_1024x768)
  {wx = 1024; wy = 768;}
  else if(plot.Resolution == RES_800x600)
  {wx = 800; wy = 600;}

  short int keycode;
  int arg1, arg2;
  xDelta = xMax - xMin;
  yDelta = yMax - yMin;

  if (xDelta > yDelta) yDelta = xDelta;
  else           xDelta = yDelta;
  cxBorder = 10;
  cyBorder = 10;
  cx = wy - 2*cxBorder - 23;
  cy = wy - 2*cyBorder - 23;
  GRInitGraphics();

#ifdef HAVE_WIN32
  GROpenGraphWin("Mesh Topology", "Mesh", -1, -1, wx, wy, Redraw2D);
#endif

#ifdef HAVE_X11
  GROpenGraphWin("Mesh Topology", "Mesh", -1, -1, wx, wy, NULL);


  while (!fDone)
  {
    switch(getevent(&keycode, &arg1, &arg2,&x,&y))
    {
    case EXPOSE:
      plot_mesh_redraw(&plot);
      flushdisplay();
      break;

    case CONFIGURENOTIFY:
      if(wx!=arg1 || wy!=arg2)
      {
        wx=arg1;
        wy=arg2;
        GRClearGraphWin();
        plot_mesh_redraw(&plot);
        flushdisplay();
      }
      break;

    case BUTTONPRESS:
      break;
    case KEYPRESS:
      if(keycode == ESCAPE)     fDone = 1;
      break;
    }
  }
  if(plot.GenerateTIFF)
  {
    gss_log.string_buf()<<"Saving Mesh Plot to TIFF File "<<plot.TIFFFileName<<"...";
    gss_log.record();
    GRSaveScreen (plot.TIFFFileName.c_str(), wx, wy);
    gss_log.string_buf()<<"ok"<<endl;
    gss_log.record();
  }

#endif

  GRFreeGraphics();
}


/*-----------------------------------------------------------------------
 * Function:     BSolver::draw_mesh_ps
 * Scope:        public
 * Description:  This function generates a postscript file.
 */
void BSolver::plot_mesh_ps(PlotDefine &plot)
{
  cx = 9360;
  cy = 9360;
  cxBorder = 1440;
  cyBorder = 3240;
  double adnode[8];

  xDelta = xMax - xMin;
  yDelta = yMax - yMin;

  if (xDelta > yDelta) yDelta = xDelta;
  else           xDelta = yDelta;
  gss_log.string_buf()<<"Saving Mesh Plot to PS File "<<plot.PSFileName<<"...";
  gss_log.record();
  GRPrintOpen(plot.PSFileName.c_str(),20);
  GRSetColor(GR_BLACK);

  //draw triangle from zone level
  for(int i = 0; i < zone_num; i++)
    for(int itri = 0; itri < zone[i].datri.size(); itri++)
    {
      Tri ptri = zone[i].datri[itri];
      adnode[0] = zone[i].danode[ptri.node[0]].x;
      adnode[1] = zone[i].danode[ptri.node[0]].y;
      adnode[2] = zone[i].danode[ptri.node[1]].x;
      adnode[3] = zone[i].danode[ptri.node[1]].y;
      adnode[4] = zone[i].danode[ptri.node[2]].x;
      adnode[5] = zone[i].danode[ptri.node[2]].y;

      int Fcolor=GR_WHITE,Bcolor;
      if(zonedata[i]->material_type==Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[i]);
        double doping = pzonedata->aux[ptri.node[0]].Nd - pzonedata->aux[ptri.node[0]].Na +
                        pzonedata->aux[ptri.node[1]].Nd - pzonedata->aux[ptri.node[1]].Na +
                        pzonedata->aux[ptri.node[2]].Nd - pzonedata->aux[ptri.node[2]].Na;
        if(doping>0) //N-Type
          Bcolor= GR_LIGHTGREEN;
        else         //P-Type
          Bcolor= GR_LIGHTMAGENTA;
      }
      else if (zonedata[i]->material==SiO2)
        Bcolor= GR_CYAN;//Rainbow_3;
      else if (zonedata[i]->material==Nitride)
        Bcolor= GR_MAGENTA;
      else if (zonedata[i]->material_type==Conductor)
        Bcolor= GR_BLUE;
      else if (zonedata[i]->material_type==Vacuum)
        Bcolor= GR_LIGHTCYAN;
      else if (zonedata[i]->material_type==PML)
        Bcolor= Rainbow_9;


      FillPoly(3, adnode, Fcolor, Bcolor);
    }
  //draw Voronoi Cell
  for(int z = 0; z < zone_num; z++)
  {
    GRSetColor(GR_BLACK);
    GRSetFillColor(GR_BLACK);

    for(int j = 0; j < zone[z].davedge.size(); j++)
    {
        VEdge Edge = zone[z].davedge[j];
        int sx1=xMAP(Edge.x1);
        int sy1=yMAP(Edge.y1);
        int sx2=xMAP(Edge.x2);
        int sy2=yMAP(Edge.y2);
        GRLine(sx1,sy1,sx2,sy2);
    }
  }
  GRPrintClose();
  gss_log.string_buf()<<"ok"<<endl;
  gss_log.record();
}

#endif
