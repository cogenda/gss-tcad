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
#ifndef _solver_h_
#define _solver_h_

#include "config.h"
#include "doping.h"
#include "mesh.h"
#include "bc.h"
#include "zonedata.h"
#include "solverdef.h"
#include "cmdbuf.h"


class BSolver : public MESH
{
public:
  char                  ModelFile[32];
  double                LatticeTemp;
  double                DeviceDepth;
  int                   Carrier;
  CmdBuf                *pcmdbuf;
  PhysicalUnitScale     scale_unit;
  vector<BZoneData*>    zonedata;
  vector<VSource*>      vsrc;                     // dynamic array for vsource pointer
  vector<ISource*>      isrc;                     // dynamic array for isource pointer
  LSource*              lsrc;                     // pointer for light source
  vector<DopingFunc *>  doping_func;
  DABC                  bc;
  SolveDefine           solve_define;
  PMISDefine            PMIS_define;
  ExportDefine          export_define;
  ImportDefine          import_define;
  AttachDefine          attach_define;
  RefineDefine          refine_define;
  PlotDefine            plot_define;
  vector<ProbeDefine>   probe_define;
  int                   static_config(CmdBuf *cmdbuf);
  void                  show_config_info();
public:
  PetscScalar           clock;
  ODE_Formula           ODE_formula;
public:
  double      set_refine_scale(Tri &, RefineDefine &);

  int         init_esource(list<Cmd> &cmdlist);
  int         init_grid(list<Cmd> &cmdlist);
  int         init_doping_func();
  int         build_zonedata();
  int         build_least_squares();

  int         import_doping_from_cgns(char* );
  int         import_mole_from_cgns(char *);

  int         setup_doping();
  int         setup_init_data();
  int         extract_ascii(char *);
  int         setup_bc();
  void        reorder();
  int         refine(RefineDefine &rf);
  void        clear_data();
  
  int         optgen_update(PetscScalar currentT); //update optical generation rate 
  
  int         probe_open(int);
  int         probe(int , PetscScalar  );
  int         probe_close();
public:
#if (defined(HAVE_X11) || defined(HAVE_WIN32))
  int         wx;                       // window size
  int         wy;                       // window size
  int         cx;                       // plot window size
  int         cy;                       // plot window size
  int         cxBorder;                 // border in x direction
  int         cyBorder;                 // border in y direction
  double      xDelta;                   // delta x spacing
  double      yDelta;                   // delta y spacing
  int         xMAP (double x) {return (int) (cx*(0.5+(x-(xMin+xMax)/2)/xDelta))+cxBorder;}
  int         yMAP (double y) {return (int) (cy*(0.5-(y-(yMax+yMin)/2)/yDelta))+cyBorder;}
  void        FillPoly (int,double*,int,int);
  void        plot_mesh_redraw(void *);
  void        plot_mesh_screen(PlotDefine &plot_define);
  void        plot_mesh_ps(PlotDefine &plot_define);
  void        show_mesh(PlotDefine &plot_define);
  void        show_mesh_redwraw();
  void        plot_data(PlotDefine &plot_define);
#endif

  int        vtk_display_data(PlotDefine & plot_define);
  int        vtk_output_file(char *);

public:
  BSolver():lsrc(0),clock(0) {}
  ~BSolver();
}
;

#endif

