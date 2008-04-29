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

#ifndef _ctrl_h_
#define _ctrl_h_
#include "ddm_nt1e.h"
#include "ddm_nt1ac.h"
#include "qddm_nt1e.h"
#include "mix1.h"
#include "mix2.h"
#include "ddm_nt2e.h"
#include "ebm_nt3e.h"
#include "em_fem.h"

class SolveControl: public DDM_Solver_L1E,
                    public DDM_Solver_L1AC,
                    public QDDM_Solver_L1E,
                    public DDM_Mix_Solver_L1E,
                    public DDM_Solver_L2E,
                    public DDM_Mix_Solver_L2E,
                    public EBM_Solver_L3E,
                    public EM_FEM_Solver
{
private:
  int  set_method(list<Cmd>::iterator );
  int  set_solve(list<Cmd>::iterator  );
  int  set_export(list<Cmd>::iterator);
  int  set_import(list<Cmd>::iterator );
  int  set_attach(list<Cmd>::iterator );
  int  set_refine(list<Cmd>::iterator );
  int  set_plot(list<Cmd>::iterator );
  int  set_plotmesh(list<Cmd>::iterator );
  int  set_vtk_plot(list<Cmd>::iterator );
  int  set_probe(list<Cmd>::iterator ) ;
  int  do_solve();
  int  do_refine();
  int  do_export();
  int  do_import();
  int  do_attach();
  int  do_plot();
  int  do_plotmesh();
  int  do_vtk_plot();

public:
  int   run_control(list<Cmd>&);
public:
  SolveControl()  {}
};

#endif
