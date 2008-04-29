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
/*  Last update: March 29, 2007                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _mix1_h_
#define _mix1_h_

#include "petscts.h"
#include "petscsnes.h"
#include "bsolver.h"
#include "ndevexch.h" 
#include "mixcomm.h"

class DDM_Mix_Solver_L1E : virtual public BSolver
{
private:
  vector<int>    zofs;         // the offset value of each zone in solution vector
  SolveDefine    *psv;
  int listener; // Our listening socket. 
  int client;   // The current client's socket. 
  sCKTinfo      CKTInfo;
  sDeviceinfo   Deviceinfo;   // the numerical device information passed from ngspice
  sPINinfo      PINinfos[7];  // the pin voltage and conduction matrix enties.
  sPINcond      PINconds[7];  // structure for calculate pin conductance matrix entries.
public:
  int            N;            // the scale if matrix
  SNES           snes;         // nonlinear solver context
  KSP            ksp;          // linear solver context
  PC             pc;           // preconditioner context
  Vec            x,r;          // solution, residual vectors
  Mat            J,JTmp;       // Jacobian matrix
  PetscInt       its;          // iteration number
  PetscReal      norm;
  SNESConvergedReason reason;
  
  PetscReal relative_toler;
  PetscReal toler_relax;
  PetscReal possion_abs_toler;
  PetscReal elec_continuty_abs_toler;
  PetscReal hole_continuty_abs_toler;
  
  PetscReal potential_norm;
  PetscReal electron_norm;
  PetscReal hole_norm;
  PetscReal possion_norm;
  PetscReal elec_continuty_norm;
  PetscReal hole_continuty_norm;
 
  
public:
  int     init_solver(SolveDefine &s) ;
  int     do_solve(SolveDefine &s)    ;
  int     destroy_solver(SolveDefine &s) ;
  
  DDM_Mix_Solver_L1E():N(0),its(0),norm(0.0) {} ;
  ~DDM_Mix_Solver_L1E() {};
  void form_function_pn_Mix1(PetscScalar *x,PetscScalar *f);
  void form_jacobian_pn_Mix1(PetscScalar *x, Mat *jac, Mat *jtmp);
  void error_norm_pn_Mix1(PetscScalar *x,PetscScalar *f);
  int  DEV_LOAD();
  int  DEV_ACCEPT();
  int  DEV_CONV_TEST();
  int  tran_solve();
  int  dc_solve();
  void solution_update();
  void time_back_recovery();
  void diverged_recovery();
};




#endif
