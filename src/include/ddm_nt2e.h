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
/*  Last update: Jan 19, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _ddm_nt2e_h_
#define _ddm_nt2e_h_

#include "petscts.h"
#include "petscsnes.h"
#include "bsolver.h"


class DDM_Solver_L2E : virtual public BSolver
{
private:
  vector<int>    zofs;         // the offset value of each zone in solution vector
  PetscTruth     mf_flg;
  SNESConvergedReason reason;
  void    solver_pre_compute();
  void    solution_update();
  void    diverged_recovery();
  void    LET_norm_estimat(PetscScalar &);
public:
  int            N;            // the scale of matrix
  SNES           snes;         // nonlinear solver context
  KSP            ksp;          // linear solver context
  PC             pc;           // preconditioner context
  Vec            x,r;          // solution, residual vectors
  Vec            x_n,x_n1,x_n2; // privious solution vector for LET evaluation
  Vec            xp;            // predict solution vector;
  Vec            LTE;           // local truncation error vector;
  Mat            J,JPrec,JTmp; // Jacobian matrix
  PetscInt       its;          // iteration number
  PetscReal      norm;
  
  PetscReal      relative_toler;
  PetscReal      toler_relax;
  PetscReal      possion_abs_toler;
  PetscReal      elec_continuty_abs_toler;
  PetscReal      hole_continuty_abs_toler;
  PetscReal      heat_equation_abs_toler;
  PetscReal      electrode_abs_toler;
  
  PetscReal      potential_norm;
  PetscReal      electron_norm;
  PetscReal      hole_norm;
  PetscReal      temperature_norm;
  PetscReal      possion_norm;
  PetscReal      elec_continuty_norm;
  PetscReal      hole_continuty_norm;
  PetscReal      heat_equation_norm;
  PetscReal      electrode_norm;

public:
  int     init_solver(SolveDefine &s)                     ;
  int     do_solve(SolveDefine &s)                        ;
  int     solve_equ(SolveDefine &s)                       ;
  int     solve_steadystate(SolveDefine &s)               ;
  int     solve_dcsweep(SolveDefine &s)                   ;
  int     solve_iv_trace(SolveDefine &s)                  ;
  int     solve_transient(SolveDefine &s)                 ;
  int     destroy_solver(SolveDefine &s)                  ;
  void    form_function_pn_2E(PetscScalar *,PetscScalar *);
  void    form_jacobian_pn_2E(PetscScalar *, Mat *, Mat *);
  void    error_norm_estimat(PetscScalar *,PetscScalar *)   ;
  DDM_Solver_L2E():N(0),its(0),norm(0.0),mf_flg(PETSC_FALSE){}
                                       
};


#endif
