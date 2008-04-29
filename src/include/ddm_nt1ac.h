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
/*  Last update: Oct 25, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _ddm_nt1ac_h_
#define _ddm_nt1ac_h_
#include "petscts.h"
#include "bsolver.h"


class DDM_Solver_L1AC : virtual public BSolver
{
private:
  vector<int>    zofs;         // the offset value of each zone in solution vector
public:
  int            N;            // the size of real jac matrix
  KSP            ksp;          // linear solver context
  PC             pc;           // preconditioner context
  PetscScalar   *x;            // DC solution;
  Vec            v,r;          // AC solution vector and right vector
  Mat            J,JTmp;       // Jacobian matrix
public:
  int     init_solver(SolveDefine &s) ;
  int     do_solve(SolveDefine &s);
  int     solve_ac(SolveDefine &s);
  void    compute_electrode_current(PetscScalar);
  int     destroy_solver(SolveDefine &s);
  DDM_Solver_L1AC():N(0){};
  void    Get_Jacobian();
  void    form_linear_system_ac1(PetscScalar, PetscScalar *, Mat *, Vec *, Mat *);
};


#endif

