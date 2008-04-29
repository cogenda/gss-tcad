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

#ifndef _em_fem_h_
#define _em_fem_h_

#include "mathfunc.h"
#include "petscksp.h"
#include "bsolver.h"
#include "mesh.h"
#include "cmdbuf.h"
#include <vector>
using namespace std;

class EM_FEM_Solver : virtual public BSolver
{
private:
  int            N;            // the scale if matrix
  int            total_nodes;  // node number of all zones
  PetscScalar    v;            //insert variable
  PetscInt       its;          // iteration number
  PetscInt       *rows,count; //variables for add dirichlet boundary condition
  PetscReal      norm;
   
  PetscErrorCode     ierr;
  KSPConvergedReason reason;  
  KSP                solver;
  PC                 prec; 
  //vector<int>    zofs;         // the offset value of each zone in solution vector
  
  int         eleNode[3];        //store nodes' global number of each element
  PetscScalar eleA[3];           //assistant variables
  PetscScalar eleB[3];
  PetscScalar eleC[3];
  PetscScalar eleMatrix[3][3];   //element matrix  
  PetscScalar eleV[3];           //element right vector  
   
  PetscScalar RealAlphaX;        //only used for PML zone
  PetscScalar RealAlphaY;
  PetscScalar RealBeta;
  PetscScalar ImageAlphaX;
  PetscScalar ImageAlphaY;
  PetscScalar ImageBeta;
  PetscScalar Electri0;          //only used for PML zone  
  PetscScalar Magne0;  
  PetscScalar TransVar; 
  PetscScalar CurrentX; 
  PetscScalar CurrentY; 
     
  PetscScalar RealPermi;         //Real part of permitivity
  PetscScalar ImagePermi;        //Image part ```````
  PetscScalar InitPhi;
  PetscScalar PhiDiff;
  PetscScalar SourceAlpha;
   
public:   
  Vec            x,b;      // solution and right vectors
  //Vectors to save Ex,Ey,Hx,Hy in complex form
  vector<PetscScalar>         TMEZ,TMHX,TMHY,TEHZ,TEEX,TEEY; 
  PetscScalar    *TM;      //Needle of TM mode for output the result 
  PetscScalar    *TE;
  Mat            A;        // matrix
  PetscScalar XLeft;       //minimun value of x in the whole model
  PetscScalar XRight;
  PetscScalar YTop;
  PetscScalar YBottom;
  PetscScalar VacXMax;     //Find max X in vacuum zone
  PetscScalar VacXMin;
  PetscScalar VacYMax;
  PetscScalar VacYMin;     //These 4 variable is used to find 4 corners of vacuum zone
  PetscScalar lamda;       //wave length of incidence light
  PetscScalar power;       //intensity of incidence light 
  PetscScalar WeightOfTE;  //weight of TE mode in incidence light
  PetscScalar WeightOfTM;
  PetscScalar OptAngle;    //incident angle of light --degree
  PetscScalar OptTheta;    //incident angle of light --radian
  PetscScalar eta;         //quantum efficiency of semiconductor material 
  PetscScalar WaveVector0; //wave vector in vacuum 
  PetscScalar SourceTM;    //Source term for TM mode
  PetscScalar SourceTE;  
  PetscScalar DeltaX;      //variable to scale grid size  
  int         m;           //constant used in PML
  PetscScalar R;              
  PetscScalar EpsInVac;    //Variables used when changing power to field
  PetscScalar MuInVac;
  bool        append;
 public:  
  int init_solver(SolveDefine &sv);
  int do_solve(SolveDefine &sv);
  int TM_solution_update(SolveDefine &sv);
  int TE_solution_update(SolveDefine &sv);
  int field_calculation(SolveDefine &sv);
  int generation_cal(SolveDefine &sv);
  int destroy_solver(SolveDefine &sv);
  EM_FEM_Solver();
};

#endif 
