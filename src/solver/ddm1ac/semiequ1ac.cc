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
/*  Last update: Nov 19, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"
#include "vec3.h"
#include <complex>
using namespace std;


void SMCZone::AC1_ddm_inner(int i, PetscScalar omega, PetscScalar *x, Mat *jac,  Vec *r, Mat *jtmp, vector<int> &zofs)
{
  Mat3      A;
  Vec3      b;
  PetscInt  IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i
  PetscScalar area = pcell->area;

  //--------------------------------
  IndexRe[0] = zofs[z]+3*i+0;
  IndexRe[1] = zofs[z]+3*i+1;
  IndexRe[2] = zofs[z]+3*i+2;
  IndexIm[0] = N + zofs[z]+3*i+0;
  IndexIm[1] = N + zofs[z]+3*i+1;
  IndexIm[2] = N + zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    ColRe[0] = zofs[z]+3*nb+0;
    ColRe[1] = zofs[z]+3*nb+1;
    ColRe[2] = zofs[z]+3*nb+2;
    ColIm[0] = N + zofs[z]+3*nb+0;
    ColIm[1] = N + zofs[z]+3*nb+1;
    ColIm[2] = N + zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,IndexRe,3,ColRe,A.m);
    A=A/area;
    MatSetValues(*jac,3,IndexRe,3,ColRe,A.m,INSERT_VALUES);
    MatSetValues(*jac,3,IndexIm,3,ColIm,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,3,IndexRe,3,IndexRe,A.m);
  A=A/area;

  A.m[1] +=  -mt->e;   //dfun(0)/dn(i)
  A.m[2] +=   mt->e;   //dfun(0)/dp(i)

  MatSetValues(*jac,3,IndexRe,3,IndexRe,A.m,INSERT_VALUES);
  MatSetValues(*jac,3,IndexIm,3,IndexIm,A.m,INSERT_VALUES);
  Set_Mat3_zero(A);
  A.m[4] = -omega;
  A.m[8] = -omega;
  MatSetValues(*jac,3,IndexIm,3,IndexRe,A.m,INSERT_VALUES);
  A = -1.0*A;
  MatSetValues(*jac,3,IndexRe,3,IndexIm,A.m,INSERT_VALUES);
  Set_Vec3_zero(b);
  VecSetValues(*r,3,IndexRe,b.v,INSERT_VALUES);
  VecSetValues(*r,3,IndexIm,b.v,INSERT_VALUES);
}


void SMCZone::AC1_ddm_ombc(int i, PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,vector<int> &zofs,DABC &bc)
{
  Mat3      A;
  Vec3      b;
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  PetscInt  IndexRe[3],IndexIm[3];
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  Set_Mat3_I(A);
  IndexRe[0] = zofs[z]+3*i+0;
  IndexRe[1] = zofs[z]+3*i+1;
  IndexRe[2] = zofs[z]+3*i+2;
  IndexIm[0] = N + zofs[z]+3*i+0;
  IndexIm[1] = N + zofs[z]+3*i+1;
  IndexIm[2] = N + zofs[z]+3*i+2;
  MatSetValues(*jac,3,IndexRe,3,IndexRe,A.m,INSERT_VALUES);
  MatSetValues(*jac,3,IndexIm,3,IndexIm,A.m,INSERT_VALUES);
  int om_equ;
  int equ_num = 3;
  int size = pzone->davcell.size();
  for(int j=0;j<electrode.size();j++)
  if(electrode[j]==pcell->bc_index-1)   {om_equ=j;break;}
  MatSetValue(*jac,zofs[z]+3*i,zofs[z]+equ_num*size+om_equ,-1,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[z]+3*i,N+zofs[z]+equ_num*size+om_equ,-1,INSERT_VALUES);
  Set_Vec3_zero(b);
  VecSetValues(*r,3,IndexRe,b.v,INSERT_VALUES);
  VecSetValues(*r,3,IndexIm,b.v,INSERT_VALUES);
}


void SMCZone::AC1_ddm_stkbc(int i, PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,vector<int> &zofs,DABC &bc)
{
  Mat3      A;
  Vec3      b;
  PetscInt  IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  int  z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i
  PetscScalar area = pcell->area;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fp_dpi = 0;
  int stk_equ;
  int equ_num = 3;
  int size = pzone->davcell.size();
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
      {stk_equ=j;break;}
  //--------------------------------
  IndexRe[0] = zofs[z]+3*i+0;
  IndexRe[1] = zofs[z]+3*i+1;
  IndexRe[2] = zofs[z]+3*i+2;
  IndexIm[0] = N + zofs[z]+3*i+0;
  IndexIm[1] = N + zofs[z]+3*i+1;
  IndexIm[2] = N + zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    ColRe[0] = zofs[z]+3*nb+0;
    ColRe[1] = zofs[z]+3*nb+1;
    ColRe[2] = zofs[z]+3*nb+2;
    ColIm[0] = N + zofs[z]+3*nb+0;
    ColIm[1] = N + zofs[z]+3*nb+1;
    ColIm[2] = N + zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,IndexRe,3,ColRe,A.m);
    A=A/area;
    A.m[0]=A.m[1]=A.m[2]=0.0;
    MatSetValues(*jac,3,IndexRe,3,ColRe,A.m,INSERT_VALUES);
    MatSetValues(*jac,3,IndexIm,3,ColIm,A.m,INSERT_VALUES);
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn(ni,fs[i].T,VB)*0.5*pcell->ilen[j];
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp(pi,fs[i].T,VB)*0.5*pcell->ilen[j];
    }
  }

  MatGetValues(*jtmp,3,IndexRe,3,IndexRe,A.m);
  A=A/area;
  A.m[0] = 1.0;
  A.m[1] = 0.0;   //dfun(0)/dn(i)
  A.m[2] = 0.0;   //dfun(0)/dp(i)


  A.m[4] +=   d_Fn_dni/area;   //dfun(1)/dn(i)
  A.m[8] += - d_Fp_dpi/area;   //dfun(2)/dp(i)

  MatSetValues(*jac,3,IndexRe,3,IndexRe,A.m,INSERT_VALUES);
  MatSetValues(*jac,3,IndexIm,3,IndexIm,A.m,INSERT_VALUES);
  Set_Mat3_zero(A);
  A.m[4] = -omega;
  A.m[8] = -omega;
  MatSetValues(*jac,3,IndexIm,3,IndexRe,A.m,INSERT_VALUES);
  A = -1.0*A;
  MatSetValues(*jac,3,IndexRe,3,IndexIm,A.m,INSERT_VALUES);

  MatSetValue(*jac,zofs[z]+3*i,zofs[z]+equ_num*size+stk_equ,-1,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[z]+3*i,N+zofs[z]+equ_num*size+stk_equ,-1,INSERT_VALUES);
  Set_Vec3_zero(b);
  VecSetValues(*r,3,IndexRe,b.v,INSERT_VALUES);
  VecSetValues(*r,3,IndexIm,b.v,INSERT_VALUES);
}


void SMCZone::AC1_ddm_insulator_gate(int i, PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,vector<int> &zofs,DABC &bc)
{
  Mat3    A;
  Vec3    b;
  PetscInt  IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int  equ_num = 3;
  int  size = pzone->davcell.size();
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);

  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_grad_P_dVapp = 0;
  int ins_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
      {ins_equ=j;break;}
  //--------------------------------
  IndexRe[0] = zofs[z]+3*i+0;
  IndexRe[1] = zofs[z]+3*i+1;
  IndexRe[2] = zofs[z]+3*i+2;
  IndexIm[0] = N + zofs[z]+3*i+0;
  IndexIm[1] = N + zofs[z]+3*i+1;
  IndexIm[2] = N + zofs[z]+3*i+2;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    ColRe[0] = zofs[z]+3*nb+0;
    ColRe[1] = zofs[z]+3*nb+1;
    ColRe[2] = zofs[z]+3*nb+2;
    ColIm[0] = N + zofs[z]+3*nb+0;
    ColIm[1] = N + zofs[z]+3*nb+1;
    ColIm[2] = N + zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,IndexRe,3,ColRe,A.m);
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar s=eps_ox/Thick;
      d_grad_P_dVi += -0.5*pcell->ilen[j]*0.25*s*3;
      d_grad_P_dVapp += 0.5*pcell->ilen[j]*s/area;
      A.m[0]+= -0.5*pcell->ilen[j]*0.25*s/area;
    }
    MatSetValues(*jac,3,IndexRe,3,ColRe,A.m,INSERT_VALUES);
    MatSetValues(*jac,3,IndexIm,3,ColIm,A.m,INSERT_VALUES);
  }
  mt->mapping(&pzone->danode[i],&aux[i],0);

  MatGetValues(*jtmp,3,IndexRe,3,IndexRe,A.m);
  A=A/area;
  A.m[0] +=  d_grad_P_dVi/area;             //dfun(0)/dP(i)
  A.m[1] +=  -mt->e;   //dfun(0)/dn(i)
  A.m[2] +=   mt->e;   //dfun(0)/dp(i)

  MatSetValues(*jac,3,IndexRe,3,IndexRe,A.m,INSERT_VALUES);
  MatSetValues(*jac,3,IndexIm,3,IndexIm,A.m,INSERT_VALUES);

  Set_Mat3_zero(A);
  A.m[4] = -omega;
  A.m[8] = -omega;
  MatSetValues(*jac,3,IndexIm,3,IndexRe,A.m,INSERT_VALUES);
  A = -1.0*A;
  MatSetValues(*jac,3,IndexRe,3,IndexIm,A.m,INSERT_VALUES);

  MatSetValue(*jac,zofs[z]+3*i,zofs[z]+equ_num*size+ins_equ,d_grad_P_dVapp,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[z]+3*i,N+zofs[z]+equ_num*size+ins_equ,d_grad_P_dVapp,INSERT_VALUES);

  Set_Vec3_zero(b);
  VecSetValues(*r,3,IndexRe,b.v,INSERT_VALUES);
  VecSetValues(*r,3,IndexIm,b.v,INSERT_VALUES);
}


void SMCZone::AC1_ddm_interface(int i, PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,
       vector<int> &zofs,DABC &bc,ISZone *pz, int n)
{
  Mat3    A;
  Vec3    b;
  PetscInt  IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int  equ_num = 3;
  int  size = pzone->davcell.size();
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);

  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;

  //--------------------------------
  IndexRe[0] = zofs[z]+3*i+0;
  IndexRe[1] = zofs[z]+3*i+1;
  IndexRe[2] = zofs[z]+3*i+2;
  IndexIm[0] = N + zofs[z]+3*i+0;
  IndexIm[1] = N + zofs[z]+3*i+1;
  IndexIm[2] = N + zofs[z]+3*i+2;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    ColRe[0] = zofs[z]+3*nb+0;
    ColRe[1] = zofs[z]+3*nb+1;
    ColRe[2] = zofs[z]+3*nb+2;
    ColIm[0] = N + zofs[z]+3*nb+0;
    ColIm[1] = N + zofs[z]+3*nb+1;
    ColIm[2] = N + zofs[z]+3*nb+2;
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    MatGetValues(*jtmp,3,IndexRe,3,ColRe,A.m);
    A=A/area;
    A.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    A.m[1] =  0;
    A.m[2] =  0;
    MatSetValues(*jac,3,IndexRe,3,ColRe,A.m,INSERT_VALUES);
    MatSetValues(*jac,3,IndexIm,3,ColIm,A.m,INSERT_VALUES);
  }
  mt->mapping(&pzone->danode[i],&aux[i],0);

  MatGetValues(*jtmp,3,IndexRe,3,IndexRe,A.m);
  A=A/area;
  A.m[0] =  d_grad_P_dVi;  //dfun(0)/dP(i)
  A.m[1] =  -mt->e*area;   //dfun(0)/dn(i)
  A.m[2] =   mt->e*area;   //dfun(0)/dp(i)
  MatSetValues(*jac,3,IndexRe,3,IndexRe,A.m,INSERT_VALUES);
  MatSetValues(*jac,3,IndexIm,3,IndexIm,A.m,INSERT_VALUES);

  Set_Mat3_zero(A);
  A.m[4] = -omega;
  A.m[8] = -omega;
  MatSetValues(*jac,3,IndexIm,3,IndexRe,A.m,INSERT_VALUES);
  A = -1.0*A;
  MatSetValues(*jac,3,IndexRe,3,IndexIm,A.m,INSERT_VALUES);

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  int column = zofs[pz->pzone->zone_index]+n;
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar Vr_n = x[zofs[pz->pzone->zone_index]+nb];     //potential of nb node
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    PetscScalar d_grad_P_dVj = pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+nb,d_grad_P_dVj,INSERT_VALUES);
    MatSetValue(*jac,N+zofs[z]+3*i+0,N+zofs[pz->pzone->zone_index]+nb,d_grad_P_dVj,INSERT_VALUES);
  }
  MatSetValue(*jac,zofs[z]+3*i+0,zofs[z]+3*i+0,d_grad_P_dVi,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[z]+3*i+0,N+zofs[z]+3*i+0,d_grad_P_dVi,INSERT_VALUES);

  Set_Vec3_zero(b);
  VecSetValues(*r,3,IndexRe,b.v,INSERT_VALUES);
  VecSetValues(*r,3,IndexIm,b.v,INSERT_VALUES);
}


void SMCZone::AC1_ddm_homojunction(int i, PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,
       vector<int> &zofs,DABC &bc,SMCZone *pz, int n)
{
  Mat3           A1,A2,AJ;
  Vec3    b;
  PetscInt       IndexRe[3],IndexIm[3],ColRe[3],ColIm[3],NbRe[3],NbIm[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  PetscScalar e  =  mt->e;
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar d_grad_P_dVi = 0;
  //--------------------------------
  IndexRe[0] = zofs[z]+3*i+0;
  IndexRe[1] = zofs[z]+3*i+1;
  IndexRe[2] = zofs[z]+3*i+2;
  IndexIm[0] = N + zofs[z]+3*i+0;
  IndexIm[1] = N + zofs[z]+3*i+1;
  IndexIm[2] = N + zofs[z]+3*i+2;
  //--------------------------------

  if(zone_index > pz->pzone->zone_index)
  {
        Set_Mat3_I(AJ);
        MatSetValues(*jac,3,IndexRe,3,IndexRe,AJ.m,INSERT_VALUES);
        MatSetValues(*jac,3,IndexIm,3,IndexIm,AJ.m,INSERT_VALUES);
        AJ = -1.0*AJ;
        NbRe[0] = zofs[pz->zone_index]+3*n+0;
        NbRe[1] = zofs[pz->zone_index]+3*n+1;
        NbRe[2] = zofs[pz->zone_index]+3*n+2;
        NbIm[0] = N+zofs[pz->zone_index]+3*n+0;
        NbIm[1] = N+zofs[pz->zone_index]+3*n+1;
        NbIm[2] = N+zofs[pz->zone_index]+3*n+2;
        MatSetValues(*jac,3,IndexRe,3,NbRe,AJ.m,INSERT_VALUES);
        MatSetValues(*jac,3,IndexIm,3,NbIm,AJ.m,INSERT_VALUES);
        return;
  }

  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[z]+3*nb+0];     //potential of nb node
    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    ColRe[0] = zofs[z]+3*nb+0;
    ColRe[1] = zofs[z]+3*nb+1;
    ColRe[2] = zofs[z]+3*nb+2;
    ColIm[0] = N+zofs[z]+3*nb+0;
    ColIm[1] = N+zofs[z]+3*nb+1;
    ColIm[2] = N+zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,IndexRe,3,ColRe,A1.m);
    //-------------------------------------
    //df(x)i/dx(r)
    A1.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];   //dfun(0)/dP(r)
    A1.m[1] =  0;                                       //dfun(0)/dn(r)
    A1.m[2] =  0;                                       //dfun(0)/dp(r)
    MatSetValues(*jac,3,IndexRe,3,ColRe,A1.m,INSERT_VALUES);
    MatSetValues(*jac,3,IndexIm,3,ColIm,A1.m,INSERT_VALUES);
  }

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  NbRe[0] = zofs[pz->zone_index]+3*n+0;
  NbRe[1] = zofs[pz->zone_index]+3*n+1;
  NbRe[2] = zofs[pz->zone_index]+3*n+2;
  NbIm[0] = N+zofs[pz->zone_index]+3*n+0;
  NbIm[1] = N+zofs[pz->zone_index]+3*n+1;
  NbIm[2] = N+zofs[pz->zone_index]+3*n+2;
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar Vj = x[zofs[pz->zone_index]+3*nb+0];     //potential of nb node
    //-------------------------------------
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    ColRe[0] = zofs[pz->zone_index]+3*nb+0;
    ColRe[1] = zofs[pz->zone_index]+3*nb+1;
    ColRe[2] = zofs[pz->zone_index]+3*nb+2;
    ColIm[0] = N+zofs[pz->zone_index]+3*nb+0;
    ColIm[1] = N+zofs[pz->zone_index]+3*nb+1;
    ColIm[2] = N+zofs[pz->zone_index]+3*nb+2;
    MatGetValues(*jtmp,3,NbRe,3,ColRe,A2.m);
    //-------------------------------------
    //df(x)i/dx(r)
    A2.m[0] =  pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];     //dfun(0)/dP(r)
    A2.m[1] =  0;                                                     //dfun(0)/dn(r)
    A2.m[2] =  0;                                                     //dfun(0)/dp(r)
    MatSetValues(*jac,3,IndexRe,3,ColRe,A2.m,INSERT_VALUES);
    MatSetValues(*jac,3,IndexIm,3,ColIm,A2.m,INSERT_VALUES);
  }

  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  PetscScalar area = pcell->area + ncell->area;
  MatGetValues(*jtmp,3,IndexRe,3,IndexRe,A1.m);
  MatGetValues(*jtmp,3,NbRe,3,NbRe,A2.m);
  AJ=A1+A2;

  AJ.m[0] =  d_grad_P_dVi;    //dfun(0)/dP(i)
  AJ.m[1] =  -e*area;  //dfun(0)/dn(i)
  AJ.m[2] =   e*area;  //dfun(0)/dp(i)

  MatSetValues(*jac,3,IndexRe,3,IndexRe,AJ.m,INSERT_VALUES);
  MatSetValues(*jac,3,IndexIm,3,IndexIm,AJ.m,INSERT_VALUES);
  Set_Mat3_zero(AJ);
  AJ.m[4] = -omega*area;
  AJ.m[8] = -omega*area;
  MatSetValues(*jac,3,IndexIm,3,IndexRe,AJ.m,INSERT_VALUES);
  AJ = -1.0*AJ;
  MatSetValues(*jac,3,IndexRe,3,IndexIm,AJ.m,INSERT_VALUES);
  Set_Vec3_zero(b);
  VecSetValues(*r,3,IndexRe,b.v,INSERT_VALUES);
  VecSetValues(*r,3,IndexIm,b.v,INSERT_VALUES);
}


void SMCZone::AC1_ddm_heterojunction(int i, PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,
       vector<int> &zofs,DABC &bc,SMCZone *pz, int n)
{
  Mat3    A;
  Vec3    b;
  PetscInt       IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  int z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  PetscScalar k = mt->kb;
  PetscScalar e  =  mt->e;
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i
  PetscScalar Ti = fs[i].T;
  PetscScalar Eci =  -e*(Vi + aux[i].affinity + mt->kb*fs[i].T/e*log(aux[i].Nc));
  PetscScalar Evi =  -e*(Vi + aux[i].affinity + aux[i].Eg/e - mt->kb*fs[i].T/e*log(aux[i].Nv));
  PetscScalar area = pcell->area;
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fp_dpi = 0;
  //--------------------------------
  IndexRe[0] = zofs[z]+3*i+0;
  IndexRe[1] = zofs[z]+3*i+1;
  IndexRe[2] = zofs[z]+3*i+2;
  IndexIm[0] = N + zofs[z]+3*i+0;
  IndexIm[1] = N + zofs[z]+3*i+1;
  IndexIm[2] = N + zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];

    PetscScalar Vj = x[zofs[z]+3*nb+0];     //potential of nb node
    PetscScalar nj = x[zofs[z]+3*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+3*nb+2];     //hole density of nb node
    //-------------------------------------
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];

    ColRe[0] = zofs[z]+3*nb+0;
    ColRe[1] = zofs[z]+3*nb+1;
    ColRe[2] = zofs[z]+3*nb+2;
    ColIm[0] = N + zofs[z]+3*nb+0;
    ColIm[1] = N + zofs[z]+3*nb+1;
    ColIm[2] = N + zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,IndexRe,3,ColRe,A.m);
    A=A/area;
    //-------------------------------------
    //df(x)i/dx(r)
    if(pzone->zone_index < pz->pzone->zone_index)
    {
        A.m[0] =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];//dfun(0)/dP(r)
        A.m[1] =  0;                                     //dfun(0)/dn(r)
        A.m[2] =  0;                                     //dfun(0)/dp(r)
    }
    else
    {
        A.m[0] =  0;                                     //dfun(0)/dP(r)
        A.m[1] =  0;                                     //dfun(0)/dn(r)
        A.m[2] =  0;                                     //dfun(0)/dp(r)
    }

    //Thermal emit current
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Ecj = -e*(Vi + pz->aux[n].affinity + k*pz->fs[n].T/e*log(pz->aux[n].Nc));
      PetscScalar Evj = -e*(Vi + pz->aux[n].affinity - k*pz->fs[n].T/e*log(pz->aux[n].Nv) + pz->aux[n].Eg/e);
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
      if(Ecj > Eci)
      {
        PetscScalar pm = pz->mt->band->EffecElecMass(Ti)/mt->band->EffecElecMass(Ti);
        d_Fn_dni += 2*e*pm*mt->band->ThermalVn(Ti)*exp(-(Ecj-Eci)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fn_dnd = -2*e*pz->mt->band->ThermalVn(Ti)*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+1,zofs[pz->pzone->zone_index]+3*n+1,d_Fn_dnd,ADD_VALUES);
        MatSetValue(*jac,N+zofs[z]+3*i+1,N+zofs[pz->pzone->zone_index]+3*n+1,d_Fn_dnd,ADD_VALUES);
      }
      else
      {
        PetscScalar pm = mt->band->EffecElecMass(Ti)/pz->mt->band->EffecElecMass(Ti);
        d_Fn_dni += -2*e*mt->band->ThermalVn(Ti)*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fn_dnd = 2*e*pm*pz->mt->band->ThermalVn(Ti)*exp(-(Eci-Ecj)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+1,zofs[pz->pzone->zone_index]+3*n+1,d_Fn_dnd,ADD_VALUES);
        MatSetValue(*jac,N+zofs[z]+3*i+1,N+zofs[pz->pzone->zone_index]+3*n+1,d_Fn_dnd,ADD_VALUES);
      }

      if(Evj > Evi)
      {
        PetscScalar pm = pz->mt->band->EffecHoleMass(Ti)/mt->band->EffecHoleMass(Ti);
        d_Fp_dpi += -2*e*pm*mt->band->ThermalVp(Ti)*exp(-(Evj-Evi)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fp_dpd = 2*e*pz->mt->band->ThermalVp(Ti)*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+2,zofs[pz->pzone->zone_index]+3*n+2,-d_Fp_dpd,ADD_VALUES);
        MatSetValue(*jac,N+zofs[z]+3*i+2,N+zofs[pz->pzone->zone_index]+3*n+2,-d_Fp_dpd,ADD_VALUES);
      }
      else
      {
        PetscScalar pm = mt->band->EffecHoleMass(Ti)/pz->mt->band->EffecHoleMass(Ti);
        d_Fp_dpi += 2*e*mt->band->ThermalVp(Ti)*0.5*pcell->ilen[j]/pcell->area;
        PetscScalar d_Fp_dpd = -2*e*pm*pz->mt->band->ThermalVp(Ti)*exp(-(Evi-Evj)/(k*Ti))*0.5*pcell->ilen[j]/pcell->area;
        MatSetValue(*jac,zofs[z]+3*i+2,zofs[pz->pzone->zone_index]+3*n+2,-d_Fp_dpd,ADD_VALUES);
        MatSetValue(*jac,N+zofs[z]+3*i+2,N+zofs[pz->pzone->zone_index]+3*n+2,-d_Fp_dpd,ADD_VALUES);
      }
      MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
      MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
    }
    MatSetValues(*jac,3,IndexRe,3,ColRe,A.m,INSERT_VALUES);
    MatSetValues(*jac,3,IndexIm,3,ColIm,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,3,IndexRe,3,IndexRe,A.m);
  A=A/area;
  if(pzone->zone_index < pz->pzone->zone_index)
  {
    A.m[0] =  d_grad_P_dVi;    //dfun(0)/dP(i)
    A.m[1] =  -e*pcell->area;  //dfun(0)/dn(i)
    A.m[2] =   e*pcell->area;  //dfun(0)/dp(i)
    for(int j=0;j<ncell->nb_num;j++)
    {
      int    nb = ncell->nb_array[j];
      A.m[0] += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
      PetscScalar value = pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
      MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*nb+0,value,INSERT_VALUES);
      MatSetValue(*jac,N+zofs[z]+3*i+0,N+zofs[pz->pzone->zone_index]+3*nb+0,value,INSERT_VALUES);
    }
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*n+1,-e*ncell->area,INSERT_VALUES);
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*n+2,e*ncell->area,INSERT_VALUES);
    MatSetValue(*jac,N+zofs[z]+3*i+0,N+zofs[pz->pzone->zone_index]+3*n+1,-e*ncell->area,INSERT_VALUES);
    MatSetValue(*jac,N+zofs[z]+3*i+0,N+zofs[pz->pzone->zone_index]+3*n+2,e*ncell->area,INSERT_VALUES);
  }
  else
  {
    A.m[0] =   1;  //dfun(0)/dP(i)
    A.m[1] =   0;  //dfun(0)/dn(i)
    A.m[2] =   0;  //dfun(0)/dp(i)
    MatSetValue(*jac,zofs[z]+3*i+0,zofs[pz->pzone->zone_index]+3*n+0,-1,INSERT_VALUES);
    MatSetValue(*jac,N+zofs[z]+3*i+0,N+zofs[pz->pzone->zone_index]+3*n+0,-1,INSERT_VALUES);
  }

  A.m[4] +=   d_Fn_dni;          //dfun(1)/dp(i)
  A.m[8] +=  -d_Fp_dpi;          //dfun(2)/dp(i)

  MatSetValues(*jac,3,IndexRe,3,IndexRe,A.m,INSERT_VALUES);
  MatSetValues(*jac,3,IndexIm,3,IndexIm,A.m,INSERT_VALUES);

  Set_Mat3_zero(A);
  A.m[4] = -omega;
  A.m[8] = -omega;
  MatSetValues(*jac,3,IndexIm,3,IndexRe,A.m,INSERT_VALUES);
  A = -1.0*A;
  MatSetValues(*jac,3,IndexRe,3,IndexIm,A.m,INSERT_VALUES);
  Set_Vec3_zero(b);
  VecSetValues(*r,3,IndexRe,b.v,INSERT_VALUES);
  VecSetValues(*r,3,IndexIm,b.v,INSERT_VALUES);
}


void SMCZone::AC1_om_electrode(int i,PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,
       vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{

  int equ_num = 3;
  int size = pzone->davcell.size();
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  Vec3    J1,J2;
  complex <PetscScalar> A[3];
  PetscInt       IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;
  complex <PetscScalar> Z1(R,omega*L);
  complex <PetscScalar> Y2(0,omega*C);
  PetscInt EquRe = zofs[zone_index]+equ_num*size+i;
  PetscInt EquIm = N + zofs[zone_index]+equ_num*size+i;
  complex <PetscScalar> I(0.0,1.0);

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    complex <PetscScalar> dJdisp_dVi=0,dJdisp_dVr=0;
    IndexRe[0] = zofs[z]+3*node+0;
    IndexRe[1] = zofs[z]+3*node+1;
    IndexRe[2] = zofs[z]+3*node+2;
    IndexIm[0] = N + zofs[z]+3*node+0;
    IndexIm[1] = N + zofs[z]+3*node+1;
    IndexIm[2] = N + zofs[z]+3*node+2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      ColRe[0] = zofs[z]+3*nb+0;
      ColRe[1] = zofs[z]+3*nb+1;
      ColRe[2] = zofs[z]+3*nb+2;
      ColIm[0] = N + zofs[z]+3*nb+0;
      ColIm[1] = N + zofs[z]+3*nb+1;
      ColIm[2] = N + zofs[z]+3*nb+2;

      MatGetValues(*jtmp,1,&IndexRe[1],3,ColRe,J1.v);
      MatGetValues(*jtmp,1,&IndexRe[2],3,ColRe,J2.v);
      //for displacement current
      dJdisp_dVi += omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;
      dJdisp_dVr = -omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;
      A[0]=(J1.v[0]-J2.v[0]+dJdisp_dVr)*DeviceDepth*Z1;
      A[1]=(J1.v[1]-J2.v[1])*DeviceDepth*Z1;
      A[2]=(J1.v[2]-J2.v[2])*DeviceDepth*Z1;
      J1.v[0] = A[0].real();
      J1.v[1] = A[1].real();
      J1.v[2] = A[2].real();
      J2.v[0] = A[0].imag();
      J2.v[1] = A[1].imag();
      J2.v[2] = A[2].imag();
      MatSetValues(*jac,1,&EquRe,3,ColRe,J1.v,ADD_VALUES);
      MatSetValues(*jac,1,&EquRe,3,ColIm,(-1.0*J2).v,ADD_VALUES);
      MatSetValues(*jac,1,&EquIm,3,ColRe,J2.v,ADD_VALUES);
      MatSetValues(*jac,1,&EquIm,3,ColIm,J1.v,ADD_VALUES);
    }

    MatGetValues(*jtmp,1,&IndexRe[1],3,IndexRe,J1.v);
    MatGetValues(*jtmp,1,&IndexRe[2],3,IndexRe,J2.v);
    A[0]=(J1.v[0]-J2.v[0]+dJdisp_dVi)*DeviceDepth*Z1;
    A[1]=(J1.v[1]-J2.v[1])*DeviceDepth*Z1;
    A[2]=(J1.v[2]-J2.v[2])*DeviceDepth*Z1;
    J1.v[0] = A[0].real();
    J1.v[1] = A[1].real();
    J1.v[2] = A[2].real();
    J2.v[0] = A[0].imag();
    J2.v[1] = A[1].imag();
    J2.v[2] = A[2].imag();
    MatSetValues(*jac,1,&EquRe,3,IndexRe,J1.v,ADD_VALUES);
    MatSetValues(*jac,1,&EquRe,3,IndexIm,(-1.0*J2).v,ADD_VALUES);
    MatSetValues(*jac,1,&EquIm,3,IndexRe,J2.v,ADD_VALUES);
    MatSetValues(*jac,1,&EquIm,3,IndexIm,J1.v,ADD_VALUES);
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

  complex <PetscScalar> K = Z1*Y2 + PetscScalar(1.0);
  MatSetValue(*jac,EquRe,EquRe,K.real(),INSERT_VALUES);
  MatSetValue(*jac,EquRe,EquIm,-K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquRe,K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquIm,K.real(),INSERT_VALUES);

  VecSetValue(*r,EquRe,pbc->Vac, INSERT_VALUES);
  VecSetValue(*r,EquIm,0, INSERT_VALUES);

}


void SMCZone::AC1_stk_electrode(int i,PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,
       vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;
  complex <PetscScalar> Z1(R,omega*L);
  complex <PetscScalar> Y2(0,omega*C);
  complex <PetscScalar> I(0,1);
  Vec3    J1,J2;
  complex <PetscScalar> A[3];
  PetscInt IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  PetscInt EquRe = zofs[zone_index]+equ_num*size+i;
  PetscInt EquIm = N + zofs[zone_index]+equ_num*size+i;

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar VB = pbc->WorkFunction-aux[node].affinity;
    PetscScalar ni = x[zofs[zone_index]+3*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+3*node+2];     //hole density of node i
    PetscScalar dJ_dni = -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB)*0.5*pcell->ilen[0]
                  -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB)*0.5*pcell->ilen[pcell->nb_num-1];
    PetscScalar dJ_dpi = -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB)*0.5*pcell->ilen[0]
                  -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB)*0.5*pcell->ilen[pcell->nb_num-1];
    complex <PetscScalar> dJdisp_dVi=0,dJdisp_dVr=0;
    IndexRe[0] = zofs[z]+3*node+0;
    IndexRe[1] = zofs[z]+3*node+1;
    IndexRe[2] = zofs[z]+3*node+2;
    IndexIm[0] = N + zofs[z]+3*node+0;
    IndexIm[1] = N + zofs[z]+3*node+1;
    IndexIm[2] = N + zofs[z]+3*node+2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      ColRe[0] = zofs[z]+3*nb+0;
      ColRe[1] = zofs[z]+3*nb+1;
      ColRe[2] = zofs[z]+3*nb+2;
      ColIm[0] = N + zofs[z]+3*nb+0;
      ColIm[1] = N + zofs[z]+3*nb+1;
      ColIm[2] = N + zofs[z]+3*nb+2;

      //for displacement current
      dJdisp_dVi += omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;
      dJdisp_dVr = -omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;

      A[0]=dJdisp_dVr*DeviceDepth*Z1;
      A[1]=0*DeviceDepth*Z1;
      A[2]=0*DeviceDepth*Z1;
      J1.v[0] = A[0].real();
      J1.v[1] = A[1].real();
      J1.v[2] = A[2].real();
      J2.v[0] = A[0].imag();
      J2.v[1] = A[1].imag();
      J2.v[2] = A[2].imag();
      MatSetValues(*jac,1,&EquRe,3,ColRe,J1.v,ADD_VALUES);
      MatSetValues(*jac,1,&EquRe,3,ColIm,(-1.0*J2).v,ADD_VALUES);
      MatSetValues(*jac,1,&EquIm,3,ColRe,J2.v,ADD_VALUES);
      MatSetValues(*jac,1,&EquIm,3,ColIm,J1.v,ADD_VALUES);
    }
    A[0] = dJdisp_dVi*DeviceDepth*Z1;
    A[1] = dJ_dni*DeviceDepth*Z1;
    A[2] = dJ_dpi*DeviceDepth*Z1;
    J1.v[0] = A[0].real();
    J1.v[1] = A[1].real();
    J1.v[2] = A[2].real();
    J2.v[0] = A[0].imag();
    J2.v[1] = A[1].imag();
    J2.v[2] = A[2].imag();
    MatSetValues(*jac,1,&EquRe,3,IndexRe,J1.v,ADD_VALUES);
    MatSetValues(*jac,1,&EquRe,3,IndexIm,(-1.0*J2).v,ADD_VALUES);
    MatSetValues(*jac,1,&EquIm,3,IndexRe,J2.v,ADD_VALUES);
    MatSetValues(*jac,1,&EquIm,3,IndexIm,J1.v,ADD_VALUES);
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

  complex <PetscScalar> K = Z1*Y2 + PetscScalar(1.0);
  MatSetValue(*jac,EquRe,EquRe,K.real(),INSERT_VALUES);
  MatSetValue(*jac,EquRe,EquIm,-K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquRe,K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquIm,K.real(),INSERT_VALUES);

  VecSetValue(*r,EquRe,pbc->Vac, INSERT_VALUES);
  VecSetValue(*r,EquIm,0, INSERT_VALUES);
}


void SMCZone::AC1_ins_electrode(int i,PetscScalar omega, PetscScalar *x,Mat *jac,Vec *r,Mat *jtmp,
       vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;

  complex <PetscScalar> Z1(R,omega*L);
  complex <PetscScalar> Y2(0,omega*C);
  complex <PetscScalar> I(0,1);
  Vec3    J1,J2;
  complex <PetscScalar> A[3];
  PetscInt IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  PetscInt EquRe = zofs[zone_index]+equ_num*size+i;
  PetscInt EquIm = N + zofs[zone_index]+equ_num*size+i;

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    complex <PetscScalar> dJdisp_dVi=0,dJdisp_dVr=0;
    IndexRe[0] = zofs[z]+3*node+0;
    IndexRe[1] = zofs[z]+3*node+1;
    IndexRe[2] = zofs[z]+3*node+2;
    IndexIm[0] = N + zofs[z]+3*node+0;
    IndexIm[1] = N + zofs[z]+3*node+1;
    IndexIm[2] = N + zofs[z]+3*node+2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      ColRe[0] = zofs[z]+3*nb+0;
      ColRe[1] = zofs[z]+3*nb+1;
      ColRe[2] = zofs[z]+3*nb+2;
      ColIm[0] = N + zofs[z]+3*nb+0;
      ColIm[1] = N + zofs[z]+3*nb+1;
      ColIm[2] = N + zofs[z]+3*nb+2;

      //for displacement current
      dJdisp_dVi += omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;
      dJdisp_dVr = -omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;

      A[0]=dJdisp_dVr*DeviceDepth*Z1;
      A[1]=0.0;
      A[2]=0.0;
      J1.v[0] = A[0].real();
      J1.v[1] = A[1].real();
      J1.v[2] = A[2].real();
      J2.v[0] = A[0].imag();
      J2.v[1] = A[1].imag();
      J2.v[2] = A[2].imag();
      MatSetValues(*jac,1,&EquRe,3,ColRe,J1.v,ADD_VALUES);
      MatSetValues(*jac,1,&EquRe,3,ColIm,(-1.0*J2).v,ADD_VALUES);
      MatSetValues(*jac,1,&EquIm,3,ColRe,J2.v,ADD_VALUES);
      MatSetValues(*jac,1,&EquIm,3,ColIm,J1.v,ADD_VALUES);
    }
    A[0] = dJdisp_dVi*DeviceDepth*Z1;
    A[1] = 0.0;
    A[2] = 0.0;
    J1.v[0] = A[0].real();
    J1.v[1] = A[1].real();
    J1.v[2] = A[2].real();
    J2.v[0] = A[0].imag();
    J2.v[1] = A[1].imag();
    J2.v[2] = A[2].imag();
    MatSetValues(*jac,1,&EquRe,3,IndexRe,J1.v,ADD_VALUES);
    MatSetValues(*jac,1,&EquRe,3,IndexIm,(-1.0*J2).v,ADD_VALUES);
    MatSetValues(*jac,1,&EquIm,3,IndexRe,J2.v,ADD_VALUES);
    MatSetValues(*jac,1,&EquIm,3,IndexIm,J1.v,ADD_VALUES);
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

  complex <PetscScalar> K = Z1*Y2 + PetscScalar(1.0);
  MatSetValue(*jac,EquRe,EquRe,K.real(),INSERT_VALUES);
  MatSetValue(*jac,EquRe,EquIm,-K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquRe,K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquIm,K.real(),INSERT_VALUES);

  VecSetValue(*r,EquRe,pbc->Vac, INSERT_VALUES);
  VecSetValue(*r,EquIm,0, INSERT_VALUES);
}


void SMCZone::AC1_om_electrode_current(int i,PetscScalar omega, PetscScalar *x,Vec *v,Mat *jtmp,
       vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth, PetscScalar &IacRe, PetscScalar &IacIm)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  Mat3    A;
  Vec3    V1,V2,b1,b2;
  PetscInt       IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  IacRe = 0;
  IacIm = 0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);

    IndexRe[0] = zofs[z]+3*node+0;
    IndexRe[1] = zofs[z]+3*node+1;
    IndexRe[2] = zofs[z]+3*node+2;
    IndexIm[0] = N + zofs[z]+3*node+0;
    IndexIm[1] = N + zofs[z]+3*node+1;
    IndexIm[2] = N + zofs[z]+3*node+2;
    PetscScalar ViRe,ViIm;
    VecGetValues(*v,1,&IndexRe[0],&ViRe);
    VecGetValues(*v,1,&IndexIm[0],&ViIm);

    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      ColRe[0] = zofs[z]+3*nb+0;
      ColRe[1] = zofs[z]+3*nb+1;
      ColRe[2] = zofs[z]+3*nb+2;
      ColIm[0] = N + zofs[z]+3*nb+0;
      ColIm[1] = N + zofs[z]+3*nb+1;
      ColIm[2] = N + zofs[z]+3*nb+2;
      MatGetValues(*jtmp,3,IndexRe,3,ColRe,A.m);
      VecGetValues(*v,3,ColRe,V1.v);
      VecGetValues(*v,3,ColIm,V2.v);
      b1 = A*V1;
      b2 = A*V2;
      IacRe += (b1.v[1]-b1.v[2])*DeviceDepth;
      IacIm += (b2.v[1]-b2.v[2])*DeviceDepth;
      //for displacement current
      PetscScalar VjRe = V1.v[0];
      PetscScalar VjIm = V2.v[0];
      IacRe += -omega*aux[node].eps*(ViIm-VjIm)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
      IacIm +=  omega*aux[node].eps*(ViRe-VjRe)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
    }
    MatGetValues(*jtmp,3,IndexRe,3,IndexRe,A.m);
    VecGetValues(*v,3,IndexRe,V1.v);
    VecGetValues(*v,3,IndexIm,V2.v);
    b1 = A*V1;
    b2 = A*V2;
    IacRe += (b1.v[1]-b1.v[2])*DeviceDepth;
    IacIm += (b2.v[1]-b2.v[2])*DeviceDepth;
  }
}


void SMCZone::AC1_stk_electrode_current(int i,PetscScalar omega, PetscScalar *x, Vec *v,Mat *jtmp,
       vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth, PetscScalar &IacRe, PetscScalar &IacIm)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));
  PetscInt       IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  IacRe = 0;
  IacIm = 0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar VB = pbc->WorkFunction-aux[node].affinity;
    PetscScalar ni = x[zofs[zone_index]+3*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+3*node+2];     //hole density of node i
    PetscScalar dJ_dni = -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB)*0.5*pcell->ilen[0]
                  -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB)*0.5*pcell->ilen[pcell->nb_num-1];
    PetscScalar dJ_dpi = -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB)*0.5*pcell->ilen[0]
                  -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB)*0.5*pcell->ilen[pcell->nb_num-1];

    IndexRe[0] = zofs[z]+3*node+0;
    IndexRe[1] = zofs[z]+3*node+1;
    IndexRe[2] = zofs[z]+3*node+2;
    IndexIm[0] = N + zofs[z]+3*node+0;
    IndexIm[1] = N + zofs[z]+3*node+1;
    IndexIm[2] = N + zofs[z]+3*node+2;
    PetscScalar ViRe,ViIm;
    PetscScalar niRe,niIm;
    PetscScalar piRe,piIm;
    VecGetValues(*v,1,&IndexRe[0],&ViRe);
    VecGetValues(*v,1,&IndexIm[0],&ViIm);
    VecGetValues(*v,1,&IndexRe[1],&niRe);
    VecGetValues(*v,1,&IndexIm[1],&niIm);
    VecGetValues(*v,1,&IndexRe[2],&piRe);
    VecGetValues(*v,1,&IndexIm[2],&piIm);
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      ColRe[0] = zofs[z]+3*nb+0;
      ColRe[1] = zofs[z]+3*nb+1;
      ColRe[2] = zofs[z]+3*nb+2;
      ColIm[0] = N + zofs[z]+3*nb+0;
      ColIm[1] = N + zofs[z]+3*nb+1;
      ColIm[2] = N + zofs[z]+3*nb+2;

      //for displacement current
      PetscScalar VjRe,VjIm;
      VecGetValues(*v,1,&ColRe[0],&VjRe);
      VecGetValues(*v,1,&ColIm[0],&VjIm);
      IacRe += -omega*aux[node].eps*(ViIm-VjIm)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
      IacIm +=  omega*aux[node].eps*(ViRe-VjRe)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
    }
    IacRe += (dJ_dni*niRe + dJ_dpi*piRe)*DeviceDepth;
    IacIm += (dJ_dni*niIm + dJ_dpi*piIm)*DeviceDepth;
  }
}

void SMCZone::AC1_ins_electrode_current(int i,PetscScalar omega, PetscScalar *x, Vec *v,Mat *jtmp,
       vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth, PetscScalar &IacRe, PetscScalar &IacIm)
{
  int equ_num = 3;
  int size = pzone->davcell.size();
  int   z = zone_index;
  int   N;
  MatGetSize(*jtmp,&N,PETSC_NULL);
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));
  PetscInt       IndexRe[3],IndexIm[3],ColRe[3],ColIm[3];
  IacRe = 0;
  IacIm = 0;
  //for displacement current
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    IndexRe[0] = zofs[z]+3*node+0;
    IndexRe[1] = zofs[z]+3*node+1;
    IndexRe[2] = zofs[z]+3*node+2;
    IndexIm[0] = N + zofs[z]+3*node+0;
    IndexIm[1] = N + zofs[z]+3*node+1;
    IndexIm[2] = N + zofs[z]+3*node+2;
    PetscScalar ViRe,ViIm;
    VecGetValues(*v,1,&IndexRe[0],&ViRe);
    VecGetValues(*v,1,&IndexIm[0],&ViIm);
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      ColRe[0] = zofs[z]+3*nb+0;
      ColRe[1] = zofs[z]+3*nb+1;
      ColRe[2] = zofs[z]+3*nb+2;
      ColIm[0] = N + zofs[z]+3*nb+0;
      ColIm[1] = N + zofs[z]+3*nb+1;
      ColIm[2] = N + zofs[z]+3*nb+2;

      PetscScalar VjRe,VjIm;
      VecGetValues(*v,1,&ColRe[0],&VjRe);
      VecGetValues(*v,1,&ColIm[0],&VjIm);
      IacRe += -omega*aux[node].eps*(ViIm-VjIm)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
      IacIm +=  omega*aux[node].eps*(ViRe-VjRe)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
    }
  }
}
