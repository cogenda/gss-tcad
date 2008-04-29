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

void ElZone::AC1_ddm_inner(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int  N = zofs[zofs.size()-1];
  PetscScalar d_grad_P_dVi = 0;
  int IndexRe = zofs[zone_index]+i;
  int IndexIm = N+zofs[zone_index]+i;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    int ColRe = zofs[zone_index]+nb;
    int ColIm = N+zofs[zone_index]+nb;
    PetscScalar d_grad_P_dVj =  pcell->elen[j]/pcell->ilen[j]/pcell->area;
    MatSetValue(*jac,IndexRe,ColRe,d_grad_P_dVj,INSERT_VALUES);
    MatSetValue(*jac,IndexIm,ColIm,d_grad_P_dVj,INSERT_VALUES);
    d_grad_P_dVi += -pcell->elen[j]/pcell->ilen[j]/pcell->area;
  }
  MatSetValue(*jac,IndexRe,IndexRe,d_grad_P_dVi,INSERT_VALUES);
  MatSetValue(*jac,IndexIm,IndexIm,d_grad_P_dVi,INSERT_VALUES);
  VecSetValue(*r,IndexRe,0,INSERT_VALUES);
  VecSetValue(*r,IndexIm,0,INSERT_VALUES);
}


void ElZone::AC1_ddm_om_contact(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
       DABC & bc,SMCZone *pz, int n)
{
  int  N = zofs[zofs.size()-1];
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
  VecSetValue(*r,zofs[zone_index]+i,0.0,INSERT_VALUES);
  VecSetValue(*r,N+zofs[zone_index]+i,0.0,INSERT_VALUES);
}

void ElZone::AC1_ddm_stk_contact(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
       DABC & bc,SMCZone *pz, int n)
{
  int  N = zofs[zofs.size()-1];
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
  VecSetValue(*r,zofs[zone_index]+i,0.0,INSERT_VALUES);
  VecSetValue(*r,N+zofs[zone_index]+i,0.0,INSERT_VALUES);
}

void ElZone::AC1_ddm_gate_contact(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
       DABC & bc,ISZone *pz, int n)
{
  int  N = zofs[zofs.size()-1];
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+1*n+0,-1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[pz->zone_index]+1*n+0,-1.0,INSERT_VALUES);
  VecSetValue(*r,zofs[zone_index]+i,0.0,INSERT_VALUES);
  VecSetValue(*r,N+zofs[zone_index]+i,0.0,INSERT_VALUES);
}

void ElZone::AC1_ddm_charge_contact(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
       DABC & bc,ISZone *pz, int n)
{
  int  N = zofs[zofs.size()-1];
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+1*n+0,-1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[pz->zone_index]+1*n+0,-1.0,INSERT_VALUES);
  VecSetValue(*r,zofs[zone_index]+i,0.0,INSERT_VALUES);
  VecSetValue(*r,N+zofs[zone_index]+i,0.0,INSERT_VALUES);
}
