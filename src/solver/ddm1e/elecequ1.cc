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
/*  Last update: June 27, 2007                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"

void ElZone::F1_ddm_inner(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar Vi = x[zofs[zone_index]+i];
  PetscScalar div_grad_P = 0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vr = x[zofs[zone_index]+nb];     //potential of nb node
    div_grad_P += pcell->elen[j]/pcell->ilen[j]*(Vr-Vi)/pcell->area;
  }
  f[zofs[zone_index]+i] =  div_grad_P;
}

void ElZone::F1_ddm_om_contact(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,SMCZone *pz, int n)
{
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[pz->zone_index]+3*n+0];
}

void ElZone::F1_ddm_stk_contact(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,SMCZone *pz, int n)
{
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[pz->zone_index]+3*n+0];
}

void ElZone::F1_ddm_gate_contact(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,ISZone *pz, int n)
{
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[pz->zone_index]+n+0];
}

void ElZone::F1_ddm_charge_contact(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,ISZone *pz, int n)
{
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[pz->zone_index]+n+0];
}

//------------------------------------------------------------------------------------------------

void ElZone::J1_ddm_inner(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar d_div_grad_P = 0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar value =  pcell->elen[j]/pcell->ilen[j]/pcell->area;
    MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+nb,value,INSERT_VALUES);
    d_div_grad_P += -pcell->elen[j]/pcell->ilen[j]/pcell->area;
  }
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,d_div_grad_P,INSERT_VALUES);
}

void ElZone::J1_ddm_om_contact(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                            SMCZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
}

void ElZone::J1_ddm_stk_contact(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                            SMCZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
}

void ElZone::J1_ddm_gate_contact(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                            ISZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+n,-1.0,INSERT_VALUES);
}

void ElZone::J1_ddm_charge_contact(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                            ISZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+n,-1.0,INSERT_VALUES);
}
