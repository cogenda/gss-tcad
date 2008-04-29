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

void ISZone::AC1_ddm_inner(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs)
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


void ISZone::AC1_ddm_electrode_insulator_interface(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
       DABC & bc,ElZone *pz, int n)
{
  int  N = zofs[zofs.size()-1];
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+n,-1.0,INSERT_VALUES);
  MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[pz->zone_index]+n,-1.0,INSERT_VALUES);
  VecSetValue(*r,zofs[zone_index]+i,0.0,INSERT_VALUES);
  VecSetValue(*r,N+zofs[zone_index]+i,0.0,INSERT_VALUES);
}


void ISZone::AC1_ddm_semiconductor_insulator_interface(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
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


void ISZone::AC1_ddm_insulator_insulator_interface(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
       DABC & bc,ISZone *pz, int n)
{
  int  N = zofs[zofs.size()-1];
  if(zone_index > pz->pzone->zone_index)
  {
    MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
    MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[zone_index]+i,1.0,INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
    MatSetValue(*jac,N+zofs[zone_index]+i,N+zofs[pz->zone_index]+3*n+0,-1.0,INSERT_VALUES);
    VecSetValue(*r,zofs[zone_index]+i,0.0,INSERT_VALUES);
    VecSetValue(*r,N+zofs[zone_index]+i,0.0,INSERT_VALUES);
    return;
  }
  PetscScalar d_grad_P_dVi = 0;
  int IndexRe = zofs[zone_index]+i;
  int IndexIm = N+zofs[zone_index]+i;
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    int ColRe = zofs[zone_index]+nb;
    int ColIm = N+zofs[zone_index]+nb;
    PetscScalar d_grad_P_dVj =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    MatSetValue(*jac,IndexRe,ColRe,d_grad_P_dVj,INSERT_VALUES);
    MatSetValue(*jac,IndexIm,ColIm,d_grad_P_dVj,INSERT_VALUES);
    d_grad_P_dVi += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
  }

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    int ColRe = zofs[pz->zone_index]+nb;
    int ColIm = N+zofs[pz->zone_index]+nb;
    PetscScalar d_grad_P_dVj =  pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,IndexRe,ColRe,d_grad_P_dVj,INSERT_VALUES);
    MatSetValue(*jac,IndexIm,ColIm,d_grad_P_dVj,INSERT_VALUES);
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
  }
  MatSetValue(*jac,IndexRe,IndexRe,d_grad_P_dVi,INSERT_VALUES);
  MatSetValue(*jac,IndexIm,IndexIm,d_grad_P_dVi,INSERT_VALUES); 
  
  VecSetValue(*r,IndexRe,0,INSERT_VALUES);
  VecSetValue(*r,IndexIm,0,INSERT_VALUES);
}


void ISZone::AC1_ddm_gatebc(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,DABC & bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int  N = zofs[zofs.size()-1];
  int IndexRe = zofs[zone_index]+i;
  int IndexIm = N+zofs[zone_index]+i;
  MatSetValue(*jac,IndexRe,IndexRe,1.0,INSERT_VALUES);
  MatSetValue(*jac,IndexIm,IndexIm,1.0,INSERT_VALUES);

  int equ_num = 1;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
     {gate_equ=j;break;}
  MatSetValue(*jac,IndexRe,zofs[zone_index]+equ_num*size+gate_equ,-1.0,INSERT_VALUES);
  MatSetValue(*jac,IndexIm,N+zofs[zone_index]+equ_num*size+gate_equ,-1.0,INSERT_VALUES);
  VecSetValue(*r,IndexRe,0.0,INSERT_VALUES);
  VecSetValue(*r,IndexIm,0.0,INSERT_VALUES);
}


void ISZone::AC1_ddm_chargebc(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,DABC & bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int  N = zofs[zofs.size()-1];
  int IndexRe = zofs[zone_index]+i;
  int IndexIm = N+zofs[zone_index]+i;
  MatSetValue(*jac,IndexRe,IndexRe,1.0,INSERT_VALUES);
  MatSetValue(*jac,IndexIm,IndexIm,1.0,INSERT_VALUES);

  int equ_num = 1;
  int size = pzone->davcell.size();
  int charge_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
     {charge_equ=j;break;}
  MatSetValue(*jac,IndexRe,zofs[zone_index]+equ_num*size+charge_equ,-1.0,INSERT_VALUES);
  MatSetValue(*jac,IndexIm,N+zofs[zone_index]+equ_num*size+charge_equ,-1.0,INSERT_VALUES);
  VecSetValue(*r,IndexRe,0.0,INSERT_VALUES);
  VecSetValue(*r,IndexIm,0.0,INSERT_VALUES);
}

void ISZone::AC1_gate_electrode(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs,
      DABC & bc, PetscScalar DeviceDepth )
{
  int  equ_num = 1;
  int  size = pzone->davcell.size();
  int  z = zone_index;
  int  N = zofs[zofs.size()-1];
  int  bc_index = electrode[i];
  GateBC *pbc = dynamic_cast <GateBC * > (bc.Get_pointer(bc_index));
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;

  complex <PetscScalar> Z1(R,omega*L);
  complex <PetscScalar> Y2(0,omega*C);
  complex <PetscScalar> A,I(0,1);
  PetscInt IndexRe,IndexIm,ColRe,ColIm;
  PetscInt EquRe = zofs[z]+equ_num*size+i;
  PetscInt EquIm = N + zofs[z]+equ_num*size+i;
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    complex <PetscScalar> dJdisp_dVi=0,dJdisp_dVr=0;
    IndexRe = zofs[z]+node;
    IndexIm = N + zofs[z]+node;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      ColRe = zofs[z]+nb;
      ColIm = N + zofs[z]+nb;

      //for displacement current
      dJdisp_dVi += omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;
      dJdisp_dVr = -omega*aux[node].eps*pcell->elen[k]/pcell->ilen[k]*I;

      A = dJdisp_dVr*DeviceDepth*Z1;
      MatSetValue(*jac,EquRe,ColRe,A.real(),ADD_VALUES);
      MatSetValue(*jac,EquRe,ColIm,-A.imag(),ADD_VALUES);
      MatSetValue(*jac,EquIm,ColRe,A.imag(),ADD_VALUES);
      MatSetValue(*jac,EquIm,ColIm,A.real(),ADD_VALUES);
    }
    A = dJdisp_dVi*DeviceDepth*Z1;
    MatSetValue(*jac,EquRe,IndexRe,A.real(),ADD_VALUES);
    MatSetValue(*jac,EquRe,IndexIm,-A.imag(),ADD_VALUES);
    MatSetValue(*jac,EquIm,IndexRe,A.imag(),ADD_VALUES);
    MatSetValue(*jac,EquIm,IndexIm,A.real(),ADD_VALUES);
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

  complex <PetscScalar> K = Z1*Y2 + PetscScalar(1.0);
  MatSetValue(*jac,EquRe,EquRe,K.real(),INSERT_VALUES);
  MatSetValue(*jac,EquRe,EquIm,-K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquRe,K.imag(),INSERT_VALUES);
  MatSetValue(*jac,EquIm,EquIm,K.real(),INSERT_VALUES);

  VecSetValue(*r,EquRe,pbc->Vac, INSERT_VALUES);
  VecSetValue(*r,EquIm,0.0, INSERT_VALUES);
}


void ISZone::AC1_charge_electrode(int i, PetscScalar omega, PetscScalar *x, Vec *r, Mat *jac, vector<int> &zofs, DABC & bc)
{
  int  equ_num = 1;
  int  size = pzone->davcell.size();
  int  z = zone_index;
  int  N = zofs[zofs.size()-1];
  int  bc_index = electrode[i];
  ChargedContactBC *pbc = dynamic_cast <ChargedContactBC * > (bc.Get_pointer(bc_index));

  PetscInt IndexRe,IndexIm,ColRe,ColIm;
  PetscInt EquRe = zofs[z]+equ_num*size+i;
  PetscInt EquIm = N + zofs[z]+equ_num*size+i;
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar dP_dVi=0;
    IndexRe = zofs[z]+node;
    IndexIm = N + zofs[z]+node;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      ColRe = zofs[z]+nb;
      ColIm = N + zofs[z]+nb;
      PetscScalar dP_dVj=aux[node].eps/pcell->ilen[k]*pcell->elen[k];
      MatSetValue(*jac,IndexRe,ColRe,dP_dVj,ADD_VALUES);
      MatSetValue(*jac,IndexIm,ColIm,dP_dVj,ADD_VALUES);
      dP_dVi+=-aux[node].eps/pcell->ilen[k]*pcell->elen[k];
    }
    MatSetValue(*jac,IndexRe,IndexRe,dP_dVi,ADD_VALUES);
    MatSetValue(*jac,IndexIm,IndexIm,dP_dVi,ADD_VALUES);
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);

  VecSetValue(*r,EquRe,0.0, INSERT_VALUES);
  VecSetValue(*r,EquIm,0.0, INSERT_VALUES);
}


void ISZone::AC1_gate_electrode_current(int i, PetscScalar omega, PetscScalar *x, Vec *v, vector<int> &zofs,
      DABC & bc, PetscScalar DeviceDepth, PetscScalar &IacRe, PetscScalar &IacIm)
{
  int   equ_num = 1;
  int   size = pzone->davcell.size();
  int   z = zone_index;
  int   N = zofs[zofs.size()-1];
  int   bc_index = electrode[i];
  GateBC *pbc = dynamic_cast <GateBC * > (bc.Get_pointer(bc_index));
  PetscInt       IndexRe,IndexIm,ColRe,ColIm;
  PetscScalar    *vv;
  VecGetArray(*v,&vv);
  IacRe = 0;
  IacIm = 0;
  //for displacement current
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    IndexRe = zofs[z]+node;
    IndexIm = N + zofs[z]+node;
    PetscScalar ViRe = vv[IndexRe];
    PetscScalar ViIm = vv[IndexIm];

    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      ColRe = zofs[z]+nb;
      ColIm = N + zofs[z]+nb;

      PetscScalar VjRe = vv[ColRe];
      PetscScalar VjIm = vv[ColIm];

      IacRe += -omega*aux[node].eps*(ViIm-VjIm)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
      IacIm +=  omega*aux[node].eps*(ViRe-VjRe)/pcell->ilen[k]*pcell->elen[k]*DeviceDepth;
    }
  }
  VecRestoreArray(*v,&vv);
}
