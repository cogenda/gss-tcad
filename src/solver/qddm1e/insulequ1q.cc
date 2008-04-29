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
/*  Last update: Nov 23, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"

void ISZone::F1Q_ddm_inner(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs)
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


void ISZone::F1Q_ddm_gatebc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  const VoronoiCell* pcell = &pzone->davcell[i];
  GateBC *pbc = dynamic_cast<GateBC * >(bc.Get_pointer(pcell->bc_index-1));
  int    gate_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {gate_equ=j;break;}
  }
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[zone_index]+equ_num*size+gate_equ] + pbc->WorkFunction;
}


void ISZone::F1Q_ddm_semiconductor_insulator_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,SMCZone *pz, int n)
{
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[pz->zone_index]+5*n+0];
}

void ISZone::F1Q_ddm_electrode_insulator_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,ElZone *pz, int n)
{
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[pz->zone_index]+n];
}

void ISZone::F1Q_ddm_insulator_insulator_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,ISZone *pz, int n)
{
  if(zone_index > pz->pzone->zone_index)
  {
        f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[pz->zone_index]+n];
        return;
  }
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar Vi = x[zofs[zone_index]+i];
  PetscScalar grad_P = 0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vr = x[zofs[zone_index]+nb];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vr-Vi);
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar Vr  = x[zofs[pz->zone_index]+nb];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vr-Vi);
  }
  f[zofs[zone_index]+i] =  grad_P;
}

void ISZone::F1Q_ddm_chargebc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  const VoronoiCell* pcell = &pzone->davcell[i];
  ChargedContactBC *pbc = dynamic_cast<ChargedContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int    charge_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {charge_equ=j;break;}
  }
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - x[zofs[zone_index]+equ_num*size+charge_equ];
}


void ISZone::F1Q_gate_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  GateBC *pbc = dynamic_cast <GateBC * > (bc.Get_pointer(bc_index));
  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+node];     //potential of node i

    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+nb];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  if(pbc->electrode_type==VoltageBC)
  {
    PetscScalar V = pbc->Vapp;
    PetscScalar R = pbc->R;
    PetscScalar C = pbc->C;
    PetscScalar L = pbc->L;
    PetscScalar In = pbc->current;
    PetscScalar Icn = pbc->cap_current;
    PetscScalar Pn = pbc->potential;
    PetscScalar P = x[zofs[zone_index]+equ_num*size+i];
    f[zofs[zone_index]+equ_num*size+i]=(L/ODE_F.dt+R)*current-V+(1+(L/ODE_F.dt+R)*C/ODE_F.dt)*P
                                       -(L/ODE_F.dt+R)*C/ODE_F.dt*Pn-L/ODE_F.dt*(In+Icn);
    pbc->Set_Current_new(current);
  }
  pbc->Set_Potential_new(x[zofs[zone_index]+equ_num*size+i]);
}


void ISZone::F1Q_charge_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  ChargedContactBC *pbc = dynamic_cast <ChargedContactBC * > (bc.Get_pointer(bc_index));
  PetscScalar grad_P=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+node];     //potential of node i
    PetscScalar bc_len = (pcell->ilen[0]+pcell->ilen[pcell->nb_num-1])/2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+nb];     //potential of nb node
      grad_P += aux[node].eps*(Vj-Vi)/pcell->ilen[k]*pcell->elen[k];	
    }
  }
  f[zofs[zone_index]+equ_num*size+i]=grad_P+pbc->QF;
}

//------------------------------------------------------------------------------------------------

void ISZone::J1Q_ddm_inner(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs)
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


void ISZone::J1Q_ddm_gatebc(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);

  int equ_num = 1;
  int size = pzone->davcell.size();
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
     {gate_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+equ_num*size+gate_equ,-1.0,INSERT_VALUES);
}


void ISZone::J1Q_ddm_semiconductor_insulator_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                            SMCZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+5*n+0,-1.0,INSERT_VALUES);
}


void ISZone::J1Q_ddm_electrode_insulator_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                            ElZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+n,-1.0,INSERT_VALUES);
}

void ISZone::J1Q_ddm_insulator_insulator_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                            ISZone *pz, int n)
{
  if(zone_index > pz->pzone->zone_index)
  {
    MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+n,-1.0,INSERT_VALUES);
    return ;
  }
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar d_div_grad_P = 0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar value =  aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+nb,value,INSERT_VALUES);
    d_div_grad_P += -aux[i].eps*pcell->elen[j]/pcell->ilen[j];
  }
  
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int  nb = ncell->nb_array[j];
    PetscScalar value =  pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+i,zofs[pz->zone_index]+nb,value,INSERT_VALUES);
    d_div_grad_P += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
  }  
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,d_div_grad_P,INSERT_VALUES);  
}

void ISZone::J1Q_ddm_chargebc(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int charge_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
     {charge_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+equ_num*size+charge_equ,-1.0,INSERT_VALUES);
}


void ISZone::J1Q_gate_electrode(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  GateBC *pbc = dynamic_cast <GateBC * > (bc.Get_pointer(bc_index));

  int matrix_row = zofs[zone_index]+equ_num*size+i;
  PetscScalar R = pbc->R;
  PetscScalar C = pbc->C;
  PetscScalar L = pbc->L;
  PetscScalar e = mt->e;

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
      int node=bc[bc_index].psegment->node_array[j];
      const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
      for(int k=0;k<pcell->nb_num;k++)
      {
        int    nb = pcell->nb_array[k];
        //for displacement current
        PetscScalar dJ_dVi =  aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
        PetscScalar dJ_dVr = -aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
        if(pbc->electrode_type==VoltageBC)
        {
          MatSetValue(*jac,matrix_row,zofs[zone_index]+node,DeviceDepth*dJ_dVi*(L/ODE_F.dt+R),ADD_VALUES);
          MatSetValue(*jac,matrix_row,zofs[zone_index]+nb,DeviceDepth*dJ_dVr*(L/ODE_F.dt+R),ADD_VALUES);
        }
    }
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  if(pbc->electrode_type==VoltageBC)
        MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP

}


void ISZone::J1Q_charge_electrode(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  ChargedContactBC *pbc = dynamic_cast <ChargedContactBC * > (bc.Get_pointer(bc_index));
  int matrix_row = zofs[zone_index]+equ_num*size+i;

  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
      int node=bc[bc_index].psegment->node_array[j];
      const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
      PetscScalar dP_dVi=0;
      for(int k=0;k<pcell->nb_num;k++)
      {
        int    nb = pcell->nb_array[k];
	PetscScalar dP_dVj=aux[node].eps/pcell->ilen[k]*pcell->elen[k];
	dP_dVi+=-aux[node].eps/pcell->ilen[k]*pcell->elen[k];
        MatSetValue(*jac,matrix_row,zofs[zone_index]+nb,dP_dVj,ADD_VALUES);
      }
      MatSetValue(*jac,matrix_row,zofs[zone_index]+node,dP_dVi,ADD_VALUES);
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
}


//----------------------------------------------------------------------------------------------
void ISZone::F1Q_efield_update(PetscScalar *x,vector<int> & zofs, DABC &bc,vector<BZoneData *>zonedata)
{
  PetscScalar a=0,b=0,c=0,P_x=0,P_y=0,w=0;
  //calculate Ex Ey with least-squares gradient construction
  for(int i=0; i<pzone->davcell.size();i++)
  {
    a=0,b=0,c=0,P_x=0,P_y=0,w=0;
    const VoronoiCell *pcell = pzone->davcell.GetPointer(i);
    for(int k=0;k<pcell->nb_num;k++)
    {
      const VoronoiCell *ncell = pzone->davcell.GetPointer(pcell->nb_array[k]);
      PetscScalar dx = ncell->x - pcell->x;
      PetscScalar dy = ncell->y - pcell->y;
      PetscScalar dP = x[zofs[zone_index]+pcell->nb_array[k]]- x[zofs[zone_index]+i];
      w=1.0/sqrt(dx*dx+dy*dy);
      a+=w*w*dx*dx;
      b+=w*w*dx*dy;
      c+=w*w*dy*dy;
      P_x+=w*w*dP*dx;
      P_y+=w*w*dP*dy;
    }
    if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
    {
      InsulatorInterfaceBC *pbc;
      pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
      int n_zone = pbc->pinterface->Find_neighbor_zone_index(zone_index);
      int n_node = pbc->pinterface->Find_neighbor_node_index(zone_index,i);
      SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
      const VoronoiCell* dcell = pz->pzone->davcell.GetPointer(n_node);
      for(int k=1;k<dcell->nb_num-1;k++)
      {
        const VoronoiCell *ncell = pz->pzone->davcell.GetPointer(dcell->nb_array[k]);
        PetscScalar dx = ncell->x - dcell->x;
        PetscScalar dy = ncell->y - dcell->y;
        PetscScalar dP = x[zofs[n_zone]+5*dcell->nb_array[k]+0]- x[zofs[n_zone]+5*n_node+0];
        w=1.0/sqrt(dx*dx+dy*dy);
        a+=w*w*dx*dx;
        b+=w*w*dx*dy;
        c+=w*w*dy*dy;
        P_x+=w*w*dP*dx;
        P_y+=w*w*dP*dy;
      }
    }
    aux[i].Ex = -(c*P_x-b*P_y)/(a*c-b*b);
    aux[i].Ey = -(a*P_y-b*P_x)/(a*c-b*b);
  }
}
