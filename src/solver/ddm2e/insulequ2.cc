/*****************************************************************************/
/*   	        8888888         88888888         88888888                    */
/*  	      8                8                8                            */
/* 	     8                 8                8                            */
/*  	     8                  88888888         88888888                    */
/* 	     8      8888                8                8                   */
/* 	      8       8                 8                8                   */
/* 	        888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS 0.4x                                                                 */
/*  Last update: June 06, 2006                                               */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"

void ISZone::F2_ddm_inner(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar Vi = x[zofs[zone_index]+2*i+0];
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar grad_P=0, grad_T=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+2*nb+0];     //potential of nb node
    grad_P += pcell->elen[j]/pcell->ilen[j]*(Vj-Vi)/pcell->area;

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid    = 0.5*(Tj+Ti);
    PetscScalar dTdx_mid = (Tj-Ti)/pcell->ilen[j];
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    grad_T += kapa*pcell->elen[j]*dTdx_mid/pcell->area;
  }
  f[zofs[zone_index]+2*i+0] =  grad_P;
  f[zofs[zone_index]+2*i+1] =  grad_T;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r  = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }

}


void ISZone::F2_ddm_neumannbc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  NeumannBC *pbc = dynamic_cast <NeumannBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Vi = x[zofs[zone_index]+2*i+0];
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar grad_P=0, grad_T=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+2*nb+0];     //potential of nb node
    grad_P += pcell->elen[j]/pcell->ilen[j]*(Vj-Vi)/pcell->area;

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid    = 0.5*(Tj+Ti);
    PetscScalar dTdx_mid = (Tj-Ti)/pcell->ilen[j];
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    grad_T += kapa*pcell->elen[j]*dTdx_mid/pcell->area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      PetscScalar r = h*pbc->T_external;
      grad_T += 0.5*pcell->ilen[j]*(r-0.25*h*(3*Ti+Tj))/pcell->area;
    }
  }
  f[zofs[zone_index]+2*i+0] =  grad_P;
  f[zofs[zone_index]+2*i+1] =  grad_T;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r  = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }

}


void ISZone::F2_ddm_gatebc_segment(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  GateBC *pbc = dynamic_cast<GateBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar grad_T=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid    = 0.5*(Tj+Ti);
    PetscScalar dTdx_mid = (Tj-Ti)/pcell->ilen[j];
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    grad_T += kapa*pcell->elen[j]*dTdx_mid/pcell->area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      PetscScalar r = h*pbc->T_external;
      grad_T += 0.5*pcell->ilen[j]*(r-0.25*h*(3*Ti+Tj))/pcell->area;
    }
  }
  int equ_num = 2;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1)	{gate_equ=j;break;}
  }
  f[zofs[zone_index]+2*i+0] =  x[zofs[zone_index]+2*i+0] - x[zofs[zone_index]+equ_num*size+gate_equ] + pbc->WorkFunction;
  f[zofs[zone_index]+2*i+1] =  grad_T;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r  = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }

}


void ISZone::F2_ddm_gatebc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc,
                                     ElZone *pz, int n)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  GateBC *pbc = dynamic_cast<GateBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar grad_T=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid    = 0.5*(Tj+Ti);
    PetscScalar dTdx_mid = (Tj-Ti)/pcell->ilen[j];
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    grad_T += kapa*pcell->elen[j]*dTdx_mid;
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int   nb = ncell->nb_array[j];
    PetscScalar  Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];     //potential of nb node
    PetscScalar kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    grad_T += kapa*ncell->elen[j]/ncell->ilen[j]*(Tj_n-Ti);
  }

  int equ_num = 2;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1)	{gate_equ=j;break;}
  }
  f[zofs[zone_index]+2*i+0] =  x[zofs[zone_index]+2*i+0] - x[zofs[zone_index]+equ_num*size+gate_equ] + pbc->WorkFunction;
  f[zofs[zone_index]+2*i+1] =  grad_T;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r  = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity1*Tt/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*Tt/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity1*(Ti-fs[i].T)/ODE_F.dt*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*(Ti-fs[i].T)/ODE_F.dt*ncell->area;
    }
  }

}


void ISZone::F2_ddm_semiconductor_insulator_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,SMCZone *pz, int n)
{
  f[zofs[zone_index]+2*i+0] = x[zofs[zone_index]+2*i+0] - x[zofs[pz->zone_index]+4*n+0];
  f[zofs[zone_index]+2*i+1] = x[zofs[zone_index]+2*i+1] - x[zofs[pz->zone_index]+4*n+3];
}


void ISZone::F2_ddm_electrode_insulator_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,ElZone *pz, int n)
{
  f[zofs[zone_index]+2*i+0] = x[zofs[zone_index]+2*i+0] - x[zofs[pz->zone_index]+2*n+0];
  f[zofs[zone_index]+2*i+1] = x[zofs[zone_index]+2*i+1] - x[zofs[pz->zone_index]+2*n+1];
}


void ISZone::F2_ddm_insulator_insulator_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,ISZone *pz, int n)
{
  if(zone_index > pz->pzone->zone_index)
  {
    f[zofs[zone_index]+2*i+0] = x[zofs[zone_index]+2*i+0] - x[zofs[pz->zone_index]+2*n+0];
    f[zofs[zone_index]+2*i+1] = x[zofs[zone_index]+2*i+1] - x[zofs[pz->zone_index]+2*n+1];
    return;
  }

  PetscScalar Vi = x[zofs[zone_index]+2*i+0];
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar grad_P=0, grad_T=0;

  const VoronoiCell* pcell = &pzone->davcell[i];
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+2*nb+0];     //potential of nb node
    grad_P += aux[i].eps*pcell->elen[j]/pcell->ilen[j]*(Vj-Vi);

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid    = 0.5*(Tj+Ti);
    PetscScalar dTdx_mid = (Tj-Ti)/pcell->ilen[j];
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    grad_T += kapa*pcell->elen[j]*dTdx_mid;
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int   nb = ncell->nb_array[j];
    PetscScalar Vr  = x[zofs[pz->zone_index]+2*nb+0];     //potential of nb node
    grad_P += pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j]*(Vr-Vi);
    PetscScalar  Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];     //potential of nb node
    PetscScalar kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    grad_T += kapa*ncell->elen[j]/ncell->ilen[j]*(Tj_n-Ti);
  }

  f[zofs[zone_index]+2*i+0] =  grad_P;
  f[zofs[zone_index]+2*i+1] =  grad_T;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r  = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity1*Tt/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*Tt/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+2*i+1] += -aux[i].density*HeatCapacity1*(Ti-fs[i].T)/ODE_F.dt*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*(Ti-fs[i].T)/ODE_F.dt*ncell->area;
    }
  }
}


void ISZone::F2_gate_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 2;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  GateBC *pbc = dynamic_cast <GateBC * > (bc.Get_pointer(bc_index));
  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+2*node+0];     //potential of node i

    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+2*nb+0];     //potential of nb node
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


void ISZone::F2_ddm_chargebc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc,
                                       ElZone *pz, int n)
{
  int equ_num = 2;
  int size = pzone->davcell.size();
  const VoronoiCell* pcell = &pzone->davcell[i];
  ChargedContactBC *pbc = dynamic_cast<ChargedContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int    charge_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {charge_equ=j;break;}
  }
  f[zofs[zone_index]+2*i+0] = x[zofs[zone_index]+2*i+0] - x[zofs[zone_index]+equ_num*size+charge_equ];
  f[zofs[zone_index]+2*i+1] = x[zofs[zone_index]+2*i+1] - x[zofs[pz->zone_index]+2*n+1];
}

void ISZone::F2_charge_electrode(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  int equ_num = 2;
  int size = pzone->davcell.size();
  int bc_index = electrode[i];
  ChargedContactBC *pbc = dynamic_cast <ChargedContactBC * > (bc.Get_pointer(bc_index));
  PetscScalar grad_P=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+2*node];     //potential of node i
    PetscScalar bc_len = (pcell->ilen[0]+pcell->ilen[pcell->nb_num-1])/2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+2*nb];     //potential of nb node
      grad_P += aux[node].eps*(Vj-Vi)/pcell->ilen[k]*pcell->elen[k];
    }
  }
  f[zofs[zone_index]+equ_num*size+i]=grad_P+pbc->QF;
}

//----------------------------------------------------------------------------------------------


void ISZone::J2_ddm_inner(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar d_grad_P_dVi=0, d_grad_T_dTi=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar d_grad_P_dVj =  pcell->elen[j]/pcell->ilen[j]/pcell->area;
    MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*nb+0,d_grad_P_dVj,INSERT_VALUES);
    d_grad_P_dVi += -pcell->elen[j]/pcell->ilen[j]/pcell->area;

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid = 0.5*(Tj+Ti);
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    PetscScalar d_grad_T_dTj = kapa*pcell->elen[j]/pcell->ilen[j]/pcell->area;
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j]/pcell->area;
  }
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,d_grad_P_dVi,INSERT_VALUES);


  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,d_grad_T_dTi,INSERT_VALUES);

}


void ISZone::J2_ddm_neumannbc(int i,PetscScalar *x,Mat *jac,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  NeumannBC *pbc = dynamic_cast <NeumannBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar d_grad_P_dVi=0, d_grad_T_dTi=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar d_grad_P_dVj =  pcell->elen[j]/pcell->ilen[j]/pcell->area;
    MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*nb+0,d_grad_P_dVj,INSERT_VALUES);
    d_grad_P_dVi += -pcell->elen[j]/pcell->ilen[j]/pcell->area;

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid = 0.5*(Tj+Ti);
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    PetscScalar d_grad_T_dTj = kapa*pcell->elen[j]/pcell->ilen[j]/pcell->area;
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j]/pcell->area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/pcell->area;
      d_grad_T_dTj += -0.5*pcell->ilen[j]*0.25*h/pcell->area;
    }
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,d_grad_P_dVi,INSERT_VALUES);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,d_grad_T_dTi,INSERT_VALUES);

}


void ISZone::J2_ddm_gatebc_segment(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  GateBC *pbc = dynamic_cast<GateBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar d_grad_T_dTi=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid = 0.5*(Tj+Ti);
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    PetscScalar d_grad_T_dTj = kapa*pcell->elen[j]/pcell->ilen[j]/pcell->area;
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j]/pcell->area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/pcell->area;
      d_grad_T_dTj += -0.5*pcell->ilen[j]*0.25*h/pcell->area;
    }
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);
  int equ_num = 2;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {gate_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+equ_num*size+gate_equ,-1.0,INSERT_VALUES);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,d_grad_T_dTi,INSERT_VALUES);
}


void ISZone::J2_ddm_gatebc_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc,
                                     ElZone *pz, int n)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  GateBC *pbc = dynamic_cast<GateBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar d_grad_T_dTi=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid = 0.5*(Tj+Ti);
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    PetscScalar d_grad_T_dTj = kapa*pcell->elen[j]/pcell->ilen[j];
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }
  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar   Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];
    PetscScalar   kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    d_grad_T_dTi += -kapa*ncell->elen[j]/ncell->ilen[j];
    PetscScalar d_grad_T_dTj = kapa*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->pzone->zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }

  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);

  int equ_num = 2;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {gate_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+equ_num*size+gate_equ,-1.0,INSERT_VALUES);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi  += -aux[i].density*HeatCapacity1*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                       -pz->aux[n].density*HeatCapacity2*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity1/ODE_F.dt*pcell->area
                      -pz->aux[n].density*HeatCapacity2/ODE_F.dt*ncell->area;
    }
  }
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,d_grad_T_dTi,INSERT_VALUES);
}



void ISZone::J2_ddm_semiconductor_insulator_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
    SMCZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+4*n+0,-1.0,INSERT_VALUES);

  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->zone_index]+4*n+3,-1.0,INSERT_VALUES);
}

void ISZone::J2_ddm_electrode_insulator_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
    ElZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+2*n+0,-1.0,INSERT_VALUES);

  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->zone_index]+2*n+1,-1.0,INSERT_VALUES);
}

void ISZone::J2_ddm_insulator_insulator_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
    ISZone *pz, int n)
{
  if(zone_index > pz->pzone->zone_index)
  {
    MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+2*n+0,-1.0,INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,1.0,INSERT_VALUES);
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->zone_index]+2*n+1,-1.0,INSERT_VALUES);
    return;
  }
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar d_grad_P_dVi=0, d_grad_T_dTi=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar d_grad_P_dVj =   aux[i].eps*pcell->elen[j]/pcell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*nb+0,d_grad_P_dVj,INSERT_VALUES);
    d_grad_P_dVi += - aux[i].eps*pcell->elen[j]/pcell->ilen[j];

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //potential of nb node
    PetscScalar T_mid = 0.5*(Tj+Ti);
    PetscScalar kapa =  mt->thermal->HeatConduction(T_mid);
    PetscScalar d_grad_T_dTj = kapa*pcell->elen[j]/pcell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j];
  }

  const VoronoiCell* ncell = pz->pzone->davcell.GetPointer(n);
  pz->mt->mapping(&pz->pzone->danode[n],&pz->aux[n],ODE_F.clock);
  for(int j=0;j<ncell->nb_num;j++)
  {
    int    nb = ncell->nb_array[j];
    PetscScalar d_grad_P_dVj =  pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+2*nb+0,d_grad_P_dVj,INSERT_VALUES);
    d_grad_P_dVi += -pz->aux[n].eps*ncell->elen[j]/ncell->ilen[j];

    PetscScalar   Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];
    PetscScalar   kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    d_grad_T_dTi += -kapa*ncell->elen[j]/ncell->ilen[j];
    PetscScalar d_grad_T_dTj = kapa*ncell->elen[j]/ncell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->pzone->zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }

  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,d_grad_P_dVi,INSERT_VALUES);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity1*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                      -pz->aux[n].density*HeatCapacity2*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*ncell->area;

    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      d_grad_T_dTi += -aux[i].density*HeatCapacity1/ODE_F.dt*pcell->area
                      -pz->aux[n].density*HeatCapacity2/ODE_F.dt*ncell->area;
    }
  }

  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,d_grad_T_dTi,INSERT_VALUES);
}


void ISZone::J2_ddm_chargebc_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc,
                                       ElZone *pz, int n)
{
  int equ_num = 2;
  int size = pzone->davcell.size();
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int charge_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
    {charge_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+equ_num*size+charge_equ,-1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->zone_index]+2*n+1,-1.0,INSERT_VALUES);
}


void ISZone::J2_gate_electrode(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int equ_num = 2;
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
        MatSetValue(*jac,matrix_row,zofs[zone_index]+2*node+0,DeviceDepth*dJ_dVi*(L/ODE_F.dt+R),ADD_VALUES);
        MatSetValue(*jac,matrix_row,zofs[zone_index]+2*nb+0,DeviceDepth*dJ_dVr*(L/ODE_F.dt+R),ADD_VALUES);
      }
    }
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
  if(pbc->electrode_type==VoltageBC)
    MatSetValue(*jac,matrix_row,matrix_row,1+(L/ODE_F.dt+R)*C/ODE_F.dt,INSERT_VALUES); //dJ/dP

}


void ISZone::J2_charge_electrode(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  int equ_num = 2;
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
      MatSetValue(*jac,matrix_row,zofs[zone_index]+2*nb,dP_dVj,ADD_VALUES);
    }
    MatSetValue(*jac,matrix_row,zofs[zone_index]+2*node,dP_dVi,ADD_VALUES);
  }
  MatAssemblyBegin(*jac,MAT_FLUSH_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FLUSH_ASSEMBLY);
}


//--------------------------------------------------------------------------------------------------
void ISZone::F2_efield_update(PetscScalar *x,vector<int> & zofs, DABC &bc,vector<BZoneData *>zonedata)
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
      PetscScalar dP = x[zofs[zone_index]+2*pcell->nb_array[k]]- x[zofs[zone_index]+2*i];
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
        PetscScalar dP = x[zofs[n_zone]+4*dcell->nb_array[k]+0]- x[zofs[n_zone]+4*n_node+0];
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
