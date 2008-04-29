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
/*  Last update: May 11, 2005                                                */
/*                                                                           */
/*  Gong Ding gdiso@ustc.edu                                                 */
/*  Xuan Chun xiaomoyu505@163.com                                            */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"

//-----------------------------------------------------------------------------
// special process of electrode bcs for mix type simulation
//-----------------------------------------------------------------------------

void ISZone::F2_mix_ddm_gatebc_segment(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
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

  f[zofs[zone_index]+2*i+0] =  x[zofs[zone_index]+2*i+0] - pbc->Vapp+ pbc->WorkFunction;
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

  //just for fill the extra equ entry 
  int equ_num = 2;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {gate_equ=j;break;}
  }
  f[zofs[zone_index]+equ_num*size+gate_equ]=0.0;
}


void ISZone::F2_mix_ddm_gatebc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc,
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
  f[zofs[zone_index]+2*i+0] =  x[zofs[zone_index]+2*i+0] - pbc->Vapp + pbc->WorkFunction;
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
  
  //just for fill the extra equ entry 
  f[zofs[zone_index]+equ_num*size+gate_equ]=0.0;
}


void ISZone::F2_mix_gate_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
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
  pbc->Set_Current_new(current);
}

//-------------------------------------------------------------------------------------------------

void ISZone::J2_mix_ddm_gatebc_segment(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
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

  
  //just for fill the matrix entry 
  int equ_num = 2;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {gate_equ=j;break;}
  }
  MatSetValue(*jac,zofs[zone_index]+equ_num*size+gate_equ,zofs[zone_index]+equ_num*size+gate_equ,1.0,INSERT_VALUES);
}


void ISZone::J2_mix_ddm_gatebc_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc,
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
    
  //just for fill the matrix entry 
  int equ_num = 2;
  int size = pzone->davcell.size();
  int gate_equ;
  for(int j=0;j<electrode.size();j++)
    if(electrode[j]==pcell->bc_index-1)
     {gate_equ=j;break;}
  MatSetValue(*jac,zofs[zone_index]+equ_num*size+gate_equ,zofs[zone_index]+equ_num*size+gate_equ,1.0,INSERT_VALUES);
}


void ISZone::F2_mix_gate_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,PetscScalar *x, Mat *jac, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  GateBC *pbc = dynamic_cast <GateBC * > (bc.Get_pointer(bc_index));
  pdI_pdV = 0;
  PetscScalar e = mt->e;
  VecZeroEntries(pdI_pdw);
  VecZeroEntries(pdF_pdV);
  PetscScalar * apdI_pdw;
  PetscScalar * apdF_pdV;
  VecGetArray(pdI_pdw,&apdI_pdw);
  VecGetArray(pdF_pdV,&apdF_pdV);
  current = pbc->Get_Current_new();

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
	apdI_pdw[zofs[zone_index]+2*node] += dJ_dVi*DeviceDepth;
	apdI_pdw[zofs[zone_index]+2*nb]   += dJ_dVr*DeviceDepth;
      }
      apdF_pdV[zofs[zone_index]+2*node] = 1.0;
      pdI_pdV = 0;
  }
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


