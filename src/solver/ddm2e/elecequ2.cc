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

void ElZone::F2_ddm_inner(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs)
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

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //Temp of nb node
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


void ElZone::F2_ddm_neumannbc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc)
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

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //Temp of nb node
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


void ElZone::F2_ddm_ombc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                   SMCZone *pz, int n)
{
  f[zofs[zone_index]+2*i+0] = x[zofs[zone_index]+2*i+0] - x[zofs[pz->zone_index]+4*n+0];
  f[zofs[zone_index]+2*i+1] = x[zofs[zone_index]+2*i+1] - x[zofs[pz->zone_index]+4*n+3];
}


void ElZone::F2_ddm_stkbc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                    SMCZone *pz, int n)
{
  f[zofs[zone_index]+2*i+0] = x[zofs[zone_index]+2*i+0] - x[zofs[pz->zone_index]+4*n+0];
  f[zofs[zone_index]+2*i+1] = x[zofs[zone_index]+2*i+1] - x[zofs[pz->zone_index]+4*n+3];
}


void ElZone::F2_ddm_gatebc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                     ISZone *pz, int n)
{
  f[zofs[zone_index]+2*i+0] = x[zofs[zone_index]+2*i+0] - x[zofs[pz->zone_index]+2*n+0];
  f[zofs[zone_index]+2*i+1] = x[zofs[zone_index]+2*i+1] - x[zofs[pz->zone_index]+2*n+1];
}


void ElZone::F2_ddm_elec_insulator_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc,
                                             ISZone *pz, int n)
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

    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //Temp of nb node
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
    PetscScalar  Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];     //Temp of nb node
    PetscScalar kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    grad_T += kapa*ncell->elen[j]/ncell->ilen[j]*(Tj_n-Ti);
  }
  
  f[zofs[zone_index]+2*i+0] =  grad_P;
  f[zofs[zone_index]+2*i+1] =  grad_T;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2 && ODE_F.BDF2_restart==false)
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

void ElZone::F2_ddm_chargebc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                     ISZone *pz, int n)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
  PetscScalar Ti = x[zofs[zone_index]+2*i+1];
  PetscScalar grad_T=0;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+2*nb+1];     //Temp of nb node
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
    PetscScalar  Tj_n = x[zofs[pz->pzone->zone_index]+2*nb+1];     //Temp of nb node
    PetscScalar kapa =  pz->mt->thermal->HeatConduction(0.5*(Ti+Tj_n));
    grad_T += kapa*ncell->elen[j]/ncell->ilen[j]*(Tj_n-Ti);
  }
  
  f[zofs[zone_index]+2*i+0] =  x[zofs[zone_index]+2*i+0] - x[zofs[pz->pzone->zone_index]+2*n+0];
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



//------------------------------------------------------------------------------------------------------

void ElZone::J2_ddm_inner(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs)
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


void ElZone::J2_ddm_neumannbc(int i,PetscScalar *x,Mat *jac,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
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




void ElZone::J2_ddm_ombc_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                   SMCZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+4*n+0,-1.0,INSERT_VALUES);

  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->zone_index]+4*n+3,-1.0,INSERT_VALUES);
}


void ElZone::J2_ddm_stkbc_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                   SMCZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+4*n+0,-1.0,INSERT_VALUES);

  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->zone_index]+4*n+3,-1.0,INSERT_VALUES);
}


void ElZone::J2_ddm_gatebc_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                     ISZone *pz, int n)
{
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[zone_index]+2*i+0,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+2*n+0,-1.0,INSERT_VALUES);

  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*i+1,1.0,INSERT_VALUES);
  MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[pz->zone_index]+2*n+1,-1.0,INSERT_VALUES);
}


void ElZone::J2_ddm_elec_insulator_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                   ISZone *pz, int n)
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
    PetscScalar d_grad_T_dTj = kapa*pcell->elen[j]/pcell->ilen[j];
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j];
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

void ElZone::J2_ddm_chargebc_interface(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
                                   ISZone *pz, int n)
{
  const VoronoiCell* pcell = &pzone->davcell[i];
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
    MatSetValue(*jac,zofs[zone_index]+2*i+1,zofs[zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j];
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
  MatSetValue(*jac,zofs[zone_index]+2*i+0,zofs[pz->zone_index]+2*n+0,-1.0,INSERT_VALUES);
  
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

