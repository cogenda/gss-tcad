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
/*  Last update: March 29, 2007                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"
#include "vec3.h"


//-----------------------------------------------------------------------------
// special process of electrode bcs for mix type simulation
//-----------------------------------------------------------------------------


void SMCZone::F1E_mix_ddm_ombc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  int z = zone_index;
  PetscScalar e  =  mt->e;
  PetscScalar kb =  mt->kb;
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar nie = mt->band->nie(fs[i].T);
  PetscScalar Nc  = mt->band->Nc(fs[i].T);
  PetscScalar Nv  = mt->band->Nv(fs[i].T);
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(pcell->bc_index-1));

  PetscScalar electron_density,hole_density;

  if(Fermi) //Fermi
  {
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) + kb*fs[i].T*log(aux[i].Nc));
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) - kb*fs[i].T*log(aux[i].Nv)+ aux[i].Eg);
    PetscScalar phin = pbc->Vapp;
    PetscScalar phip = pbc->Vapp;
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    f[zofs[z]+3*i+0] = Nc*fermi_half(etan) - Nv*fermi_half(etap)+ Na - Nd;
    f[zofs[z]+3*i+1] = ni - Nc*fermi_half(etan);
    f[zofs[z]+3*i+2] = pi - Nv*fermi_half(etap);
  }
  else     //Boltzmann
  {
    f[zofs[z]+3*i+0] = Vi - kb*fs[i].T/e*asinh((Nd-Na)/(2*nie)) + mt->kb*fs[i].T/e*log(Nc/nie)
                       + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) - pbc->Vapp;
    PetscScalar electron_density,hole_density;
    if(Na>Nd)   //p-type
    {
      hole_density = (-(Nd-Na)+sqrt((Nd-Na)*(Nd-Na)+4*nie*nie))/2.0;
      electron_density = nie*nie/hole_density;
    }
    else        //n-type
    {
      electron_density = ((Nd-Na)+sqrt((Nd-Na)*(Nd-Na)+4*nie*nie))/2.0;
      hole_density = nie*nie/electron_density;
    }
    f[zofs[z]+3*i+1] = ni - electron_density;  //electron density
    f[zofs[z]+3*i+2] = pi - hole_density;      //hole density
  }
}


void SMCZone::F1E_mix_ddm_stkbc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int>& zofs,DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);

  //Schottky current
  PetscScalar Fn = mt->band->SchottyJsn(ni,fs[i].T,pbc->WorkFunction-aux[i].affinity)*L;
  PetscScalar Fp = mt->band->SchottyJsp(pi,fs[i].T,pbc->WorkFunction-aux[i].affinity)*L;
  //Schotty Barrier Lowerring
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));

  f[zofs[zone_index]+3*i+0] = x[zofs[zone_index]+3*i+0] + pbc->WorkFunction - deltaVB - pbc->Vapp;
  f[zofs[zone_index]+3*i+1] = (f[zofs[zone_index]+3*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+3*i+2] = (f[zofs[zone_index]+3*i+2]-Fp)/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}


void SMCZone::F1E_mix_ddm_insulator_gate(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));

  PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+3*i+2];     //hole density of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);

  PetscScalar grad_P = 0;

  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
    //the poisson's equation on third boundary type
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar vgate = pbc->Vapp - pbc->WorkFunction;
      PetscScalar q = mt->e*pbc->QF; //sigma is the surface change density
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar r=q + eps_ox/Thick*vgate;
      PetscScalar s=eps_ox/Thick;
      grad_P += 0.5*pcell->ilen[j]*(r-0.25*s*(3*Vi+Vj));
    }
  }
  f[zofs[zone_index]+3*i+0] = (f[zofs[zone_index]+3*i+0]+grad_P)/pcell->area
                              + mt->e*((pi-aux[i].Na)+(aux[i].Nd-ni));
  f[zofs[zone_index]+3*i+1] = f[zofs[zone_index]+3*i+1]/pcell->area;
  f[zofs[zone_index]+3*i+2] = f[zofs[zone_index]+3*i+2]/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      f[zofs[zone_index]+3*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+3*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+3*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+3*i+2] += -(pi-fs[i].p)/ODE_F.dt;
    }
  }
}




void SMCZone::J1E_mix_ddm_ombc(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  PetscScalar    A[9];
  PetscInt       index[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(pcell->bc_index-1));
  index[0] = zofs[zone_index]+3*i+0;
  index[1] = zofs[zone_index]+3*i+1;
  index[2] = zofs[zone_index]+3*i+2;
  if(Fermi)  //Fermi
  {
    PetscScalar e  =  mt->e;
    PetscScalar kb =  mt->kb;
    PetscScalar Vt =  kb*fs[i].T/e;
    mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
    PetscScalar Nc  = mt->band->Nc(fs[i].T);
    PetscScalar Nv  = mt->band->Nv(fs[i].T);
    PetscScalar Vi = x[zofs[zone_index]+3*i+0];     //potential of node i
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) + kb*fs[i].T*log(aux[i].Nc));
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) - kb*fs[i].T*log(aux[i].Nv)+ aux[i].Eg);
    PetscScalar phin = pbc->Vapp;
    PetscScalar phip = pbc->Vapp;
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    PetscScalar detan_dVi =  1.0/Vt;
    PetscScalar detap_dVi = -1.0/Vt;
    A[0] = Nc*fermi_mhalf(etan)*detan_dVi - Nv*fermi_mhalf(etap)*detap_dVi;
    A[1] = 0;
    A[2] = 0;
    A[3] = -Nc*fermi_mhalf(etan)*detan_dVi;
    A[4] = 1;
    A[5] = 0;
    A[6] = -Nv*fermi_mhalf(etap)*detap_dVi;
    A[7] = 0;
    A[8] = 1;
    MatSetValues(*jac,3,index,3,index,A,INSERT_VALUES);

    //f[zofs[z]+3*i+0] = Nc*fermi_half(etan) - Nv*fermi_half(etap)+ Na - Nd;
    //f[zofs[z]+3*i+1] = ni - Nc*fermi_half(etan);
    //f[zofs[z]+3*i+2] = pi - Nv*fermi_half(etap);
  }
  else
  {
    A[0] = 1;  A[1] = 0;  A[2] = 0;
    A[3] = 0;  A[4] = 1;  A[5] = 0;
    A[6] = 0;  A[7] = 0;  A[8] = 1;
    MatSetValues(*jac,3,index,3,index,A,INSERT_VALUES);
  }
}


void SMCZone::J1E_mix_ddm_stkbc(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat3    A;
  PetscInt       index[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int  z = zone_index;

  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  PetscScalar Vi = x[zofs[z]+3*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar area = pcell->area;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fp_dpi = 0;
  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,index,3,col,A.m);
    A.m[0] = 0.0;
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn(ni,fs[i].T,VB-deltaVB)*0.5*pcell->ilen[j];
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp(pi,fs[i].T,VB+deltaVB)*0.5*pcell->ilen[j];
    }
    MatSetValues(*jac,3,index,3,col,A.m,INSERT_VALUES);
  }

  MatGetValues(*jtmp,3,index,3,index,A.m);
  A=A/area;
  A.m[0] =  1.0;                      //dfun(0)/dP(i)
  A.m[1] =  0.0;                      //dfun(0)/dn(i)
  A.m[2] =  0.0;                      //dfun(0)/dp(i)

  A.m[4] +=   d_Fn_dni/area;          //dfun(1)/dn(i)
  A.m[8] +=  -d_Fp_dpi/area;          //dfun(2)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[4] += -1/ODE_F.dt;
      A.m[8] += -1/ODE_F.dt;
    }
  }
  MatSetValues(*jac,3,index,3,index,A.m,INSERT_VALUES);
}


void SMCZone::J1E_mix_ddm_insulator_gate(int i,PetscScalar *x, Mat *jac, Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  Mat3    A;
  PetscInt       index[3],col[3];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;

  PetscScalar ni = x[zofs[z]+3*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+3*i+2];     //hole density of node i

  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;

  //--------------------------------
  index[0] = zofs[z]+3*i+0;
  index[1] = zofs[z]+3*i+1;
  index[2] = zofs[z]+3*i+2;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    col[0] = zofs[z]+3*nb+0;
    col[1] = zofs[z]+3*nb+1;
    col[2] = zofs[z]+3*nb+2;
    MatGetValues(*jtmp,3,index,3,col,A.m);
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar s=eps_ox/Thick;
      d_grad_P_dVi += -0.5*pcell->ilen[j]*0.25*s*3;
      A.m[0]+= -0.5*pcell->ilen[j]*0.25*s/area;
    }
    MatSetValues(*jac,3,index,3,col,A.m,INSERT_VALUES);
  }

  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);

  MatGetValues(*jtmp,3,index,3,index,A.m);
  A=A/area;
  A.m[0] +=  d_grad_P_dVi/area;             //dfun(0)/dP(i)
  A.m[1] +=  -mt->e;   //dfun(0)/dn(i)
  A.m[2] +=   mt->e;   //dfun(0)/dp(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[4] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[8] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[4] += -1/ODE_F.dt;
      A.m[8] += -1/ODE_F.dt;
    }
  }
  MatSetValues(*jac,3,index,3,index,A.m,INSERT_VALUES);
}


void SMCZone::F1E_mix_om_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int bc_index = electrode[i];
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  PetscScalar e  =  mt->e;

  //calculate the total current of ohmic electrode
  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node i
    //conduction current
    current += DeviceDepth*(f[zofs[zone_index]+3*node+1]-f[zofs[zone_index]+3*node+2]);
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  pbc->Set_Current_new(current);

}

void SMCZone::F1E_mix_stk_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));
  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node i
    PetscScalar ni = x[zofs[zone_index]+3*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+3*node+2];     //hole density of node i
    mt->mapping(&pzone->danode[node],&aux[node],ODE_F.clock);
    PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[node].eps,sqrt(aux[node].Ex*aux[node].Ex+aux[node].Ey*aux[node].Ey));
    PetscScalar Jsn = mt->band->SchottyJsn(ni,fs[node].T,pbc->WorkFunction-aux[node].affinity-deltaVB);
    PetscScalar Jsp = mt->band->SchottyJsp(pi,fs[node].T,pbc->WorkFunction-aux[node].affinity+deltaVB);
    current += -Jsn*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsn*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  pbc->Set_Current_new(current);

}


void SMCZone::F1E_mix_ins_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));
  PetscScalar current=0;

  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node i
    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }

  pbc->Set_Current_new(current);
}


void SMCZone::F1E_mix_om_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
                                        PetscScalar *x, Mat *jtmp, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(bc_index));
  pdI_pdV = 0;
  PetscScalar e = mt->e;
  PetscScalar kb =  mt->kb;
  VecZeroEntries(pdI_pdw);
  VecZeroEntries(pdF_pdV);
  PetscScalar * apdI_pdw;
  PetscScalar * apdF_pdV;
  VecGetArray(pdI_pdw,&apdI_pdw);
  VecGetArray(pdF_pdV,&apdF_pdV);
  current = pbc->Get_Current_new();
  PetscScalar    A1[3],A2[3];
  PetscInt       index[3],col[3];
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node i
    PetscScalar dJdisp_dVi=0,dJdisp_dVr=0;
    index[0] = zofs[zone_index]+3*node+0;
    index[1] = zofs[zone_index]+3*node+1;
    index[2] = zofs[zone_index]+3*node+2;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+3*nb+0];     //potential of nb node
      col[0] = zofs[zone_index]+3*nb+0;
      col[1] = zofs[zone_index]+3*nb+1;
      col[2] = zofs[zone_index]+3*nb+2;

      MatGetValues(*jtmp,1,&index[1],3,col,A1);
      MatGetValues(*jtmp,1,&index[2],3,col,A2);

      //for displacement current
      dJdisp_dVi += aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
      dJdisp_dVr = -aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];

      //
      apdI_pdw[col[0]] += (A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
      apdI_pdw[col[1]] += (A1[1]-A2[1])*DeviceDepth;
      apdI_pdw[col[2]] += (A1[2]-A2[2])*DeviceDepth;
    }

    MatGetValues(*jtmp,1,&index[1],3,index,A1);
    MatGetValues(*jtmp,1,&index[2],3,index,A2);

    apdI_pdw[index[0]] += (A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
    apdI_pdw[index[1]] += (A1[1]-A2[1])*DeviceDepth;
    apdI_pdw[index[2]] += (A1[2]-A2[2])*DeviceDepth;
    
    if(Fermi)
    {
      PetscScalar Vt =  kb*fs[node].T/e;
      mt->mapping(&pzone->danode[node],&aux[node],ODE_F.clock);
      PetscScalar Nc  = mt->band->Nc(fs[node].T);
      PetscScalar Nv  = mt->band->Nv(fs[node].T);
      PetscScalar Vi = x[zofs[zone_index]+3*node+0];     //potential of node 
      PetscScalar Ec =  -(e*Vi + aux[node].affinity + mt->band->EgNarrowToEc(fs[node].T) + kb*fs[node].T*log(aux[node].Nc));
      PetscScalar Ev =  -(e*Vi + aux[node].affinity - mt->band->EgNarrowToEv(fs[node].T) - kb*fs[node].T*log(aux[node].Nv)+ aux[node].Eg);
      PetscScalar phin = pbc->Vapp;
      PetscScalar phip = pbc->Vapp;
      PetscScalar etan = (-e*phin-Ec)/kb/fs[node].T;
      PetscScalar etap = (Ev+e*phip)/kb/fs[node].T;
      PetscScalar detan_dVi =  1.0/Vt;
      PetscScalar detap_dVi = -1.0/Vt;
      apdF_pdV[index[0]] = Nc*fermi_mhalf(etan)*detan_dVi - Nv*fermi_mhalf(etap)*detap_dVi;
      apdF_pdV[index[1]] = -Nc*fermi_mhalf(etan)*detan_dVi;
      apdF_pdV[index[2]] = -Nv*fermi_mhalf(etap)*detap_dVi;
    }
    else
      apdF_pdV[index[0]] = 1.0;
  }
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F1E_mix_stk_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
    PetscScalar *x, Mat *jtmp, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));
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
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar VB = pbc->WorkFunction-aux[node].affinity;
    PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[node].eps,sqrt(aux[node].Ex*aux[node].Ex+aux[node].Ey*aux[node].Ey));
    PetscScalar ni = x[zofs[zone_index]+3*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+3*node+2];     //hole density of node i
    PetscScalar dJ_dni = -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB-deltaVB)*0.5*pcell->ilen[0]
                         -mt->band->pdSchottyJsn_pdn(ni,fs[node].T,VB-deltaVB)*0.5*pcell->ilen[pcell->nb_num-1];
    PetscScalar dJ_dpi = -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB+deltaVB)*0.5*pcell->ilen[0]
                         -mt->band->pdSchottyJsp_pdp(pi,fs[node].T,VB+deltaVB)*0.5*pcell->ilen[pcell->nb_num-1];

    //for displacement current
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar dI_dVi =  DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
      PetscScalar dI_dVr = -DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
      apdI_pdw[zofs[zone_index]+3*node+0] += dI_dVi;
      apdI_pdw[zofs[zone_index]+3*nb+0]   += dI_dVr;
    }
    apdI_pdw[zofs[zone_index]+3*node+1] = DeviceDepth*dJ_dni;
    apdI_pdw[zofs[zone_index]+3*node+2] = DeviceDepth*dJ_dpi;

    apdF_pdV[zofs[zone_index]+3*node+0] = 1.0;
  }

  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F1E_mix_ins_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
    PetscScalar *x, Mat *jtmp, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));
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
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);

    //for displacement current
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar dI_dVi =  DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
      PetscScalar dI_dVr = -DeviceDepth*pcell->elen[k]*aux[node].eps/pcell->ilen[k]/ODE_F.dt;
      apdI_pdw[zofs[zone_index]+3*node+0] += dI_dVi;
      apdI_pdw[zofs[zone_index]+3*nb+0]   += dI_dVr;
      if(k==0||k==pcell->nb_num-1)
      {
        PetscScalar Thick = pbc->Thick;
        PetscScalar eps_ox = mt->eps0*pbc->eps;
        PetscScalar s=eps_ox/Thick;
        apdF_pdV[zofs[zone_index]+3*node+0] += 0.5*pcell->ilen[k]*s/pcell->area;
      }
    }
  }

  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}

