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
/*  Last update: May 11, 2007                                                */
/*                                                                           */
/*  Gong Ding gdiso@ustc.edu                                                 */
/*  Xuan Chun xiaomoyu505@163.com                                            */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "zonedata.h"
#include "vec4.h"


//-----------------------------------------------------------------------------
// special process of electrode bcs for mix type simulation
//-----------------------------------------------------------------------------


void SMCZone::F2E_mix_ddm_ombc_segment(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar kb = mt->kb;
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = fabs(x[zofs[z]+4*i+3]);
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar nie = mt->band->nie(fs[i].T);
  PetscScalar Nc  = mt->band->Nc(fs[i].T);
  PetscScalar Nv  = mt->band->Nv(fs[i].T);
  PetscScalar H = 0;
  PetscScalar electron_density,hole_density;

  if(Fermi) //Fermi
  {
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) + kb*fs[i].T*log(aux[i].Nc));
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) - kb*fs[i].T*log(aux[i].Nv)+ aux[i].Eg);
    PetscScalar phin = pbc->Vapp;
    PetscScalar phip = pbc->Vapp;
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    f[zofs[z]+4*i+0] = Nc*fermi_half(etan) - Nv*fermi_half(etap)+ Na - Nd;
    f[zofs[z]+4*i+1] = ni - Nc*fermi_half(etan);
    f[zofs[z]+4*i+2] = pi - Nv*fermi_half(etap);
  }
  else
  {
    f[zofs[z]+4*i+0] =  Vi - kb*fs[i].T/mt->e*asinh((Nd-Na)/(2*nie))  + mt->kb*fs[i].T/mt->e*log(Nc/nie)
                        + aux[i].affinity  + mt->band->EgNarrowToEc(fs[i].T) - pbc->Vapp;

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
    f[zofs[z]+4*i+1] = x[zofs[z]+4*i+1] - electron_density;  //electron density
    f[zofs[z]+4*i+2] = x[zofs[z]+4*i+2] - hole_density;      //hole density
  }
  PetscScalar grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];            //lattice temperature of nb node
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
  f[zofs[zone_index]+4*i+3] =   grad_T + H;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }

}


void SMCZone::F2E_mix_ddm_ombc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc,
    ElZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar e  =  mt->e;
  PetscScalar kb = mt->kb;
  PetscScalar Na = aux[i].Na;
  PetscScalar Nd = aux[i].Nd;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = fabs(x[zofs[z]+4*i+3]);
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar nie = mt->band->nie(fs[i].T);
  PetscScalar Nc  = mt->band->Nc(fs[i].T);
  PetscScalar Nv  = mt->band->Nv(fs[i].T);
  PetscScalar H = 0;
  PetscScalar electron_density,hole_density;

  if(Fermi) //Fermi
  {
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) + kb*fs[i].T*log(aux[i].Nc));
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) - kb*fs[i].T*log(aux[i].Nv)+ aux[i].Eg);
    PetscScalar phin = pbc->Vapp;
    PetscScalar phip = pbc->Vapp;
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    f[zofs[z]+4*i+0] = Nc*fermi_half(etan) - Nv*fermi_half(etap)+ Na - Nd;
    f[zofs[z]+4*i+1] = ni - Nc*fermi_half(etan);
    f[zofs[z]+4*i+2] = pi - Nv*fermi_half(etap);
  }
  {
    f[zofs[z]+4*i+0] =  Vi - kb*fs[i].T/mt->e*asinh((Nd-Na)/(2*nie))  + mt->kb*fs[i].T/mt->e*log(Nc/nie)
                        + aux[i].affinity  + mt->band->EgNarrowToEc(fs[i].T) - pbc->Vapp;

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
    f[zofs[z]+4*i+1] = x[zofs[z]+4*i+1] - electron_density;  //electron density
    f[zofs[z]+4*i+2] = x[zofs[z]+4*i+2] - hole_density;      //hole density
  }
  PetscScalar grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];            //lattice temperature of nb node
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
  f[zofs[zone_index]+4*i+3] =   grad_T + H;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*Tt/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*Tt/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*(Ti-fs[i].T)/ODE_F.dt*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*(Ti-fs[i].T)/ODE_F.dt*ncell->area;
    }
  }

}


void SMCZone::F2E_mix_ddm_stkbc_segment(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int>& zofs,DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);

  //Schotty Barrier Lowerring
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  //potential
  f[zofs[z]+4*i+0] = x[zofs[z]+4*i+0] + pbc->WorkFunction - deltaVB - pbc->Vapp;

  PetscScalar Fn=0,Fp=0,grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];     //lattice temperature of nb node
    if(j==0||j==pcell->nb_num-1)
    {
      Fn += mt->band->SchottyJsn(ni,Ti,pbc->WorkFunction-aux[i].affinity - deltaVB)*0.5*pcell->ilen[j];
      Fp += mt->band->SchottyJsp(pi,Ti,pbc->WorkFunction-aux[i].affinity + deltaVB)*0.5*pcell->ilen[j];
      PetscScalar h = pbc->heat_transfer;
      PetscScalar r = h*pbc->T_external;
      grad_T += 0.5*pcell->ilen[j]*(r-0.25*h*(3*Ti+Tj));
    }
  }

  f[zofs[zone_index]+4*i+1] = (f[zofs[zone_index]+4*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+4*i+2] = (f[zofs[zone_index]+4*i+2]-Fp)/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T)/pcell->area;

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }
}


void SMCZone::F2E_mix_ddm_stkbc_interface(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int>& zofs,DABC &bc,
    ElZone *pz, int n)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);

  //Schotty Barrier Lowerring
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  //potential
  f[zofs[z]+4*i+0] = x[zofs[z]+4*i+0] + pbc->WorkFunction - deltaVB - pbc->Vapp;

  PetscScalar Fn=0,Fp=0,grad_T=0;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int  nb = pcell->nb_array[j];
    if(j==0||j==pcell->nb_num-1)
    {
      Fn += mt->band->SchottyJsn(ni,Ti,pbc->WorkFunction-aux[i].affinity - deltaVB)*0.5*pcell->ilen[j];
      Fp += mt->band->SchottyJsp(pi,Ti,pbc->WorkFunction-aux[i].affinity + deltaVB)*0.5*pcell->ilen[j];
    }
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
  f[zofs[zone_index]+4*i+1] = (f[zofs[zone_index]+4*i+1]+Fn)/pcell->area;
  f[zofs[zone_index]+4*i+2] = (f[zofs[zone_index]+4*i+2]-Fp)/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T);

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*Tt/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*Tt/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity1*(Ti-fs[i].T)/ODE_F.dt*pcell->area
                                   -pz->aux[n].density*HeatCapacity2*(Ti-fs[i].T)/ODE_F.dt*ncell->area;
    }
  }
}


void SMCZone::F2E_mix_ddm_insulator_gate(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));

  int size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e  = mt->e;
  PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[zone_index]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[zone_index]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[zone_index]+4*i+3];      //lattice temperature of node i
  mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar grad_P = 0, grad_T = 0;

  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
    PetscScalar Tj = x[zofs[zone_index]+4*nb+3];            //lattice temperature of nb node

    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar vgate = pbc->Vapp - pbc->WorkFunction;
      PetscScalar q = e*pbc->QF; //sigma is the surface change density
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar eps = aux[i].eps;
      PetscScalar r=q/eps + eps_ox/eps/Thick*vgate;
      PetscScalar s=eps_ox/eps/Thick;
      grad_P += 0.5*pcell->ilen[j]*(r-0.25*s*(3*Vi+Vj));

      PetscScalar h = pbc->heat_transfer;
      grad_T += 0.5*pcell->ilen[j]*(h*pbc->T_external-0.25*h*(3*Ti+Tj));
    }
  }

  f[zofs[zone_index]+4*i+0] = (f[zofs[zone_index]+4*i+0]+grad_P)/pcell->area
                              + e/aux[i].eps*((pi-aux[i].Na)+(aux[i].Nd-ni));
  f[zofs[zone_index]+4*i+1] = f[zofs[zone_index]+4*i+1]/pcell->area;
  f[zofs[zone_index]+4*i+2] = f[zofs[zone_index]+4*i+2]/pcell->area;
  f[zofs[zone_index]+4*i+3] = (f[zofs[zone_index]+4*i+3]+grad_T)/pcell->area;


  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar Tn = (2-r)/(1-r)*ni-1.0/(r*(1-r))*fs[i].n+(1-r)/r*fs[i].n_last;
      PetscScalar Tp = (2-r)/(1-r)*pi-1.0/(r*(1-r))*fs[i].p+(1-r)/r*fs[i].p_last;
      PetscScalar Tt = (2-r)/(1-r)*Ti-1.0/(r*(1-r))*fs[i].T+(1-r)/r*fs[i].T_last;
      f[zofs[zone_index]+4*i+1] += -Tn/(ODE_F.dt_last+ODE_F.dt);
      f[zofs[zone_index]+4*i+2] += -Tp/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*Tt/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      f[zofs[zone_index]+4*i+1] += -(ni-fs[i].n)/ODE_F.dt;
      f[zofs[zone_index]+4*i+2] += -(pi-fs[i].p)/ODE_F.dt;
      PetscScalar HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      f[zofs[zone_index]+4*i+3] += -aux[i].density*HeatCapacity*(Ti-fs[i].T)/ODE_F.dt;
    }
  }
}


//------------------------------------------------------------------------------------------

void SMCZone::J2E_mix_ddm_ombc_segment(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i
  PetscScalar Vt = kb*Ti/e;
  PetscScalar area = pcell->area;
  PetscScalar d_grad_T_dTi = 0;
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    PetscScalar kapa =  mt->thermal->HeatConduction(0.5*(Ti+Tj));
    //-------------------------------------
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j]/area;
    //-------------------------------------
    PetscScalar d_grad_T_dTj =  kapa*pcell->elen[j]/pcell->ilen[j]/area; //dfun(3)/dT(r)
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/area;
      d_grad_T_dTj += -0.5*pcell->ilen[j]*0.25*h/area;
    }
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[z]+4*nb+3,d_grad_T_dTj,INSERT_VALUES);
  }
  if(Fermi)  //Fermi
  {
    PetscScalar Vt =  kb*fs[i].T/e;
    mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
    PetscScalar Nc  = mt->band->Nc(fs[i].T);
    PetscScalar Nv  = mt->band->Nv(fs[i].T);
    PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) + kb*fs[i].T*log(aux[i].Nc));
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) - kb*fs[i].T*log(aux[i].Nv)+ aux[i].Eg);
    PetscScalar phin = pbc->Vapp;
    PetscScalar phip = pbc->Vapp;
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    PetscScalar detan_dVi =  1.0/Vt;
    PetscScalar detap_dVi = -1.0/Vt;
    Set_Mat4_zero(A);
    A.m[0] = Nc*fermi_mhalf(etan)*detan_dVi - Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[4] = -Nc*fermi_mhalf(etan)*detan_dVi;  
    A.m[5] = 1;  
    A.m[8] = -Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[10] = 1;
    A.m[15] =  d_grad_T_dTi;
    MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
  }
  else  //Boltzmann
  {
    Set_Mat4_I(A);
    A.m[15] =  d_grad_T_dTi;     //dfun(3)/dT(i)
  }
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }

  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}


void SMCZone::J2E_mix_ddm_ombc_interface(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc,
    ElZone *pz, int n)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  OhmicBC *pbc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i
  PetscScalar Vt = kb*Ti/e;
  PetscScalar area = pcell->area;
  PetscScalar d_grad_T_dTi = 0;
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    PetscScalar kapa =  mt->thermal->HeatConduction(0.5*(Ti+Tj));
    //-------------------------------------
    d_grad_T_dTi += -kapa*pcell->elen[j]/pcell->ilen[j];
    //-------------------------------------
    PetscScalar d_grad_T_dTj =  kapa*pcell->elen[j]/pcell->ilen[j]; //dfun(3)/dT(r)
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[z]+4*nb+3,d_grad_T_dTj,INSERT_VALUES);
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
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[pz->pzone->zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }

  if(Fermi)  //Fermi
  {
    PetscScalar Vt =  kb*fs[i].T/e;
    mt->mapping(&pzone->danode[i],&aux[i],ODE_F.clock);
    PetscScalar Nc  = mt->band->Nc(fs[i].T);
    PetscScalar Nv  = mt->band->Nv(fs[i].T);
    PetscScalar Vi = x[zofs[zone_index]+4*i+0];     //potential of node i
    PetscScalar Ec =  -(e*Vi + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) + kb*fs[i].T*log(aux[i].Nc));
    PetscScalar Ev =  -(e*Vi + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) - kb*fs[i].T*log(aux[i].Nv)+ aux[i].Eg);
    PetscScalar phin = pbc->Vapp;
    PetscScalar phip = pbc->Vapp;
    PetscScalar etan = (-e*phin-Ec)/kb/fs[i].T;
    PetscScalar etap = (Ev+e*phip)/kb/fs[i].T;
    PetscScalar detan_dVi =  1.0/Vt;
    PetscScalar detap_dVi = -1.0/Vt;
    Set_Mat4_zero(A);
    A.m[0] = Nc*fermi_mhalf(etan)*detan_dVi - Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[4] = -Nc*fermi_mhalf(etan)*detan_dVi;  
    A.m[5] = 1;  
    A.m[8] = -Nv*fermi_mhalf(etap)*detap_dVi;  
    A.m[10] = 1;
    A.m[15] =  d_grad_T_dTi;
    MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
  }
  else  //Boltzmann
  {
    Set_Mat4_I(A);
    A.m[15] =  d_grad_T_dTi;     //dfun(3)/dT(i)
  }
  
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*pcell->area
               -pz->aux[n].density*HeatCapacity2*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1/ODE_F.dt*pcell->area
               -pz->aux[n].density*HeatCapacity2/ODE_F.dt*ncell->area;
    }
  }

  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}


void SMCZone::J2E_mix_ddm_stkbc_segment(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc)
{
  Mat4   A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int size = pzone->davcell.size();

  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  PetscScalar area = pcell->area;
  PetscScalar d_grad_T_dTi = 0;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fn_dTi = 0;
  PetscScalar d_Fp_dpi = 0;
  PetscScalar d_Fp_dTi = 0;
  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar nj = x[zofs[z]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    MatGetValues(*jtmp,4,index,4,col,A.m);
    A.m[0] = 0.0;
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn (ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fn_dTi += mt->band->pdSchottyJsn_pdTl(ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp (pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dTi += mt->band->pdSchottyJsp_pdTl(pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/pcell->area;
      A.m[15] += -0.5*pcell->ilen[j]*0.25*h/pcell->area;
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);
  }
  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/area;

  A.m[0] =  1;

  A.m[5] +=  d_Fn_dni;                              //dfun(1)/dn(i)
  A.m[7] +=  d_Fn_dTi;

  A.m[10]+=  -d_Fp_dpi;                                       //dfun(2)/dp(i)
  A.m[11]+=  -d_Fp_dTi;

  A.m[15]+=  d_grad_T_dTi; //dfun(3)/dT(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}


void SMCZone::J2E_mix_ddm_stkbc_interface(int i,PetscScalar *x,Mat *jac,Mat *jtmp,ODE_Formula &ODE_F, vector<int> &zofs,DABC &bc,
    ElZone *pz, int n)
{
  Mat4   A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  SchottkyBC *pbc = dynamic_cast<SchottkyBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int equ_num = 4;
  int size = pzone->davcell.size();
  PetscScalar area = pcell->area;
  Vec4 scale;
  Set_Vec4(scale,area,area,area,1.0);
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar VB = pbc->WorkFunction-aux[i].affinity;
  PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[i].eps,sqrt(aux[i].Ex*aux[i].Ex+aux[i].Ey*aux[i].Ey));
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  PetscScalar d_grad_T_dTi = 0;
  PetscScalar d_Fn_dni = 0;
  PetscScalar d_Fn_dTi = 0;
  PetscScalar d_Fp_dpi = 0;
  PetscScalar d_Fp_dTi = 0;
  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar nj = x[zofs[z]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    MatGetValues(*jtmp,4,index,4,col,A.m);
    A.m[0] = 0.0;
    A=A/scale;
    if(j==0||j==pcell->nb_num-1)
    {
      d_Fn_dni += mt->band->pdSchottyJsn_pdn (ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fn_dTi += mt->band->pdSchottyJsn_pdTl(ni,Ti,VB-deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dpi += mt->band->pdSchottyJsp_pdp (pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
      d_Fp_dTi += mt->band->pdSchottyJsp_pdTl(pi,Ti,VB+deltaVB)*0.5*pcell->ilen[j]/pcell->area;
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);
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
    MatSetValue(*jac,zofs[z]+4*i+3,zofs[pz->pzone->zone_index]+2*nb+1,d_grad_T_dTj,INSERT_VALUES);
  }
  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole
  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/scale;

  A.m[0] =  1;

  A.m[5] +=  d_Fn_dni;                              //dfun(1)/dn(i)
  A.m[7] +=  d_Fn_dTi;

  A.m[10]+=  -d_Fp_dpi;                                       //dfun(2)/dp(i)
  A.m[11]+=  -d_Fp_dTi;

  A.m[15]+=  d_grad_T_dTi; //dfun(3)/dT(i)

  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*pcell->area
                 -pz->aux[n].density*HeatCapacity2*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt)*ncell->area;
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar HeatCapacity1 =  mt->thermal->HeatCapacity(Ti);
      PetscScalar HeatCapacity2 =  pz->mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity1/ODE_F.dt*pcell->area
                 -pz->aux[n].density*HeatCapacity2/ODE_F.dt*ncell->area;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}


void SMCZone::J2E_mix_ddm_insulator_gate(int i,PetscScalar *x, Mat *jac, Mat *jtmp, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc)
{
  Mat4    A;
  PetscInt       index[4],col[4];
  const VoronoiCell* pcell = pzone->davcell.GetPointer(i);
  InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC * >(bc.Get_pointer(pcell->bc_index-1));
  int z = zone_index;
  int  size = pzone->davcell.size();
  PetscScalar kb = mt->kb;
  PetscScalar e = mt->e;
  PetscScalar Vi = x[zofs[z]+4*i+0];     //potential of node i
  PetscScalar ni = x[zofs[z]+4*i+1];     //electron density of node i
  PetscScalar pi = x[zofs[z]+4*i+2];     //hole density of node i
  PetscScalar Ti = x[zofs[z]+4*i+3];   //lattice temperature of node i

  PetscScalar area = pcell->area;
  PetscScalar L = 0.5*(pcell->ilen[0]+pcell->ilen[pcell->nb_num-1]);
  PetscScalar d_grad_P_dVi = 0;
  PetscScalar d_grad_P_dVapp = 0;
  PetscScalar d_grad_T_dTi = 0;

  //--------------------------------
  index[0] = zofs[z]+4*i+0;
  index[1] = zofs[z]+4*i+1;
  index[2] = zofs[z]+4*i+2;
  index[3] = zofs[z]+4*i+3;
  //--------------------------------
  for(int j=0;j<pcell->nb_num;j++)
  {
    int    nb = pcell->nb_array[j];
    PetscScalar Vj = x[zofs[z]+4*nb+0];     //potential of nb node
    PetscScalar nj = x[zofs[z]+4*nb+1];     //electron density of nb node
    PetscScalar pj = x[zofs[z]+4*nb+2];     //hole density of nb node
    PetscScalar Tj = x[zofs[z]+4*nb+3];     //lattice temperature of nb node
    col[0] = zofs[z]+4*nb+0;
    col[1] = zofs[z]+4*nb+1;
    col[2] = zofs[z]+4*nb+2;
    col[3] = zofs[z]+4*nb+3;
    MatGetValues(*jtmp,4,index,4,col,A.m);
    A=A/area;
    if(j==0||j==pcell->nb_num-1)
    {
      PetscScalar Thick = pbc->Thick;
      PetscScalar eps_ox = mt->eps0*pbc->eps;
      PetscScalar eps = aux[i].eps;
      PetscScalar s=eps_ox/eps/Thick;
      d_grad_P_dVi += -0.5*pcell->ilen[j]*0.25*s*3/area;
      d_grad_P_dVapp += 0.5*pcell->ilen[j]*s/area;
      A.m[0]+= -0.5*pcell->ilen[j]*0.25*s/area;
      PetscScalar h = pbc->heat_transfer;
      d_grad_T_dTi += -0.5*pcell->ilen[j]*0.25*h*3/area;
      A.m[15] += -0.5*pcell->ilen[j]*0.25*h/area;
    }
    MatSetValues(*jac,4,index,4,col,A.m,INSERT_VALUES);

  }
  //fun(0) is the poisson's equation
  //fun(1) is the continuous equation of electron
  //fun(2) is the continuous equation of hole

  MatGetValues(*jtmp,4,index,4,index,A.m);
  A=A/area;
  A.m[0] +=  d_grad_P_dVi;
  A.m[1] +=  -e/aux[i].eps;                                //dfun(0)/dn(i)
  A.m[2] +=   e/aux[i].eps;                                //dfun(0)/dp(i)

  A.m[15]+=  d_grad_T_dTi; //dfun(3)/dT(i)
  if(ODE_F.TimeDependent==true)
  {
    //second order
    if(ODE_F.BDF_Type==BDF2&&ODE_F.BDF2_restart==false)
    {
      PetscScalar r = ODE_F.dt_last/(ODE_F.dt_last+ODE_F.dt);
      A.m[5] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      A.m[10] += -(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity*(2-r)/(1-r)/(ODE_F.dt_last+ODE_F.dt);
    }
    else //first order
    {
      A.m[5] += -1/ODE_F.dt;
      A.m[10] += -1/ODE_F.dt;
      PetscScalar       HeatCapacity =  mt->thermal->HeatCapacity(Ti);
      A.m[15] += -aux[i].density*HeatCapacity/ODE_F.dt;
    }
  }
  MatSetValues(*jac,4,index,4,index,A.m,INSERT_VALUES);
}


void SMCZone::F2E_mix_om_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
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
    PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //potential of node i
    current += DeviceDepth*(f[zofs[zone_index]+4*node+1]-f[zofs[zone_index]+4*node+2]);
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  pbc->Set_Current_new(current);

}


void SMCZone::F2E_mix_stk_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int bc_index = electrode[i];
  SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc.Get_pointer(bc_index));

  PetscScalar current=0;
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //electron density of node i
    PetscScalar ni = x[zofs[zone_index]+4*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+4*node+2];     //hole density of node i
    PetscScalar Ti = x[zofs[zone_index]+4*node+3];
    mt->mapping(&pzone->danode[node],&aux[node],ODE_F.clock);
    PetscScalar deltaVB=mt->band->SchottyBarrierLowerring(aux[node].eps,sqrt(aux[node].Ex*aux[node].Ex+aux[node].Ey*aux[node].Ey));
    PetscScalar Jsn = mt->band->SchottyJsn(ni,Ti,pbc->WorkFunction-aux[node].affinity-deltaVB);
    PetscScalar Jsp = mt->band->SchottyJsp(pi,Ti,pbc->WorkFunction-aux[node].affinity+deltaVB);
    current += -Jsn*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[0]*DeviceDepth;
    current += -Jsn*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    current += -Jsp*0.5*pcell->ilen[pcell->nb_num-1]*DeviceDepth;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  pbc->Set_Current_new(current);

}


void SMCZone::F2E_mix_ins_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
  int bc_index = electrode[i];
  InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc.Get_pointer(bc_index));
  PetscScalar current=0;

  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node = bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //electron density of node i
    for(int k=0;k<pcell->nb_num;k++)
    {
      int  nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
      //displacement current
      current += DeviceDepth*pcell->elen[k]*aux[node].eps*((Vi-Vj)-(fs[node].P-fs[nb].P))/pcell->ilen[k]/ODE_F.dt;
    }
  }
  pbc->Set_Current_new(current);
}


void SMCZone::F2E_mix_om_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,PetscScalar *x, Mat *jtmp, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
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
  PetscScalar    A1[4],A2[4];
  PetscInt       index[4],col[4];
  for(int j=0;j<bc[bc_index].psegment->node_array.size();j++)
  {
    int node=bc[bc_index].psegment->node_array[j];
    const VoronoiCell* pcell = pzone->davcell.GetPointer(node);
    PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //potential of node i
    PetscScalar dJdisp_dVi=0,dJdisp_dVr=0;
    index[0] = zofs[zone_index]+4*node+0;
    index[1] = zofs[zone_index]+4*node+1;
    index[2] = zofs[zone_index]+4*node+2;
    index[3] = zofs[zone_index]+4*node+3;
    for(int k=0;k<pcell->nb_num;k++)
    {
      int    nb = pcell->nb_array[k];
      PetscScalar Vj = x[zofs[zone_index]+4*nb+0];     //potential of nb node
      col[0] = zofs[zone_index]+4*nb+0;
      col[1] = zofs[zone_index]+4*nb+1;
      col[2] = zofs[zone_index]+4*nb+2;
      col[3] = zofs[zone_index]+4*nb+3;

      MatGetValues(*jtmp,1,&index[1],4,col,A1);
      MatGetValues(*jtmp,1,&index[2],4,col,A2);

      //for displacement current
      dJdisp_dVi += aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];
      dJdisp_dVr = -aux[node].eps/pcell->ilen[k]/ODE_F.dt*pcell->elen[k];

      //
      apdI_pdw[col[0]] += (A1[0]-A2[0]+dJdisp_dVr)*DeviceDepth;
      apdI_pdw[col[1]] += (A1[1]-A2[1])*DeviceDepth;
      apdI_pdw[col[2]] += (A1[2]-A2[2])*DeviceDepth;
      apdI_pdw[col[3]] += (A1[3]-A2[3])*DeviceDepth;
    }

    MatGetValues(*jtmp,1,&index[1],4,index,A1);
    MatGetValues(*jtmp,1,&index[2],4,index,A2);

    apdI_pdw[index[0]] += (A1[0]-A2[0]+dJdisp_dVi)*DeviceDepth;
    apdI_pdw[index[1]] += (A1[1]-A2[1])*DeviceDepth;
    apdI_pdw[index[2]] += (A1[2]-A2[2])*DeviceDepth;
    apdI_pdw[index[3]] += (A1[3]-A2[3])*DeviceDepth;

    if(Fermi)
    {
      PetscScalar Vt =  kb*fs[node].T/e;
      mt->mapping(&pzone->danode[node],&aux[node],ODE_F.clock);
      PetscScalar Nc  = mt->band->Nc(fs[node].T);
      PetscScalar Nv  = mt->band->Nv(fs[node].T);
      PetscScalar Vi = x[zofs[zone_index]+4*node+0];     //potential of node 
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
    else  //Boltzmann
      apdF_pdV[index[0]] = 1.0;
    pdI_pdV = 0;
  }
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F2E_mix_stk_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,PetscScalar *x, Mat *jtmp, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
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
    PetscScalar ni = x[zofs[zone_index]+4*node+1];     //electron density of node i
    PetscScalar pi = x[zofs[zone_index]+4*node+2];     //hole density of node i
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
      apdI_pdw[zofs[zone_index]+4*node+0] += dI_dVi;
      apdI_pdw[zofs[zone_index]+4*nb+0]   += dI_dVr;
    }
    apdI_pdw[zofs[zone_index]+4*node+1] = DeviceDepth*dJ_dni;
    apdI_pdw[zofs[zone_index]+4*node+2] = DeviceDepth*dJ_dpi;

    apdF_pdV[zofs[zone_index]+4*node+0] = 1.0;
    pdI_pdV = 0;
  }

  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


void SMCZone::F2E_mix_ins_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,PetscScalar *x, Mat *jtmp, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
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
      apdI_pdw[zofs[zone_index]+4*node+0] += dI_dVi;
      apdI_pdw[zofs[zone_index]+4*nb+0]   += dI_dVr;
      if(k==0||k==pcell->nb_num-1)
      {
        PetscScalar Thick = pbc->Thick;
        PetscScalar eps_ox = mt->eps0*pbc->eps;
        PetscScalar s=eps_ox/Thick;
        apdF_pdV[zofs[zone_index]+4*node+0] += 0.5*pcell->ilen[k]*s/pcell->area;
      }
    }
    pdI_pdV = 0;
  }

  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}

