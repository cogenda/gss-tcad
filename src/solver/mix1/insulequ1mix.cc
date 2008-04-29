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

//-----------------------------------------------------------------------------
// special process of electrode bcs for mix type simulation
//-----------------------------------------------------------------------------

void ISZone::F1_mix_ddm_gatebc(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  int gate_equ;
  const VoronoiCell* pcell = &pzone->davcell[i];
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {gate_equ=j;break;}
  }
  GateBC *pbc = dynamic_cast<GateBC * >(bc.Get_pointer(pcell->bc_index-1));
  f[zofs[zone_index]+i] = x[zofs[zone_index]+i] - pbc->Vapp + pbc->WorkFunction;
  //just for fill the equ entry 
  f[zofs[zone_index]+equ_num*size+gate_equ]=0.0;
}


void ISZone::F1_mix_gate_electrode_current(int i,PetscScalar *x,PetscScalar *f, ODE_Formula &ODE_F, vector<int> & zofs, DABC &bc, PetscScalar DeviceDepth)
{
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
  
  pbc->Set_Current_new(current);
}


void ISZone::J1_mix_ddm_gatebc(int i,PetscScalar *x,Mat *jac, ODE_Formula &ODE_F, vector<int> & zofs,DABC &bc )
{
  int equ_num = 1;
  int size = pzone->davcell.size();
  int gate_equ;
  const VoronoiCell* pcell = &pzone->davcell[i];
  for(int j=0;j<electrode.size();j++)
  {
    if(electrode[j]==pcell->bc_index-1) {gate_equ=j;break;}
  }
  MatSetValue(*jac,zofs[zone_index]+i,zofs[zone_index]+i,1.0,INSERT_VALUES);
  //just for fill the matrix entry 
  MatSetValue(*jac,zofs[zone_index]+equ_num*size+gate_equ,zofs[zone_index]+equ_num*size+gate_equ,1.0,INSERT_VALUES);
}


void ISZone::F1_mix_gate_electrode_Load(int bc_index, double & current, PetscScalar & pdI_pdV, Vec & pdI_pdw, Vec & pdF_pdV,
         PetscScalar *x, Mat *jac, ODE_Formula &ODE_F, vector<int> &zofs, DABC & bc, PetscScalar DeviceDepth)
{
  GateBC *pbc = dynamic_cast <GateBC * > (bc.Get_pointer(bc_index));
  pdI_pdV = 0.0;
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
	apdI_pdw[zofs[zone_index]+node] += dJ_dVi*DeviceDepth;
	apdI_pdw[zofs[zone_index]+nb]   += dJ_dVr*DeviceDepth;
      }
      apdF_pdV[zofs[zone_index]+node] = 1.0;
  }
  VecRestoreArray(pdI_pdw,&apdI_pdw);
  VecRestoreArray(pdF_pdV,&apdF_pdV);
}


