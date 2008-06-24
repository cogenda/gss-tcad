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
/*  Last update: Jan 19, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#include "ddm_nt1e.h"
#include "private/kspimpl.h"
#include "private/snesimpl.h"
#include "src/snes/impls/tr/tr.h"
#include "log.h"
#include "opt_update.h"
#include "gen_update.h"



/* ----------------------------------------------------------------------------
 * form_function_pn:  This function setup DDM equation F(x)=0
 */
void DDM_Solver_L1E::form_function_pn_1E(PetscScalar *x,PetscScalar *f)
{
  // compute flux along triangle edges. semiconductor zone only
  for(int z=0;z<zone_num;z++)
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].datri.size();i++)
      {
        Tri *ptri = &zone[z].datri[i];
        pzonedata->F1E_Tri_ddm(ptri,x,f,zofs);
      }
    }

  for(int z=0;z<zone_num;z++)
  {
    //semiconductor zone
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[z]);

      // process electrode.
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==OhmicContact)
          pzonedata->F1E_om_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==SchottkyContact)
          pzonedata->F1E_stk_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==InsulatorContact)
          pzonedata->F1E_ins_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
      }

      // process cell variables and boundaries.
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index||bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F1E_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
          pzonedata->F1E_ddm_ombc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
          pzonedata->F1E_ddm_stkbc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorContact)
          pzonedata->F1E_ddm_insulator_gate(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F1E_ddm_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HomoInterface)
        {
          HomoInterfaceBC *pbc;
          pbc = dynamic_cast<HomoInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F1E_ddm_homojunction(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HeteroInterface)
        {
          HeteroInterfaceBC *pbc;
          pbc = dynamic_cast<HeteroInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F1E_ddm_heterojunction(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->F1E_ddm_inner(i,x,f,ODE_formula,zofs);
      }
    }

    // Insulator zone
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata= dynamic_cast< ISZone * >(zonedata[z]);
      pzonedata->F1_efield_update(x,zofs,bc,zonedata);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index || bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F1_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==GateContact)
          pzonedata->F1_ddm_gatebc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
          pzonedata->F1_ddm_chargebc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
        {
          ElectrodeInsulatorBC *pbc;
          pbc = dynamic_cast<ElectrodeInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
          pzonedata->F1_ddm_electrode_insulator_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F1_ddm_semiconductor_insulator_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Insulator_Insulator)
        {
          InsulatorInsulatorBC *pbc;
          pbc = dynamic_cast<InsulatorInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F1_ddm_insulator_insulator_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->F1_ddm_inner(i,x,f,ODE_formula,zofs);
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==GateContact)
          pzonedata->F1_gate_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          pzonedata->F1_charge_electrode(i,x,f,ODE_formula,zofs,bc);
      }
    }

    // Electrode zone
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata= dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index || bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F1_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
          pzonedata->F1_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==OhmicContact)
        {
          OhmicBC *pbc = dynamic_cast<OhmicBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F1_ddm_om_contact(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==SchottkyContact)
        {
          SchottkyBC *pbc = dynamic_cast<SchottkyBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F1_ddm_om_contact(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==GateContact)
        {
          GateBC *pbc = dynamic_cast<GateBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F1_ddm_gate_contact(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
        {
          ChargedContactBC *pbc = dynamic_cast<ChargedContactBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F1_ddm_charge_contact(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->F1_ddm_inner(i,x,f,ODE_formula,zofs);
      }
    }

  }
}


/* ----------------------------------------------------------------------------
 * form_jacobian_pn:  This function setup jacobian matrix F'(x) of DDM equation F(x)=0
 * dual carrier edition
 */
void DDM_Solver_L1E::form_jacobian_pn_1E(PetscScalar *x, Mat *jac, Mat *jtmp)
{

  //Compute Jacobian entries for flux along triangle edges.
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[z]);

      // compute flux Jacobian along triangle edges
      for(int i=0;i<zone[z].datri.size();i++)
      {
        Tri *ptri = &zone[z].datri[i];
        pzonedata->J1E_Tri_ddm(ptri,x,jtmp,zofs);
      }
    }
  }
  MatAssemblyBegin(*jtmp,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jtmp,MAT_FINAL_ASSEMBLY);

  //Compute Jacobian entries and insert into matrix.
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[z]);

      // process electrode.
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==OhmicContact)
          pzonedata->J1E_om_electrode(i,x,jac,&JTmp,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==SchottkyContact)
          pzonedata->J1E_stk_electrode(i,x,jac,&JTmp,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==InsulatorContact)
          pzonedata->J1E_ins_electrode(i,x,jac,&JTmp,ODE_formula,zofs,bc,DeviceDepth);
      }

      // process cell variables and boundaries.
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        SemiData *pfs = &pzonedata->fs[i];
        if(!pcell->bc_index||bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->J1E_ddm_inner(i,x,jac,&JTmp,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
          pzonedata->J1E_ddm_ombc(i,x,jac,&JTmp,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
          pzonedata->J1E_ddm_stkbc(i,x,jac,&JTmp,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorContact)
          pzonedata->J1E_ddm_insulator_gate(i,x,jac,&JTmp,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J1E_ddm_interface(i,x,jac,&JTmp,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HomoInterface)
        {
          HomoInterfaceBC *pbc;
          pbc = dynamic_cast<HomoInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J1E_ddm_homojunction(i,x,jac,&JTmp,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HeteroInterface)
        {
          HeteroInterfaceBC *pbc;
          pbc = dynamic_cast<HeteroInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J1E_ddm_heterojunction(i,x,jac,&JTmp,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->J1E_ddm_inner(i,x,jac,&JTmp,ODE_formula,zofs);
      }
    }

    // Insulator zone
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata= dynamic_cast< ISZone * >(zonedata[z]);

      for(int i=0;i<zone[z].davcell.size();i++)
      {

        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index||bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->J1_ddm_inner(i,x,jac,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==GateContact)
          pzonedata->J1_ddm_gatebc(i,x,jac,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==ChargedContact)
          pzonedata->J1_ddm_chargebc(i,x,jac,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
        {
          ElectrodeInsulatorBC *pbc;
          pbc = dynamic_cast<ElectrodeInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
          pzonedata->J1_ddm_electrode_insulator_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J1_ddm_semiconductor_insulator_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Insulator_Insulator)
        {
          InsulatorInsulatorBC *pbc;
          pbc = dynamic_cast<InsulatorInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J1_ddm_insulator_insulator_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->J1_ddm_inner(i,x,jac,ODE_formula,zofs);
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==GateContact)
          pzonedata->J1_gate_electrode(i,x,jac,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          pzonedata->J1_charge_electrode(i,x,jac,ODE_formula,zofs,bc);
      }
    }

    // Electrode zone
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata= dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index || bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->J1_ddm_inner(i,x,jac,ODE_formula,zofs);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
          pzonedata->J1_ddm_inner(i,x,jac,ODE_formula,zofs);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==OhmicContact)
        {
          OhmicBC *pbc = dynamic_cast<OhmicBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J1_ddm_om_contact(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==SchottkyContact)
        {
          SchottkyBC *pbc = dynamic_cast<SchottkyBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J1_ddm_stk_contact(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==GateContact)
        {
          GateBC *pbc = dynamic_cast<GateBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J1_ddm_gate_contact(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
        {
          ChargedContactBC *pbc = dynamic_cast<ChargedContactBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J1_ddm_charge_contact(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->J1_ddm_inner(i,x,jac,ODE_formula,zofs);
      }
    }

  }
}




/* ----------------------------------------------------------------------------
 * SNES_form_function_pn_1:  wrap function for petsc nonlinear solver
 */
PetscErrorCode SNES_form_function_pn_1E(SNES snes, Vec x,Vec f,void *dummy)
{
  PetscScalar    *xx,*ff;
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)dummy;
  VecZeroEntries(f);
  //Get pointers to vector data.
  VecGetArray(x,&xx);
  VecGetArray(f,&ff);

  ps->form_function_pn_1E(xx,ff);
  ps->error_norm_estimat(xx,ff);

  //Restore vectors
  VecRestoreArray(x,&xx);
  VecRestoreArray(f,&ff);
  VecNorm(f,NORM_2,&ps->norm);
  return 0;
}


/* ----------------------------------------------------------------------------
 * SNES_form_jacobian_pn_1:  wrap function for petsc nonlinear solver
 */
PetscErrorCode SNES_form_jacobian_pn_1E(SNES snes, Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
  PetscScalar    *xx;
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)dummy;
  //Get pointer to vector data
  VecGetArray(x,&xx);
  //clear old matrix
  MatZeroEntries(*jac);
  MatZeroEntries(ps->JTmp);
  //build matrix here
  ps->form_jacobian_pn_1E(xx,jac,&ps->JTmp);

  *flag = SAME_NONZERO_PATTERN;
  //Restore vector
  VecRestoreArray(x,&xx);
  //Assemble matrix
  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  //MatView(*jac,PETSC_VIEWER_DRAW_WORLD);
  return 0;
}

/* ----------------------------------------------------------------------------
 * SNESMF_form_jacobian_pn_1:  wrap function for petsc nonlinear solver with
 * matrix-free method
 */
PetscErrorCode SNESMF_form_jacobian_pn_1E(SNES snes, Vec x,Mat *A,Mat *B,MatStructure *flag,void *dummy)
{
  PetscScalar    *xx;
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)dummy;
  //Get pointer to vector data
  VecGetArray(x,&xx);
  //clear old matrix
  MatZeroEntries(*B);
  MatZeroEntries(ps->JTmp);
  //build matrix here
  ps->form_jacobian_pn_1E(xx,B,&ps->JTmp);
  *flag = SAME_NONZERO_PATTERN;
  //Restore vector
  VecRestoreArray(x,&xx);
  //Assemble matrix
  MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
  //MatView(*B,PETSC_VIEWER_DRAW_WORLD);
  return 0;
}


/* ----------------------------------------------------------------------------
 * error_norm_estimat:  This function compute X and RHS error norms. 
 */
void DDM_Solver_L1E:: error_norm_estimat(PetscScalar *x,PetscScalar *f)
{
  // do clear
  potential_norm = 0;
  electron_norm  = 0;
  hole_norm      = 0;

  possion_norm        = 0;
  elec_continuty_norm = 0;
  hole_continuty_norm = 0;
  electrode_norm      = 0;

  int offset=0;
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        potential_norm += x[offset+0]*x[offset+0];
        electron_norm  += x[offset+1]*x[offset+1];
        hole_norm      += x[offset+2]*x[offset+2];
        possion_norm   += f[offset+0]*f[offset+0];
        elec_continuty_norm += f[offset+1]*f[offset+1];
        hole_continuty_norm += f[offset+2]*f[offset+2];
        offset += 3;
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        potential_norm += x[offset]*x[offset];
        electrode_norm += f[offset]*f[offset];
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        potential_norm += x[offset]*x[offset];
        possion_norm   += f[offset]*f[offset];
        offset += 1;
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        potential_norm += x[offset]*x[offset];
        electrode_norm += f[offset]*f[offset];
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        potential_norm += x[offset]*x[offset];
        possion_norm   += f[offset]*f[offset];
        offset += 1;
      }
    }
  }

  potential_norm = sqrt(potential_norm);
  electron_norm  = sqrt(electron_norm);
  hole_norm      = sqrt(hole_norm);

  possion_norm        = sqrt(possion_norm);
  elec_continuty_norm = sqrt(elec_continuty_norm);
  hole_continuty_norm = sqrt(hole_continuty_norm);
  electrode_norm      = sqrt(electrode_norm);

}


PetscErrorCode KSPConvergenceTest_L1E(KSP ksp,PetscInt n,PetscReal rnorm,KSPConvergedReason* reason,void *dummy)
{
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)dummy;
  ksp->rtol = 1e-12*ps->N;
  ksp->abstol= PetscMax(1e-6*ps->norm,1e-15*ps->N);
  return KSPDefaultConverged(ksp,n,rnorm,reason,0);
}


PetscErrorCode LineSearch_ConvergenceTest_L1E(SNES snes,PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm,
    SNESConvergedReason *reason,void *dummy)
{
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)dummy;

  *reason = SNES_CONVERGED_ITERATING;
  if (!it)
  {
    snes->ttol = fnorm*snes->rtol;
    gss_log.string_buf()<<" "<<"its\t"<<"| Eq(V) |   "<<"| Eq(n) |   "<<"| Eq(p) |   "<<"|delta x|\n"
    <<"----------------------------------------------------------------------\n";
    gss_log.record();
  }

  gss_log.string_buf().precision(3);
  gss_log.string_buf()<<"  "<<it<<"\t"
  <<ps->possion_norm<<"   "
  <<ps->elec_continuty_norm<<"   "
  <<ps->hole_continuty_norm<<"   "
  <<pnorm<<"\n";
  gss_log.record();
  gss_log.string_buf().precision(6);

  if (fnorm != fnorm || 
      ps->possion_norm != ps->possion_norm || 
      ps->elec_continuty_norm != ps->elec_continuty_norm ||
      ps->hole_continuty_norm != ps->hole_continuty_norm)
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  else if (ps->possion_norm < ps->possion_abs_toler &&
           ps->elec_continuty_norm < ps->elec_continuty_abs_toler &&
           ps->hole_continuty_norm < ps->hole_continuty_abs_toler &&
           ps->electrode_norm      < ps->electrode_abs_toler)
  {
    *reason = SNES_CONVERGED_FNORM_ABS;
  }
  else if (snes->nfuncs >= snes->max_funcs)
  {
    *reason = SNES_DIVERGED_FUNCTION_COUNT;
  }
    
  if (it && !*reason)
  {
    if (fnorm <= snes->ttol)
    {
      *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }
    else if (pnorm < ps->relative_toler &&
             ps->possion_norm        < ps->toler_relax*ps->possion_abs_toler &&
             ps->elec_continuty_norm < ps->toler_relax*ps->elec_continuty_abs_toler &&
             ps->hole_continuty_norm < ps->toler_relax*ps->hole_continuty_abs_toler &&
             ps->electrode_norm      < ps->toler_relax*ps->electrode_abs_toler)
    {
      *reason = SNES_CONVERGED_PNORM_RELATIVE;
    }
  }

  if(*reason>0)
  {
    gss_log.string_buf()<<"----------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[*reason]<<"    *residual norm = "<<fnorm<<"\n\n\n";
    gss_log.record();
  }


  return(0);
}


PetscErrorCode TrustRegion_ConvergenceTest_L1E(SNES snes,PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm,
    SNESConvergedReason *reason,void *dummy)
{
  SNES_TR *neP = (SNES_TR *)snes->data;
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)dummy;
  if (!ps->its)
  {
    gss_log.string_buf()<<" "<<"its\t"<<"| Eq(V) |   "<<"| Eq(n) |   "<<"| Eq(p) |   "<<"|delta x|\n"
    <<"----------------------------------------------------------------------\n";
    gss_log.record();
  }

  gss_log.string_buf().precision(3);
  gss_log.string_buf()<<"  "<<ps->its++<<"\t"
  <<ps->possion_norm<<"   "
  <<ps->elec_continuty_norm<<"   "
  <<ps->hole_continuty_norm<<"   "
  <<pnorm<<"\n";
  gss_log.record();
  gss_log.string_buf().precision(6);

  *reason = SNES_CONVERGED_ITERATING;

  if (!it)
  {
    snes->ttol = fnorm*snes->rtol;
  }

  if (fnorm != fnorm)
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  else if (neP->delta < xnorm * snes->deltatol)
  {
    if( ps->possion_norm        < ps->toler_relax*ps->possion_abs_toler &&
        ps->elec_continuty_norm < ps->toler_relax*ps->elec_continuty_abs_toler &&
        ps->hole_continuty_norm < ps->toler_relax*ps->hole_continuty_abs_toler &&
        ps->electrode_norm      < ps->toler_relax*ps->electrode_abs_toler)
    {
      *reason = SNES_CONVERGED_TR_DELTA;
    }
    else
    {
      *reason = SNES_DIVERGED_LOCAL_MIN;
    }
  }
  else if (ps->possion_norm < ps->possion_abs_toler &&
           ps->elec_continuty_norm < ps->elec_continuty_abs_toler &&
           ps->hole_continuty_norm < ps->hole_continuty_abs_toler &&
           ps->electrode_norm      < ps->electrode_abs_toler)
  {
    *reason = SNES_CONVERGED_FNORM_ABS;
  }
  else if (snes->nfuncs >= snes->max_funcs)
  {
    *reason = SNES_DIVERGED_FUNCTION_COUNT;
  }

  if (it && !*reason)
  {
    if (fnorm <= snes->ttol)
    {
      *reason = SNES_CONVERGED_FNORM_RELATIVE;
    }
    else if (pnorm < ps->relative_toler &&
             ps->possion_norm        < ps->toler_relax*ps->possion_abs_toler &&
             ps->elec_continuty_norm < ps->toler_relax*ps->elec_continuty_abs_toler &&
             ps->hole_continuty_norm < ps->toler_relax*ps->hole_continuty_abs_toler &&
             ps->electrode_norm      < ps->toler_relax*ps->electrode_abs_toler)
    {
      *reason = SNES_CONVERGED_PNORM_RELATIVE;
    }
  }

  if (snes->nfuncs >= snes->max_funcs)
  {
    *reason = SNES_DIVERGED_FUNCTION_COUNT;
  }
  if(*reason>0)
  {
    gss_log.string_buf()<<"----------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[*reason]<<"    *residual norm = "<<fnorm<<"\n\n\n";
    gss_log.record();
  }

  return(0);
}



PetscErrorCode LimitorByBankRose_L1E(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  int it;
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)checkctx;
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;
  VecGetArray(x,&xx);
  VecGetArray(y,&yy);
  VecGetArray(w,&ww);
  
  SNESGetIterationNumber(snes,&it);
  static PetscScalar K;
  if(it==0) K=0;
  
  Vec gk0,gk1;
  int flag=0;
  VecDuplicate(x,&gk0);
  VecDuplicate(x,&gk1);
  SNESComputeFunction(snes,x,gk0);
  PetscScalar gk0_norm;
  VecNorm(gk0,NORM_2,&gk0_norm); 
  while(1)
  {
    PetscScalar tk=1.0/(1+K*gk0_norm);
    if(tk<1e-4) break;
    //w=x-tk*y;
    flag=1;
    VecWAXPY(w,-tk,y,x);
    SNESComputeFunction(snes,w,gk1);
    PetscScalar gk1_norm;
    VecNorm(gk1,NORM_2,&gk1_norm); 
    if((1-gk1_norm/gk0_norm)<0.9*tk)
    {
      if(K==0) K=1;
      else K*=10;
    }
    else
    {
      K/=10;
      break;
    }
  }
  
  //prevent negative carrier density
  int offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
	if(fabs(yy[offset])>1.0)  {ww[offset] = xx[offset] - dsign(yy[offset])*1.0;flag=1;}
	if (ww[offset+1] <0 ) {ww[offset+1]=1.0*pow(ps->scale_unit.s_centimeter,-3);flag=1;}
	if (ww[offset+2] <0 ) {ww[offset+2]=1.0*pow(ps->scale_unit.s_centimeter,-3);flag=1;}
        offset += 3;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
    }
  }
    
  if(flag)
  {
    *changed_y = PETSC_FALSE;
    *changed_w = PETSC_TRUE;
  }  
  VecDestroy(gk0);
  VecDestroy(gk1);
  VecRestoreArray(x,&xx);
  VecRestoreArray(y,&yy);
  VecRestoreArray(w,&ww);
  return(0);
}


PetscErrorCode LimitorByPotential_L1E(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;
  VecGetArray(x,&xx);
  VecGetArray(y,&yy);
  VecGetArray(w,&ww);
  PetscScalar dV_max=0;
  
  //search for dV_max;
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)checkctx;
  int    offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        if(fabs(yy[offset])>dV_max) {dV_max=fabs(yy[offset]); }
        offset += 3;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
    }
  }
  if(dV_max<1e-2) return(0);
  
  //do logarithmic potential damping;
  PetscScalar Vt=0.026*ps->LatticeTemp;
  PetscScalar f = log(1+dV_max/Vt)/(dV_max/Vt);
  
  //prevent negative carrier density
  offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        ww[offset] = xx[offset] - f*yy[offset];
	if (ww[offset+1] <0 ) ww[offset+1]=1.0*pow(ps->scale_unit.s_centimeter,-3);
	if (ww[offset+2] <0 ) ww[offset+2]=1.0*pow(ps->scale_unit.s_centimeter,-3);
        offset += 3;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        ww[offset] = xx[offset] - f*yy[offset];
        offset += 1;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        ww[offset] = xx[offset] - f*yy[offset];
        offset += 1;
      }
    }
  }
 
  VecRestoreArray(x,&xx);
  VecRestoreArray(y,&yy);
  VecRestoreArray(w,&ww);

  *changed_y = PETSC_FALSE;
  *changed_w = PETSC_TRUE;
  return(0);
}


PetscErrorCode LimitorNonNegativeCarrier_L1E(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;
  VecGetArray(x,&xx);
  VecGetArray(y,&yy);
  VecGetArray(w,&ww);
  
  DDM_Solver_L1E *ps = (DDM_Solver_L1E*)checkctx;
  int flag=0;
  int offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      { 
        if(fabs(yy[offset])>1.0)  {ww[offset] = xx[offset] - dsign(yy[offset])*1.0;flag=1;}
        if (ww[offset+1] <0 ) {ww[offset+1]=1.0*pow(ps->scale_unit.s_centimeter,-3);flag=1;}
	if (ww[offset+2] <0 ) {ww[offset+2]=1.0*pow(ps->scale_unit.s_centimeter,-3);flag=1;}
        offset += 3;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
    }
  }
  VecRestoreArray(x,&xx);
  VecRestoreArray(y,&yy);
  VecRestoreArray(w,&ww);
  
  if(flag) 
  {
    *changed_y = PETSC_FALSE;
    *changed_w = PETSC_TRUE;
  }
  return(0);
}

void ProjectionNonNegativeCarrier_L1E(Vec x, Vec xo, DDM_Solver_L1E *ps)
{
  PetscScalar    *xx;
  PetscScalar    *oo;
  VecGetArray(x,&xx);
  VecGetArray(xo,&oo);
  
  int offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      { 
        if (xx[offset+1] <0 )
          xx[offset+1] = 0.01 * oo[offset+1];
        if (xx[offset+2] <0 )
          xx[offset+2] = 0.01 * oo[offset+2];
        offset += 3;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      offset += pzonedata->pzone->davcell.size();
    }
  }
  VecRestoreArray(x,&xx);
  VecRestoreArray(xo,&oo);
}



/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::init_solver:  This function do initial setup for nonlinear solver
 */
int DDM_Solver_L1E::init_solver(SolveDefine &sv)
{
  gss_log.string_buf()<<"DDM solver Level 1E init...";

  //set Tolerances
  relative_toler           = sv.relative_toler;
  toler_relax              = sv.toler_relax;
  possion_abs_toler        = sv.possion_abs_toler;
  elec_continuty_abs_toler = sv.elec_continuty_abs_toler;
  hole_continuty_abs_toler = sv.hole_continuty_abs_toler;
  electrode_abs_toler      = sv.electrode_abs_toler;

  //if user defined a matrix-free method
  PetscOptionsHasName(PETSC_NULL,"-snes_mf_operator",&mf_flg);

  // compute the scale of problem
  zofs.resize(zone_num+1);
  zofs[0] = 0;

  for(int i=0;i<zone_num;i++)
  {
    if(zonedata[i]->material_type==Semiconductor) //semiconductor zone
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      N += 3*pzonedata->pzone->davcell.size();
      N += pzonedata->electrode.size();  //additional equation for electrode
      zofs[i+1] = N;
    }
    else if(zonedata[i]->material_type==Insulator) //Insulator zone
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      N += zone[i].davcell.size();
      N += pzonedata->electrode.size();  //additional equation for electrode
      zofs[i+1] = N;
    }
    else if(zonedata[i]->material_type==Conductor) //Electrode zone
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      N += zone[i].davcell.size();
      zofs[i+1] = N;
    }
    else
    {
      zofs[i+1] = N;
    }
  }

  VecCreateSeq(PETSC_COMM_SELF,N,&x);
  VecDuplicate(x,&r);
  VecDuplicate(x,&x_n);
  VecDuplicate(x,&x_n1);
  VecDuplicate(x,&x_n2);
  VecDuplicate(x,&xp);
  VecDuplicate(x,&LTE);
  //Evaluate initial guess,importnat for newton solver
  PetscScalar init_value[3];
  int    index[3];
  int    offset=0;
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        index[0] = offset+0;
        index[1] = offset+1;
        index[2] = offset+2;
        init_value[0] = pzonedata->fs[i].P;
        init_value[1] = pzonedata->fs[i].n;
        init_value[2] = pzonedata->fs[i].p;
        VecSetValues(x,3,index,init_value,INSERT_VALUES);
        offset += 3;
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        VecSetValue(x,offset,bc.Get_pointer(pzonedata->electrode[j])->Get_Potential(),INSERT_VALUES);
        offset += 1;
      }
    }

    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        VecSetValue(x,offset,pzonedata->fs[i].P,INSERT_VALUES);
        offset += 1;
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        VecSetValue(x,offset,bc.Get_pointer(pzonedata->electrode[j])->Get_Potential(),INSERT_VALUES);
        offset += 1;
      }
    }

    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        VecSetValue(x,offset,pzonedata->fs[i].P,INSERT_VALUES);
        offset += 1;
      }
    }

  }
  //
  VecAssemblyBegin(x);
  VecAssemblyEnd(x);

  // Create Jacobian matrix data structure
  // pre-alloc approximate memory
  int *nnz = new int[N];
  int *p = nnz;
  for(int i=0;i<zone_num;i++)
  {
    //semiconductor zone
    if(zonedata[i]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      for(int j=0;j<pzonedata->pzone->davcell.size();j++)
      {
        const VoronoiCell* pcell = pzonedata->pzone->davcell.GetPointer(j);
        //inner node, alloc exact memory.
        if(!pcell->bc_index || bc[pcell->bc_index-1].psegment->interface==-1)
        {
          *p++ = 3*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 3*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 3*(pzonedata->pzone->davcell[j].nb_num+1);
        }
        else //interface node, slightly more than needed.
        {
          *p++ = 3*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 3*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 3*(2*pzonedata->pzone->davcell[j].nb_num+4);
        }
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
        *p++ = bc[pzonedata->electrode[j]].psegment->node_array.size()+1;

    }
    //Insulator zones, poisson equation only
    if(zonedata[i]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      for(int j=0;j<zone[i].davcell.size();j++)
        *p++ = 1*(zone[i].davcell[j].nb_num+1);
      for(int j=0;j<pzonedata->electrode.size();j++)
        *p++ = bc[pzonedata->electrode[j]].psegment->node_array.size()+1;
    }
    //Electrode zones, poisson equation only
    if(zonedata[i]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      for(int j=0;j<zone[i].davcell.size();j++)
        *p++ = 1*(zone[i].davcell[j].nb_num+1);
    }
  }

  MatCreate(PETSC_COMM_SELF,&J);
  MatSetSizes(J,N,N,N,N);
  if(sv.LS=="superlu")
  {
#ifdef PETSC_HAVE_SUPERLU  
    MatSetType(J,MATSUPERLU);
#else
    MatSetType(J,MATSEQAIJ);    
#endif    
  }
  else if(sv.LS=="umfpack")
  {
#ifdef PETSC_HAVE_UMFPACK
    MatSetType(J,MATUMFPACK);
#else
    MatSetType(J,MATSEQAIJ);    
#endif    
  }
  else
  {
    MatSetType(J,MATSEQAIJ);
  }
  MatSeqAIJSetPreallocation(J,0,nnz);
  
  MatCreateSeqAIJ(PETSC_COMM_SELF,N,N,0,nnz,&JTmp);
  if(mf_flg)
    MatCreateSeqAIJ(PETSC_COMM_SELF,N,N,0,nnz,&JPrec);
  
  delete [] nnz;


  SNESCreate(PETSC_COMM_WORLD,&snes);
  SNESSetFunction(snes,r,SNES_form_function_pn_1E,this);
  if(mf_flg)
  {
    gss_log.string_buf()<<"Matrix-Free...";
    SNESSetJacobian(snes,J,JPrec,SNESMF_form_jacobian_pn_1E,this);
  }
  else
  {
    SNESSetJacobian(snes,J,J,SNES_form_jacobian_pn_1E,this);
  }

  // the maximum number of iterations
  PetscInt maxit = sv.maxit;
  if(sv.Type==EQUILIBRIUM)  maxit = 10*sv.maxit;

  SNESSetType(snes,SNESLS);  //default method
  if(sv.NS==LineSearch || sv.NS==Basic)
  {
    if(sv.NS==LineSearch)
      SNESLineSearchSet(snes,SNESLineSearchCubic,PETSC_NULL);
    if(sv.NS==Basic)
      SNESLineSearchSet(snes,SNESLineSearchNo,PETSC_NULL);
    
    if(sv.Damping==DampingBankRose)
       SNESLineSearchSetPostCheck(snes,LimitorByBankRose_L1E,this);
    else if(sv.Damping==DampingPotential)
       SNESLineSearchSetPostCheck(snes,LimitorByPotential_L1E,this);
    else
       SNESLineSearchSetPostCheck(snes,LimitorNonNegativeCarrier_L1E,this);
    
    SNESSetTolerances(snes,1e-12*N,1e-14,1e-9,maxit,1000);
    SNESSetConvergenceTest(snes,LineSearch_ConvergenceTest_L1E,this);
  }
  else if(sv.NS==TrustRegion)
  {
    SNESSetType(snes,SNESTR);
    SNESSetTolerances(snes,1e-12*N,1e-14,1e-9,maxit,1000);
    SNESSetTrustRegionTolerance(snes,1e-30);
    SNESSetConvergenceTest(snes,TrustRegion_ConvergenceTest_L1E,this);
    //must set TR delta0 to a sufficient large value, or TR can't find real solution.
    SNES_TR  *neP = (SNES_TR*)snes->data;
    neP->delta0 = N;
  }

  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp,&pc);
  //KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL,PETSC_NULL);
  if(sv.LS=="lu"||sv.LS=="superlu"||sv.LS=="umfpack")
  {
    KSPSetType(ksp,KSPPREONLY);
    PCSetType(pc,PCLU);
    PCFactorSetReuseOrdering(pc,PETSC_TRUE);
    PCFactorSetPivoting(pc,1.0);
    PCFactorReorderForNonzeroDiagonal(pc,1e-14);
    PCFactorSetShiftNonzero(pc,1e-12);
  }
  else
  {
    KSPSetType(ksp,sv.LS.c_str());
    if(sv.LS=="gmres")  KSPGMRESSetRestart(ksp,150);
    KSPSetTolerances(ksp,1e-10*N,1e-20*N,PETSC_DEFAULT,N/10);
    PCSetType(pc,PCILU);
    PCFactorSetLevels(pc,3);
    PCFactorSetShiftNonzero(pc,1e-14);
    PCFactorReorderForNonzeroDiagonal(pc,1e-12);
  }
  KSPSetFromOptions(ksp);
  SNESSetFromOptions(snes);
  gss_log.string_buf()<<"done\n";
  gss_log.record();
  return 0;
}



/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::do_solve:  This function solve the problem
 */
int DDM_Solver_L1E::do_solve(SolveDefine &sv)
{
  if(sv.Type==EQUILIBRIUM)
    solve_equ(sv);

  else if(sv.Type==STEADYSTATE)
    solve_steadystate(sv);

  else if(sv.Type==DCSWEEP)
    solve_dcsweep(sv);

  else if(sv.Type==TRANSIENT)
    solve_transient(sv);

  else if(sv.Type==TRACE)
    solve_iv_trace(sv);

  else
  {
    gss_log.string_buf()<<"DDML1E not support this solver type!\n";
    gss_log.record();
  }
  return 0;

}

/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solve_equ:  This function compute equilibrium
 * all the stimulate source(s) are zero. time step set to inf
 */
int DDM_Solver_L1E::solve_equ(SolveDefine &sv)
{
  gss_log.string_buf()<<"Compute equilibrium\n";
  gss_log.record();

  bc.Clear_Vapp();
  bc.Clear_Iapp();
  ODE_formula.TimeDependent = false;
  ODE_formula.dt = 1e100;
  ODE_formula.clock = 0.0;
  its = 0;

  for(int i=0;i<zone_num;i++)
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->HighFieldMobility= false;
      pzonedata->ImpactIonization = false;
      pzonedata->IIType           = sv.IIType;
      pzonedata->BandBandTunneling = false;
      pzonedata->IncompleteIonization = false;
      pzonedata->QuantumMechanical = false;
      pzonedata->Fermi             = sv.Fermi;
      pzonedata->EJModel           = sv.EJModel;
    }
  // for equilibrium compute, low field mobility is ok.
  solver_pre_compute();
  probe_open(EQUILIBRIUM);
  SNESSolve(snes,PETSC_NULL,x);
  solution_update();
  probe(EQUILIBRIUM,0);
  probe_close();
  SNESGetConvergedReason(snes,&reason);
  if(reason<0)
  {
    gss_log.string_buf()<<"----------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason]<<"    !residual norm = "<<norm<<"\n\n\n";
    gss_log.record();
  }
  return 0;
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solve_steadystate:  This function compute steadystate
 * set all the stimulate source(s) at there value of t=0. time step set to inf
 */
int DDM_Solver_L1E::solve_steadystate(SolveDefine &sv)
{

  gss_log.string_buf()<<"Compute steady-state\n";
  gss_log.record();

  bc.Update_Vapp(0);
  bc.Update_Iapp(0);
  optgen_update(0);
  ODE_formula.TimeDependent = false;
  ODE_formula.dt = 1e100;
  ODE_formula.clock = 0.0;
  its = 0;

  for(int i=0;i<zone_num;i++)
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->HighFieldMobility= sv.HighFieldMobility;
      pzonedata->ImpactIonization = sv.ImpactIonization;
      pzonedata->IIType           = sv.IIType;
      pzonedata->BandBandTunneling = sv.BandBandTunneling;
      pzonedata->IncompleteIonization = sv.IncompleteIonization;
      pzonedata->QuantumMechanical = sv.QuantumMechanical;
      pzonedata->Fermi             = sv.Fermi;
      pzonedata->EJModel           = sv.EJModel;
    }
  solver_pre_compute();
  probe_open(STEADYSTATE);
  SNESSolve(snes,PETSC_NULL,x);
  solution_update();
  probe(STEADYSTATE,0);
  probe_close();
  SNESGetConvergedReason(snes,&reason);
  if(reason<0)
  {
    gss_log.string_buf()<<"----------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[reason]<<"    !residual norm = "<<norm<<"\n\n\n";
    gss_log.record();
  }
  return 0;

}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solve_dcsweep:  This function compute IV curve
 * time step set to inf
 */
int DDM_Solver_L1E::solve_dcsweep(SolveDefine &sv)
{
  PetscScalar I=0;
  FILE *fiv;
  PetscScalar current_scale_mA = scale_unit.s_mA;
  PetscScalar voltage_scale_V =scale_unit.s_volt;

  bc.Update_Vapp(0);
  bc.Update_Iapp(0);
  optgen_update(0);
  ODE_formula.TimeDependent = false;
  ODE_formula.dt = 1e100;
  ODE_formula.clock = 0.0;

  //set which electrodes are required to record IV
  for(int i=0;i<zone_num;i++)
  {
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->HighFieldMobility= sv.HighFieldMobility;
      pzonedata->ImpactIonization = sv.ImpactIonization;
      pzonedata->IIType           = sv.IIType;
      pzonedata->BandBandTunneling = sv.BandBandTunneling;
      pzonedata->IncompleteIonization = sv.IncompleteIonization;
      pzonedata->QuantumMechanical = sv.QuantumMechanical;
      pzonedata->Fermi             = sv.Fermi;
      pzonedata->EJModel           = sv.EJModel;
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
          {  
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
            sv.Electrode_Record_Index.push_back(k);
          }
    }
    if(zonedata[i]->material_type == Insulator)
    {
      ISZone * pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
          {  
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
            sv.Electrode_Record_Index.push_back(k);
          }
    }
  }


  // output DC Scan information
  if(sv.Electrode_VScan!=-1)
  {
    if(!sv.Electrode_VScan_Name.size())
    {
      gss_log.string_buf()<<"No VScan Electrode Specified.\n";
      gss_log.record();
      return 1;
    } 
    gss_log.string_buf()<<"DC voltage scan from "<<sv.VStart
    <<" step "<<sv.VStep
    <<" to "<<sv.VStop<<"\n";
    gss_log.record();
  }
  else
  {
    gss_log.string_buf()<<"DC current scan from "<<sv.IStart/current_scale_mA
    <<" step "<<sv.IStep/current_scale_mA
    <<" to "<<sv.IStop/current_scale_mA<<"\n";
    gss_log.record();
  }

  // prepare IV file. If no file given, output to screen.
  if(!sv.IVFile.empty())
  {
    if(sv.IVFileAppend)
    {
      fiv=fopen(sv.IVFile.c_str(),"a");
      if(!fiv) fiv=stdout;
    }
    else
    {
      fiv=fopen(sv.IVFile.c_str(),"w");
      if(!fiv) fiv=stdout;
      fprintf(fiv,"#");
      for(int j=0;j<sv.Electrode_Record.size();j++)
        fprintf(fiv,"Vp(%s)[V]  I(%s)[mA]   ",
                sv.Electrode_Record_Name[sv.Electrode_Record_Index[j]].c_str(),
                sv.Electrode_Record_Name[sv.Electrode_Record_Index[j]].c_str()
               );
      fprintf(fiv,"\n");
    }
  }
  else
    fiv = stdout;
    
  // begin
  if(sv.Electrode_VScan!=-1)
  {
    PetscScalar V = sv.VStart;
    PetscScalar V_n, V_n1, V_n2;
    for(int i=0;i<sv.Electrode_VScan_Name.size();i++)
      bc.Set_electrode_type(sv.Electrode_VScan_Name[i].c_str(),VoltageBC);
    int diverged_retry=0; 
    stack<PetscScalar> V_retry;
    probe_open(DCSWEEP_VSCAN);
    for (int step=0; V*sv.VStep < sv.VStop *sv.VStep *(1+1e-7);)
    {
      its = 0;
      gss_log.string_buf()<<"DC Scan: V("<<sv.Electrode_VScan_Name[0];
      for(int i=1;i<sv.Electrode_VScan_Name.size();i++) 
        gss_log.string_buf()<<", "<<sv.Electrode_VScan_Name[i];
      gss_log.string_buf()<<") = "<<V/voltage_scale_V<<" V"<<"\n";
      gss_log.record();
      for(int i=0;i<sv.Electrode_VScan_Name.size();i++)
        bc.Set_Vapp(sv.Electrode_VScan_Name[i].c_str(),V);
      solver_pre_compute();
      SNESSolve(snes,PETSC_NULL,x);
      SNESGetConvergedReason(snes,&reason);

      if(reason>0 && norm==norm) //ok, converged.
      {
        diverged_retry=0;
        solution_update();
        probe(DCSWEEP_VSCAN,V);
        // save solution
        VecCopy(x_n1,x_n2); V_n2=V_n1;
        VecCopy(x_n,x_n1); V_n1=V_n;
        VecCopy(x,x_n); V_n=V;
        step++;
        
        if (V_retry.empty())
          V+=sv.VStep;
        else
        {
          V=V_retry.top();
          V_retry.pop();
        }

        if(fabs(V-sv.VStop)<1e-10)
          V=sv.VStop;

        for(int j=0;j<sv.Electrode_Record.size();j++)
        {
          fprintf(fiv,"%lf\t%e\t",
              double(bc.Get_pointer(sv.Electrode_Record[j])->Get_Potential_new()/voltage_scale_V),
              double(bc.Get_pointer(sv.Electrode_Record[j])->Get_Current_new()/current_scale_mA));
        }
        if(sv.Electrode_Record.size())
        {
          fprintf(fiv,"\n");
          fflush(fiv);
        }
      }
      else // oh, diverged... reduce step and try again
      {
        if(step==0)
        {
          gss_log.string_buf()<<"------> Failed in the first step.\n\n\n";
          gss_log.record();
          break;
        }
        if(V_retry.size()>=8)
        {
          gss_log.string_buf()<<"------> Too many failed steps, give up trying.\n\n\n";
          gss_log.record();
          break;
        }
        // recover last solution
        VecCopy(x_n,x);
        V_retry.push(V);
        V=(V+V_n)/2.0;
        gss_log.string_buf()<<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        gss_log.record();
      }

      if (sv.Projection)
      {
        // solution projection
        PetscScalar hn = V-V_n;
        PetscScalar hn1 = V_n-V_n1;
        PetscScalar hn2 = V_n1-V_n2;
        if (step>=3)
        {
          // quadradic projection
          PetscScalar cn=hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
          PetscScalar cn1=-hn*(hn+hn1+hn2)/(hn1*hn2);
          PetscScalar cn2=hn*(hn+hn1)/(hn2*(hn1+hn2));

          VecAXPY(x,cn,x_n);
          VecAXPY(x,cn1,x_n1);
          VecAXPY(x,cn2,x_n2);
          ProjectionNonNegativeCarrier_L1E(x,x_n,this);
        }
        else if (step>=2)
        {
          // linear projection
          VecAXPY(x, hn/hn1,x_n);
          VecAXPY(x,-hn/hn1,x_n1);
          ProjectionNonNegativeCarrier_L1E(x,x_n,this);
        }
      }
    }
  }

  if(sv.Electrode_IScan!=-1)
  {
    bc.Set_electrode_type(sv.Electrode_IScan_Name.c_str(),CurrentBC);
    PetscScalar I = sv.IStart;
    PetscScalar I_n, I_n1, I_n2;
    int diverged_retry=0;
    stack<PetscScalar> I_retry;
    probe_open(DCSWEEP_ISCAN);
    for (int step=0; I*sv.IStep < sv.IStop *sv.IStep *(1+1e-7); )
    {
      its = 0;
      gss_log.string_buf()<<"DC Scan: I("<<sv.Electrode_IScan_Name<<") = "<<I/current_scale_mA<<" mA"<<"\n";
      gss_log.record();
      bc.Set_Iapp(sv.Electrode_IScan_Name.c_str(),I);
      solver_pre_compute();
      SNESSolve(snes,PETSC_NULL,x);
      SNESGetConvergedReason(snes,&reason);

      if(reason>0 && norm==norm) //ok, converged.
      {
        diverged_retry=0;
        solution_update();
        probe(DCSWEEP_ISCAN,I);
        // save solution
        VecCopy(x_n1,x_n2); I_n2=I_n1;
        VecCopy(x_n,x_n1); I_n1=I_n;
        VecCopy(x,x_n); I_n=I;
        step++;
        
        if (I_retry.empty())
          I+=sv.IStep;
        else
        {
          I=I_retry.top();
          I_retry.pop();
        }

        if(fabs(I-sv.IStop)<1e-10)
          I=sv.IStop;
        for(int j=0;j<sv.Electrode_Record.size();j++)
        {
          fprintf(fiv,"%lf\t%e\t",
              double(bc.Get_pointer(sv.Electrode_Record[j])->Get_Potential_new()/voltage_scale_V),
              double(bc.Get_pointer(sv.Electrode_Record[j])->Get_Current_new()/current_scale_mA));
        }
        if(sv.Electrode_Record.size())
        {
          fprintf(fiv,"\n");
          fflush(fiv);
        }
      }
      else // oh, diverged... reduce step and try again
      {
        if(step==0)
        {
          gss_log.string_buf()<<"------> Failed in the first step.\n\n\n";
          gss_log.record();
          break;
        }
        if(I_retry.size()>=8) //failed 8 times, stop tring
        {
          gss_log.string_buf()<<"------> Too many failed steps, give up tring.\n\n\n";
          gss_log.record();
          break;
        }
        VecCopy(x_n,x);
        I_retry.push(I);
        I=(I+I_n)/2.0;
        gss_log.string_buf()<<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        gss_log.record();
      }

      if (sv.Projection)
      {
        // solution projection
        PetscScalar hn = I-I_n;
        PetscScalar hn1 = I_n-I_n1;
        PetscScalar hn2 = I_n1-I_n2;
        if (step>=3)
        {
          // quadradic projection
          PetscScalar cn=hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
          PetscScalar cn1=-hn*(hn+hn1+hn2)/(hn1*hn2);
          PetscScalar cn2=hn*(hn+hn1)/(hn2*(hn1+hn2));

          VecAXPY(x,cn,x_n);
          VecAXPY(x,cn1,x_n1);
          VecAXPY(x,cn2,x_n2);
          ProjectionNonNegativeCarrier_L1E(x,x_n,this);
        }
        else if (step>=2)
        {
          // linear projection
          VecAXPY(x, hn/hn1,x_n);
          VecAXPY(x,-hn/hn1,x_n1);
          ProjectionNonNegativeCarrier_L1E(x,x_n,this);
        }
      }
    }
  }

  if(!sv.IVFile.empty())        fclose(fiv);
  probe_close();
  
  return 0;

}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solve_iv_trace:  This function use continuation method to trace 
 * IV curve automatically
 */
int DDM_Solver_L1E::solve_iv_trace(SolveDefine &sv)
{
  int         error=0;
  int         first_step=1;
  int         slope_flag=0;
  Vec         pdI_pdw;
  Vec         pdF_pdV;
  Vec         pdw_pdV;
  PetscScalar pdI_pdV;
  PetscScalar dI_dV;
  PetscScalar I=0;
  PetscScalar V = sv.VStart;
  PetscScalar VStep = sv.VStep;
  PetscScalar r,Rload=0,Rload_new,Rref;
  PetscScalar angle,slope=PI/2,slope_new,slope_chord;
  FILE *fiv;
  PetscScalar current_scale_mA = scale_unit.s_mA;
  PetscScalar voltage_scale_V =scale_unit.s_volt;
  int      bc_index;
  BaseBC * pbc = NULL;
  SMCZone *pz = NULL;
  int      electrode_index;
  KSP      kspc;
  PC       pcc;
  VecCreateSeq(PETSC_COMM_SELF,N,&pdI_pdw);
  VecDuplicate(pdI_pdw,&pdF_pdV);
  VecDuplicate(pdI_pdw,&pdw_pdV);
  KSPCreate(MPI_COMM_SELF,&kspc);
  KSPSetType(kspc,KSPPREONLY);
  KSPGetPC(kspc,&pcc);
  PCSetType(pcc,PCLU);
  PCFactorReorderForNonzeroDiagonal(pcc,1e-14);
  PCFactorSetShiftNonzero(pcc,1e-12);

  bc.Update_Vapp(0);
  bc.Update_Iapp(0);
  bc.Set_electrode_type(sv.Electrode_VScan_Name[0].c_str(),VoltageBC);
  bc.Set_Vapp(sv.Electrode_VScan_Name[0].c_str(),V);
  ODE_formula.TimeDependent = false;
  ODE_formula.dt = 1e100;

  // output TRACE information
  gss_log.string_buf()<<"IV automatically trace by continuation method\n";
  gss_log.record();

  //set which electrodes are required to do trace or record IV
  for(int i=0;i<zone_num;i++)
  {
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->HighFieldMobility= sv.HighFieldMobility;
      pzonedata->ImpactIonization = sv.ImpactIonization;
      pzonedata->IIType           = sv.IIType;
      pzonedata->BandBandTunneling = sv.BandBandTunneling;
      pzonedata->IncompleteIonization = sv.IncompleteIonization;
      pzonedata->QuantumMechanical = sv.QuantumMechanical;
      pzonedata->Fermi             = sv.Fermi;
      pzonedata->EJModel           = sv.EJModel;
      
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
          {  
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
            sv.Electrode_Record_Index.push_back(k);
          }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        if(bc.is_electrode_label(pzonedata->electrode[j],sv.Electrode_VScan_Name[0].c_str()))
        {
          pz = pzonedata;
          electrode_index = j;
          bc_index = pzonedata->electrode[j];
          pbc = bc.Get_pointer(bc_index);
        }
      }
    }
    if(zonedata[i]->material_type == Insulator)
    {
      ISZone * pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
          {  
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
            sv.Electrode_Record_Index.push_back(k);
          }
    }
  }

  //check error
  if(!pbc || !pz )
  {
    gss_log.string_buf()<<"I can't find suitable electrode to do IV trace.\n";
    gss_log.record();
    error = 1;
    goto trace_end;
  }
  if(pbc->BCType!=OhmicContact && pbc->BCType!=SchottkyContact)
  {
    gss_log.string_buf()<<"IV trace can only be performanced to Ohmic or Schottky contact.\n";
    gss_log.record();
    error = 1;
    goto trace_end;
  }

  // prepare IV file. If no file given, output to screen.
  if(!sv.IVFile.empty())
  {
    if(sv.IVFileAppend)
    {
      fiv=fopen(sv.IVFile.c_str(),"a");
      if(!fiv) fiv=stdout;
    }
    else
    {
      fiv=fopen(sv.IVFile.c_str(),"w");
      if(!fiv) fiv=stdout;
      fprintf(fiv,"#");
      for(int j=0;j<sv.Electrode_Record.size();j++)
        fprintf(fiv,"Vp(%s)[V]  I(%s)[mA]   ",
                sv.Electrode_Record_Name[sv.Electrode_Record_Index[j]].c_str(),
                sv.Electrode_Record_Name[sv.Electrode_Record_Index[j]].c_str()
               );
      fprintf(fiv,"\n");
    }
  }
  else
    fiv = stdout;

  // initial condition check
  pbc->Set_R(Rload);
  solver_pre_compute();
  SNESSolve(snes,PETSC_NULL,x);
  SNESGetConvergedReason(snes,&reason);
  if(reason<0)
  {
    gss_log.string_buf()<<"I can't get convergence even at initial point, please give a good initial condition.\n\n";
    gss_log.record();
    error = 1;
    goto trace_end;
  }
  solution_update();
  // output IV curve
  for(int j=0;j<sv.Electrode_Record.size();j++)
  {
    I = bc.Get_pointer(sv.Electrode_Record[j])->Get_Current_new();
    fprintf(fiv,"%lf\t%e\t",
            double(bc.Get_pointer(sv.Electrode_Record[j])->Get_Potential()/voltage_scale_V),
            double(I/current_scale_mA));
  }
  if(sv.Electrode_Record.size())
  {
    fprintf(fiv,"\n");
    fflush(fiv);
  }

  //loop here
  while(pbc->Get_Potential()*sv.VStep < sv.VStop*sv.VStep && fabs(pbc->Get_Current())<sv.IStop)
  {
    gss_log.string_buf()<<"Trace for V("<<sv.Electrode_VScan_Name[0]<<")="<<pbc->Get_Potential()<<"(V)\n";
    gss_log.record();

    // for last Rload, increase VStep
    int recovery=0;
    do
    {
      its = 0;
      if(!slope_flag && first_step)
        V+=VStep;
      if(!slope_flag && !first_step)
        V*=1+VStep*cos(slope)/pbc->Get_Potential();
      bc.Set_Vapp(sv.Electrode_VScan_Name[0].c_str(),V);
      SNESSolve(snes,PETSC_NULL,x);
      SNESGetConvergedReason(snes,&reason);
      if(reason<0)
      {
        gss_log.string_buf()<<"I can't get convergence at this step, do recovery...\n\n";
        gss_log.record();
        diverged_recovery();
        V/=1+VStep*cos(slope)/pbc->Get_Potential();
        VStep/=2;
        recovery++;
      }
      if(recovery>8)
      {
        gss_log.string_buf()<<"Too many failed steps, give up tring.\n\n\n";
        gss_log.record();
        error = 1;
        goto trace_end;
      }
    }
    while(reason<0);

    // for the new bias point, recompute slope
    // calculate the dynamic resistance of IV curve by different approximation
    //r = (pbc->Get_Potential_new()-pbc->Get_Potential())/(pbc->Get_Current_new()-pbc->Get_Current());

    if(pbc->BCType==OhmicContact)
      pz->F1E_om_electrode_Trace(electrode_index,pdI_pdV,pdI_pdw,pdF_pdV,x,&JTmp,&J,ODE_formula,zofs,bc,DeviceDepth);
    if(pbc->BCType==SchottkyContact)
      pz->F1E_stk_electrode_Trace(electrode_index,pdI_pdV,pdI_pdw,pdF_pdV,x,&JTmp,&J,ODE_formula,zofs,bc,DeviceDepth);
    KSPSetOperators(kspc,J,J,SAME_NONZERO_PATTERN);
    KSPSetUp(kspc);
    KSPSolve(kspc,pdF_pdV,pdw_pdV);
    VecDot(pdI_pdw,pdw_pdV,&dI_dV);
    r=1/dI_dV;
    //printf("dI/dV=%e pdI/pdV=%e VStep=%e\n",1/r,dI_dV,VStep);

    if(first_step)
    {
      Rref = pbc->Get_Potential_new()/pbc->Get_Current_new();
      slope=slope_new=atan(Rref/r);
    }
    else
    {
      Rref = pbc->Get_Potential()/pbc->Get_Current();     // Rref, for scaling the dynamic resistance
      slope_new = atan(Rref/r);
    }

    //resolve turning point
    if(tan(slope)*tan(slope_new)<0) //the sign is changed
    {
      PetscScalar chord = (pbc->Get_Potential_new()-pbc->Get_Potential())/(pbc->Get_Current_new()-pbc->Get_Current());
      slope_chord = atan(Rref/chord);
      //printf("slope_chord=%e slope=%e slope_new=%e \n",slope_chord/PI*180, slope/PI*180, slope_new/PI*180);
      //for vertical turning point
      if(fabs(slope)/PI*180 >70 && fabs(slope_new)/PI*180 >70 &&
          ((slope_chord*slope>0 && fabs(slope_chord)>fabs(slope)) || (slope_chord*slope_new>0 && fabs(slope_chord)>fabs(slope_new))))
      {
        VStep*=-1.0;
        //gss_log.string_buf()<<"vertical turning point...\n";
        //gss_log.record();
      }
      //for horizontal turning point
      if( fabs(slope)/PI*180 <20 && fabs(slope_new)/PI*180 <20 &&
          ((slope_chord*slope>0 && fabs(slope_chord)<fabs(slope)) || (slope_chord*slope_new>0 && fabs(slope_chord)<fabs(slope_new))))
      {
        VStep*=1.0;
        //gss_log.string_buf()<<"horizontal turning point...\n";
        //gss_log.record();
      }
    }

    // check the slope change
    angle=fabs(slope-slope_new);
    if(angle>PI/2) angle=PI-angle;
    if(angle<PI/36)      VStep*=1.5;   // slope change less than 5 degree
    else if(angle<PI/24) VStep=VStep;  // slope change less than 10 degree, but greater than 5  degree
    else if(angle<PI/12) VStep/=2;     // slope change less than 15 degree, but greater than 10 degree
    else                                               // slope greater than 15 degree, reject this solution
    {
      //printf("slope change old=%e new=%e angle=%e \n",slope/PI*180, slope_new/PI*180,  angle/PI*180);
      gss_log.string_buf()<<"Slope of IV curve changes too quickly, do recovery...\n\n";
      gss_log.record();
      //V/=1+VStep*cos(slope)/pbc->Get_Potential();
      VStep/=2;
      //diverged_recovery();
      //continue;
    }
    if(fabs(VStep*cos(slope))>sv.VStepMax) VStep=VStep/fabs(VStep)*sv.VStepMax/cos(slope);
    // ok, update solutions
    solution_update();

    // since Rload will be changed, we should change bias voltage to meet the truncation point.

    Rload_new=Rref/r;
    V+=pbc->Get_Current()*(Rload_new-Rload);
    pbc->Set_R(Rload_new);
    bc.Set_Vapp(sv.Electrode_VScan_Name[0].c_str(),V);
    Rload=Rload_new;
    slope = slope_new;
    // output IV curve
    for(int j=0;j<sv.Electrode_Record.size();j++)
    {
      I = bc.Get_pointer(sv.Electrode_Record[j])->Get_Current_new();
      fprintf(fiv,"%lf\t%e\t",
              double(bc.Get_pointer(sv.Electrode_Record[j])->Get_Potential_new()/voltage_scale_V),
              double(I/current_scale_mA));
    }
    if(sv.Electrode_Record.size())
    {
      fprintf(fiv,"\n");
      fflush(fiv);
    }

    first_step=0;
  }

trace_end:
  KSPDestroy(kspc);
  VecDestroy(pdI_pdw);
  VecDestroy(pdF_pdV);
  VecDestroy(pdw_pdV);
  if(!sv.IVFile.empty())        fclose(fiv);
  return error;
}



void DDM_Solver_L1E::LET_norm_estimat(PetscScalar & r)
{
  PetscScalar hn= ODE_formula.dt;
  PetscScalar hn1= ODE_formula.dt_last;
  PetscScalar hn2= ODE_formula.dt_last_last;
  PetscScalar eps_r=1e-3,eps_a=1e-4;
  
  VecZeroEntries(xp);
  VecZeroEntries(LTE);
    
  if(ODE_formula.BDF_Type==BDF1)
  {
    VecAXPY(xp,1+hn/hn1,x_n);
    VecAXPY(xp,-hn/hn1,x_n1);
    VecAXPY(LTE,hn/(hn+hn1),x);
    VecAXPY(LTE,-hn/(hn+hn1),xp);
  }
  else if(ODE_formula.BDF_Type==BDF2)
  {
    PetscScalar cn=1+hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
    PetscScalar cn1=-hn*(hn+hn1+hn2)/(hn1*hn2);
    PetscScalar cn2=hn*(hn+hn1)/(hn2*(hn1+hn2));
    
    VecAXPY(xp,cn,x_n);
    VecAXPY(xp,cn1,x_n1);
    VecAXPY(xp,cn2,x_n2);
    VecAXPY(LTE, hn/(hn+hn1+hn2),x);
    VecAXPY(LTE,-hn/(hn+hn1+hn2),xp);
  }    
    
  PetscScalar    *xx,*ll;
  VecGetArray(x,&xx);
  VecGetArray(LTE,&ll);
  int offset=0;
  int N=0;
  for(int z=0;z<zone_num;z++)
  {
      if(zonedata[z]->material_type==Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
        for(int i=0;i<pzonedata->pzone->davcell.size();i++)
        {
          ll[offset+0]=0;
	  ll[offset+1]=ll[offset+1]/(eps_r*xx[offset+1]+eps_a);
	  ll[offset+2]=ll[offset+2]/(eps_r*xx[offset+2]+eps_a);
          offset += 3;
	  N+=2;
        }
	for(int i=0;i<pzonedata->electrode.size();i++)
        {
	  ll[offset+0]=0;
	  offset += 1;
	}  
      }
      else if(zonedata[z]->material_type==Insulator)
      {
        ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
        for(int i=0;i<pzonedata->pzone->davcell.size();i++)
        {
          ll[offset]=0;
          offset += 1;
        }
        for(int i=0;i<pzonedata->electrode.size();i++)
        {
	  ll[offset+0]=0;
	  offset += 1;
	}
      }
      else if(zonedata[z]->material_type==Conductor)
      {
        ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
        for(int i=0;i<pzonedata->pzone->davcell.size();i++)
        {
          ll[offset]=0;
          offset += 1;
        }
      }
  }
  VecRestoreArray(x,&xx);
  VecRestoreArray(LTE,&ll);
  VecNorm(LTE,NORM_2,&r);
      
  r/=sqrt((PetscScalar)N);
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solve_transient:  This function does transient simulation.
 */
int DDM_Solver_L1E::solve_transient(SolveDefine &sv)
{

  PetscScalar I=0;
  FILE *fiv;

  PetscScalar current_scale_mA = scale_unit.s_mA;
  PetscScalar voltage_scale_V =scale_unit.s_volt;
  int diverged_retry=0;
  int step=0;
  PetscScalar r=1.0;
  clock = sv.TStart+sv.TStep;
  ODE_formula.TimeDependent = true;
  ODE_formula.BDF_Type = sv.BDF_Type;
  ODE_formula.dt = sv.TStep;
  if(sv.BDF_Type==BDF2)
  {
    ODE_formula.BDF2_restart = true;
  }

  for(int i=0;i<zone_num;i++)
  {
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->HighFieldMobility= sv.HighFieldMobility;
      pzonedata->ImpactIonization = sv.ImpactIonization;
      pzonedata->IIType           = sv.IIType;
      pzonedata->BandBandTunneling = sv.BandBandTunneling;
      pzonedata->IncompleteIonization = sv.IncompleteIonization;
      pzonedata->QuantumMechanical = sv.QuantumMechanical;
      pzonedata->Fermi             = sv.Fermi;
      pzonedata->EJModel           = sv.EJModel;
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
          { 
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
            sv.Electrode_Record_Index.push_back(k);
          }
    }
    if(zonedata[i]->material_type == Insulator)
    {
      ISZone * pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
          {
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
            sv.Electrode_Record_Index.push_back(k);
          }
    }
  }

  gss_log.string_buf()<<"Transient compute from "<<sv.TStart<<" ps step "<<sv.TStep<<" ps to "<<sv.TStop<<" ps"<<"\n";
  gss_log.record();

  //open record file
  if(!sv.IVFile.empty())
  {
    if(sv.IVFileAppend)
    {
      fiv=fopen(sv.IVFile.c_str(),"a");
      if(!fiv) fiv=stdout;
    }
    else
    {
      fiv=fopen(sv.IVFile.c_str(),"w");
      if(!fiv) fiv=stdout;
      fprintf(fiv,"#");
      fprintf(fiv,"time[ps]  ");
      for(int j=0;j<sv.Electrode_Record.size();j++)
        fprintf(fiv,"Vp(%s)[V]  I(%s)[mA]   ",
                sv.Electrode_Record_Name[sv.Electrode_Record_Index[j]].c_str(),
                sv.Electrode_Record_Name[sv.Electrode_Record_Index[j]].c_str()
               );
      fprintf(fiv,"\n");
    }
  }
  else
    fiv = stdout;
  
  probe_open(TRANSIENT);

  // the main loop of transient solver.
  do 
  {
    gss_log.string_buf()<<"t = "<<clock<<" ps"<<"\n";
    gss_log.record();
    ODE_formula.clock = clock;
    
    its = 0; //clear nonlinear iteration counter
    //update sources
    bc.Update_Vapp(clock);
    bc.Update_Iapp(clock);
    optgen_update(clock);
    
    //we do solve here!
    solver_pre_compute();
    SNESSolve(snes,PETSC_NULL,x);
    SNESGetConvergedReason(snes,&reason);

    //nonlinear solution diverged?
    if(reason<0 || !(norm==norm))
    {
      if(++diverged_retry>=8) //failed 8 times, stop tring
      {
        gss_log.string_buf()<<"------> Too many failed steps, give up tring.\n\n\n";
        gss_log.record();
        break;
      }
      diverged_recovery();
      gss_log.string_buf()<<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
      gss_log.record();
      ODE_formula.dt/=2.0;  // reduce time step by a factor of two
      clock-=ODE_formula.dt;
      if(clock<sv.TStart) clock=sv.TStart;
      continue;
    }
    
    //ok, nonlinear solution converged.
    diverged_retry=0;
    
    //do LTE estimation and auto time step control
    if(sv.AutoStep && ((ODE_formula.BDF_Type==BDF1 && step>=3) || (ODE_formula.BDF_Type==BDF2 && step>=4)))
    {
      LET_norm_estimat(r);

      if(ODE_formula.BDF_Type==BDF1)  
        r = pow(r,PetscScalar(-1.0/2));
      else if(ODE_formula.BDF_Type==BDF2)
        r = pow(r,PetscScalar(-1.0/3));
      
      if(r<0.9) //reject this solution
      {
        diverged_recovery();
        gss_log.string_buf()<<"------> LTE too large, time step rejected...\n\n\n";
        gss_log.record();
	clock-=ODE_formula.dt;
	ODE_formula.dt*=r;  // reduce time step by a factor of r
        clock+=ODE_formula.dt;
        continue;
      }
      ODE_formula.dt=ODE_formula.dt * (r > 2.0 ? 2.0 : r);
      if(ODE_formula.dt>4*sv.TStep) ODE_formula.dt=4*sv.TStep;
    }
    else
    { 
       if(fabs(ODE_formula.dt) < fabs(sv.TStep))
         ODE_formula.dt*=1.1;
    }
    
    //update solution
    step++;
    solution_update();
    VecCopy(x_n1,x_n2);
    VecCopy(x_n,x_n1);
    VecCopy(x,x_n);
    ODE_formula.dt_last_last = ODE_formula.dt_last;
    ODE_formula.dt_last = ODE_formula.dt;

    //predict next solution
    if(sv.Predict && (ODE_formula.BDF_Type==BDF1 && step>=3) || (ODE_formula.BDF_Type==BDF2 && step>=4))
    {
      PetscScalar hn= ODE_formula.dt;
      PetscScalar hn1= ODE_formula.dt_last;
      PetscScalar hn2= ODE_formula.dt_last_last;
      VecZeroEntries(x);
      if(ODE_formula.BDF_Type==BDF1)
      {
        VecAXPY(x,1+hn/hn1,x_n);
        VecAXPY(x,-hn/hn1,x_n1);
      }
      else if(ODE_formula.BDF_Type==BDF2)
      {
        PetscScalar cn=1+hn*(hn+2*hn1+hn2)/(hn1*(hn1+hn2));
        PetscScalar cn1=-hn*(hn+hn1+hn2)/(hn1*hn2);
        PetscScalar cn2=hn*(hn+hn1)/(hn2*(hn1+hn2));
    
        VecAXPY(x,cn,x_n);
        VecAXPY(x,cn1,x_n1);
        VecAXPY(x,cn2,x_n2);
      }    
    }
    
    //save solution to file
    probe(TRANSIENT,clock);
    if(sv.Electrode_Record.size())
    {
      fprintf(fiv,"%e\t",double(clock));
      for(int j=0;j<sv.Electrode_Record.size();j++)
      {
        I = bc.Get_pointer(sv.Electrode_Record[j])->Get_Current_new();
        fprintf(fiv,"%lf\t%e\t",
                double(bc.Get_pointer(sv.Electrode_Record[j])->Get_Potential_new()/voltage_scale_V),
                double(I/current_scale_mA));
      }
      fprintf(fiv,"%lf\t",double(r));
      if(sv.Electrode_Record.size())
      {
        fprintf(fiv,"\n");
        fflush(fiv);
      }
    }
    
    //make sure we can terminat near TStop, the relative error should less than 1e-10.  
    clock+=ODE_formula.dt;
    if(clock > sv.TStop && clock < (sv.TStop + ODE_formula.dt - 1e-10*ODE_formula.dt))
    {
        ODE_formula.dt -= clock - sv.TStop; 
        clock=sv.TStop;
    }
    
    //clear the first step flag of BDF2     
    if(ODE_formula.BDF_Type==BDF2)
      ODE_formula.BDF2_restart = false;

      
  }while(clock<sv.TStop+0.5*ODE_formula.dt);

  //close record file
  if(!sv.IVFile.empty())
    fclose(fiv);
  probe_close();
  
  return 0;

}



/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solver_pre_compute:  This function precompute some parameters before
 * each newton iteration.
 */
void DDM_Solver_L1E::solver_pre_compute()
{
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        SemiData *pfs = &(pzonedata->fs[i]);
        SemiAuxData *paux = &(pzonedata->aux[i]);
        pzonedata->mt->mapping(&zone[z].danode[i],paux,0);
        pzonedata->aux[i].mun =  pzonedata->mt->mob->ElecMob(pfs->p,pfs->n,pfs->T,0,0,pfs->T);
        pzonedata->aux[i].mup =  pzonedata->mt->mob->HoleMob(pfs->p,pfs->n,pfs->T,0,0,pfs->T);
      }
    }
  }
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solver_post_compute:  This function do post process after each 
 * DC sweep or transient simulation
 */
void DDM_Solver_L1E::solver_post_compute()
{}



/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::solution_update:  This function restore solution data from SNES
 */
void DDM_Solver_L1E::solution_update()
{
  //------------------------------------------------------------
  PetscScalar    *xx;
  VecGetArray(x,&xx);
  int    offset=0;
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      PetscScalar kb = pzonedata->mt->kb;
      PetscScalar e = pzonedata->mt->e;
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        pzonedata->fs[i].P_last = pzonedata->fs[i].P;
        pzonedata->fs[i].n_last = pzonedata->fs[i].n;
        pzonedata->fs[i].p_last = pzonedata->fs[i].p;

        pzonedata->fs[i].P = xx[offset+0];
        pzonedata->fs[i].n = fabs(xx[offset+1]);
        pzonedata->fs[i].p = fabs(xx[offset+2]);

        pzonedata->mt->mapping(&pzonedata->pzone->danode[i],&pzonedata->aux[i],0);
        PetscScalar nie = pzonedata->mt->band->nie(pzonedata->fs[i].T);
        pzonedata->aux[i].Ec = -(e*pzonedata->fs[i].P + pzonedata->aux[i].affinity + pzonedata->mt->band->EgNarrowToEc(pzonedata->fs[i].T));//conduction band energy level
        pzonedata->aux[i].Ev = -(e*pzonedata->fs[i].P + pzonedata->aux[i].affinity - pzonedata->mt->band->EgNarrowToEv(pzonedata->fs[i].T) + pzonedata->mt->band->Eg(pzonedata->fs[i].T));//valence band energy level
        pzonedata->aux[i].phi_intrinsic = -0.5*( pzonedata->aux[i].Ec+pzonedata->aux[i].Ev + kb*pzonedata->fs[i].T*log(pzonedata->aux[i].Nv/pzonedata->aux[i].Nc))/e;
        if (pzonedata->Fermi)
        {
          pzonedata->aux[i].phin = -(pzonedata->aux[i].Ec + kb*pzonedata->fs[i].T*fermi_mhalf(fabs(pzonedata->fs[i].n)/pzonedata->aux[i].Nc))/e;
          pzonedata->aux[i].phip = -(pzonedata->aux[i].Ev - kb*pzonedata->fs[i].T*fermi_mhalf(fabs(pzonedata->fs[i].p)/pzonedata->aux[i].Nv))/e;
        }
        else
        {
          pzonedata->aux[i].phin = -(pzonedata->aux[i].Ec + kb*pzonedata->fs[i].T*log(fabs(pzonedata->fs[i].n)/pzonedata->aux[i].Nc))/e;
          pzonedata->aux[i].phip = -(pzonedata->aux[i].Ev - kb*pzonedata->fs[i].T*log(fabs(pzonedata->fs[i].p)/pzonedata->aux[i].Nv))/e;
        }

        pzonedata->fs[i].Eqc = -e*(pzonedata->fs[i].P+pzonedata->aux[i].affinity);
        pzonedata->fs[i].Eqv = -e*(pzonedata->fs[i].P+pzonedata->aux[i].affinity+pzonedata->aux[i].Eg);
        offset += 3;
      }
      pzonedata->F1E_efield_update(xx,zofs,bc,zonedata);
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar P = xx[offset];
        PetscScalar Pn = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        PetscScalar C =  bc.Get_pointer(pzonedata->electrode[i])->Get_C();
        PetscScalar I =  bc.Get_pointer(pzonedata->electrode[i])->Get_Current_new();
        PetscScalar In=  bc.Get_pointer(pzonedata->electrode[i])->Get_Current();
        bc.Get_pointer(pzonedata->electrode[i])->Set_Cap_Current(C*(P-Pn)/ODE_formula.dt);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current(I);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current_old(In);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Potential(P);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Potential_old(Pn);
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        pzonedata->fs[i].P_last = pzonedata->fs[i].P;
        pzonedata->fs[i].P = xx[offset];
        offset += 1;
      }
      pzonedata->F1_efield_update(xx,zofs,bc,zonedata);
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar P = xx[offset];
        PetscScalar Pn = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        PetscScalar C =  bc.Get_pointer(pzonedata->electrode[i])->Get_C();
        PetscScalar I =  bc.Get_pointer(pzonedata->electrode[i])->Get_Current_new();
        PetscScalar In=  bc.Get_pointer(pzonedata->electrode[i])->Get_Current();
        bc.Get_pointer(pzonedata->electrode[i])->Set_Cap_Current(C*(P-Pn)/ODE_formula.dt);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current(I);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current_old(In);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Potential(P);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Potential_old(Pn);
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        pzonedata->fs[i].P_last = pzonedata->fs[i].P;
        pzonedata->fs[i].P = xx[offset];
        offset += 1;
      }
    }
  }
  VecRestoreArray(x,&xx);
}

/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::diverged_recovery:  This function recovery latest solution data
 * if SNES diverged.
 */
void DDM_Solver_L1E::diverged_recovery()
{
  //------------------------------------------------------------
  PetscScalar    *xx;
  VecGetArray(x,&xx);
  int    offset=0;
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        xx[offset+0] = pzonedata->fs[i].P;
        xx[offset+1] = pzonedata->fs[i].n;
        xx[offset+2] = pzonedata->fs[i].p;
        offset += 3;
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar P = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        xx[offset] = P;
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        xx[offset] = pzonedata->fs[i].P;
        offset += 1;
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar P = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        xx[offset] = P;
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        xx[offset] = pzonedata->fs[i].P;
        offset += 1;
      }
    }
  }
  VecRestoreArray(x,&xx);
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1E::destroy_solver:  This function do destroy the nonlinear solver
 */
int DDM_Solver_L1E::destroy_solver(SolveDefine &sv)
{
  // free work space
  N = 0;
  zofs.clear();
  VecDestroy(x);
  VecDestroy(r);
  VecDestroy(x_n);
  VecDestroy(x_n1);
  VecDestroy(x_n2);
  VecDestroy(xp);
  VecDestroy(LTE);
  MatDestroy(J);
  MatDestroy(JTmp);
  if(mf_flg)
    MatDestroy(JPrec);
  SNESDestroy(snes);
  return 0;
}
