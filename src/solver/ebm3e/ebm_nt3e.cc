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
/*  Last update: May 9, 2007                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "mathfunc.h"
#include "ebm_nt3e.h"
#include "private/kspimpl.h"
#include "private/snesimpl.h"
#include "src/snes/impls/tr/tr.h"
#include "log.h"



/* ----------------------------------------------------------------------------
 * form_function_pn:  This function setup DDM equation F(x)=0
 */
void EBM_Solver_L3E::form_function_pn_3E(PetscScalar *x,PetscScalar *f)
{
  // compute flux along triangle edges. semiconductor zone only
  for(int z=0;z<zone_num;z++)
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].datri.size();i++)
      {
        Tri *ptri = &zone[z].datri[i];
        pzonedata->F3E_Tri_ddm(ptri,x,f,zofs);
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
          pzonedata->F3E_om_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==SchottkyContact)
          pzonedata->F3E_stk_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==InsulatorContact)
          pzonedata->F3E_ins_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
      }

      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index)
          pzonedata->F3E_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F3E_ddm_neumannbc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
        {
          OhmicBC *pbc = dynamic_cast<OhmicBC*>(bc.Get_pointer(pcell->bc_index-1));
          if(pbc->pinterface)//the ohmic bc is an electrode region
          {
            int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
            int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
            ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
            pzonedata->F3E_ddm_ombc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
          }
          else //the ohmic bc is a segment
            pzonedata->F3E_ddm_ombc_segment(i,x,f,ODE_formula,zofs,bc);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
        {
          SchottkyBC *pbc = dynamic_cast<SchottkyBC*>(bc.Get_pointer(pcell->bc_index-1));
          if(pbc->pinterface)//the Schottky bc is an electrode region
          {
            int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
            int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
            ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
            pzonedata->F3E_ddm_stkbc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
          }
          else //the Schottky bc is a segment
            pzonedata->F3E_ddm_stkbc_segment(i,x,f,ODE_formula,zofs,bc);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorContact)
          pzonedata->F3E_ddm_insulator_gate(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F3E_ddm_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HomoInterface)
        {
          HomoInterfaceBC *pbc;
          pbc = dynamic_cast<HomoInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F3E_ddm_homojunction(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->F3E_ddm_inner(i,x,f,ODE_formula,zofs);
      }

    }

    // Insulator zone
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata= dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index)
          pzonedata->F3_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F3_ddm_neumannbc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
        {
          ChargedContactBC *pbc = dynamic_cast<ChargedContactBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_chargebc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==GateContact)
        {
          GateBC *pbc = dynamic_cast<GateBC*>(bc.Get_pointer(pcell->bc_index-1));
          if(pbc->pinterface)//the gate bc is an electrode region
          {
            int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
            int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
            ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
            pzonedata->F3_ddm_gatebc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
          }
          else //the gate bc is a segment
            pzonedata->F3_ddm_gatebc_segment(i,x,f,ODE_formula,zofs,bc);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
        {
          ElectrodeInsulatorBC *pbc = dynamic_cast<ElectrodeInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_electrode_insulator_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_semiconductor_insulator_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Insulator_Insulator)
        {
          InsulatorInsulatorBC *pbc;
          pbc = dynamic_cast<InsulatorInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_insulator_insulator_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->F3_ddm_inner(i,x,f,ODE_formula,zofs);
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==GateContact)
          pzonedata->F3_gate_electrode(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          pzonedata->F3_charge_electrode(i,x,f,ODE_formula,zofs,bc);
      }
    }

    // Electrode zone
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata= dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index)
          pzonedata->F3_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F3_ddm_neumannbc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
        {
          ChargedContactBC *pbc = dynamic_cast<ChargedContactBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_chargebc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==GateContact)
        {
          GateBC *pbc = dynamic_cast<GateBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_gatebc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
        {
          ElectrodeInsulatorBC *pbc = dynamic_cast<ElectrodeInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_elec_insulator_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
        {
          OhmicBC *pbc = dynamic_cast<OhmicBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_ombc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
        {
          SchottkyBC *pbc = dynamic_cast<SchottkyBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->F3_ddm_stkbc_interface(i,x,f,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->F3_ddm_inner(i,x,f,ODE_formula,zofs);
      }
    }

  }

}


/* ----------------------------------------------------------------------------
 * form_jacobian_pn:  This function setup jacobian matrix F'(x) of DDM equation F(x)=0
 * dual carrier edition
 */
void EBM_Solver_L3E::form_jacobian_pn_3E(PetscScalar *x, Mat *jac, Mat *jtmp)
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
        pzonedata->J3E_Tri_ddm(ptri,x,jtmp,zofs);
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
          pzonedata->J3E_om_electrode(i,x,jac,jtmp,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==SchottkyContact)
          pzonedata->J3E_stk_electrode(i,x,jac,jtmp,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==InsulatorContact)
          pzonedata->J3E_ins_electrode(i,x,jac,jtmp,ODE_formula,zofs,bc,DeviceDepth);
      }

      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        SemiData *pfs = &pzonedata->fs[i];
        if(!pcell->bc_index)
          pzonedata->J3E_ddm_inner(i,x,jac,jtmp,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->J3E_ddm_neumannbc(i,x,jac,jtmp,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
        {
          OhmicBC *pbc = dynamic_cast<OhmicBC*>(bc.Get_pointer(pcell->bc_index-1));
          if(pbc->pinterface)//the ohmic bc is an electrode region
          {
            int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
            int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
            ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
            pzonedata->J3E_ddm_ombc_interface(i,x,jac,jtmp,ODE_formula,zofs,bc,pz,n_node);
          }
          else //the ohmic bc is a segment
            pzonedata->J3E_ddm_ombc_segment(i,x,jac,jtmp,ODE_formula,zofs,bc);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
        {
          SchottkyBC *pbc = dynamic_cast<SchottkyBC*>(bc.Get_pointer(pcell->bc_index-1));
          if(pbc->pinterface)//the Schottky bc is an electrode region
          {
            int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
            int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
            ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
            pzonedata->J3E_ddm_stkbc_interface(i,x,jac,jtmp,ODE_formula,zofs,bc,pz,n_node);
          }
          else //the Schottky bc is a segment
            pzonedata->J3E_ddm_stkbc_segment(i,x,jac,jtmp,ODE_formula,zofs,bc);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorContact)
          pzonedata->J3E_ddm_insulator_gate(i,x,jac,jtmp,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J3E_ddm_interface(i,x,jac,jtmp,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HomoInterface)
        {
          HomoInterfaceBC *pbc;
          pbc = dynamic_cast<HomoInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J3E_ddm_homojunction(i,x,jac,jtmp,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->J3E_ddm_inner(i,x,jac,jtmp,ODE_formula,zofs);
      }
    }

    // Insulator zone
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata= dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index)
          pzonedata->J3_ddm_inner(i,x,jac,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->J3_ddm_neumannbc(i,x,jac,ODE_formula,zofs,bc);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
        {
          ChargedContactBC *pbc = dynamic_cast<ChargedContactBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_chargebc_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==GateContact)
        {
          GateBC *pbc = dynamic_cast<GateBC*>(bc.Get_pointer(pcell->bc_index-1));
          if(pbc->pinterface)//the gate bc is an electrode region
          {
            int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
            int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
            ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
            pzonedata->J3_ddm_gatebc_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
          }
          else //the gate bc is a segment
            pzonedata->J3_ddm_gatebc_segment(i,x,jac,ODE_formula,zofs,bc);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
        {
          ElectrodeInsulatorBC *pbc = dynamic_cast<ElectrodeInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_electrode_insulator_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_semiconductor_insulator_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Insulator_Insulator)
        {
          InsulatorInsulatorBC *pbc;
          pbc = dynamic_cast<InsulatorInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_insulator_insulator_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->J3_ddm_inner(i,x,jac,ODE_formula,zofs);
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==GateContact)
          pzonedata->J3_gate_electrode(i,x,jac,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          pzonedata->J3_charge_electrode(i,x,jac,ODE_formula,zofs,bc);
      }
    }

    // Electrode zone
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata= dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index)
          pzonedata->J3_ddm_inner(i,x,jac,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->J3_ddm_neumannbc(i,x,jac,ODE_formula,zofs,bc);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
        {
          ChargedContactBC *pbc = dynamic_cast<ChargedContactBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_chargebc_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==GateContact)
        {
          GateBC *pbc = dynamic_cast<GateBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_gatebc_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
        {
          ElectrodeInsulatorBC *pbc = dynamic_cast<ElectrodeInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_elec_insulator_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
        {
          OhmicBC *pbc = dynamic_cast<OhmicBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_ombc_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
        {
          SchottkyBC *pbc = dynamic_cast<SchottkyBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->J3_ddm_stkbc_interface(i,x,jac,ODE_formula,zofs,bc,pz,n_node);
        }
        else //prevent some mistakes caused by inexact bc
          pzonedata->J3_ddm_inner(i,x,jac,ODE_formula,zofs);
      }
    }

  }

}



/* ----------------------------------------------------------------------------
 * SNES_form_function_pn_3E:  wrap function for petsc nonlinear solver
 */
PetscErrorCode SNES_form_function_pn_3E(SNES snes, Vec x,Vec f,void *dummy)
{
  PetscScalar    *xx,*ff;
  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)dummy;
  VecZeroEntries(f);
  //Get pointers to vector data.
  VecGetArray(x,&xx);
  VecGetArray(f,&ff);

  ps->form_function_pn_3E(xx,ff);
  ps->error_norm_pn_3E(xx,ff);

  //Restore vectors
  VecRestoreArray(x,&xx);
  VecRestoreArray(f,&ff);
  VecNorm(f,NORM_2,&ps->norm);

  return 0;
}


/* ----------------------------------------------------------------------------
 * SNES_form_jacobian_pn_2E:  wrap function for petsc nonlinear solver
 */
PetscErrorCode SNES_form_jacobian_pn_3E(SNES snes, Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
  PetscScalar    *xx;
  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)dummy;
  //Get pointer to vector data
  VecGetArray(x,&xx);
  //clear old matrix
  MatZeroEntries(*jac);
  MatZeroEntries(ps->JTmp);
  //build matrix here
  ps->form_jacobian_pn_3E(xx,jac,&ps->JTmp);
  *flag = SAME_NONZERO_PATTERN;
  //Restore vector
  VecRestoreArray(x,&xx);
  //Assemble matrix
  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  //MatView(*jac,PETSC_VIEWER_DRAW_WORLD);
  //getchar();
  return 0;
}


/* ----------------------------------------------------------------------------
 * SNESMF_form_jacobian_pn_3E:  wrap function for petsc nonlinear solver with
 * matrix-free method
 */
PetscErrorCode SNESMF_form_jacobian_pn_3E(SNES snes, Vec x,Mat *A,Mat *B,MatStructure *flag,void *dummy)
{
  PetscScalar    *xx;
  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)dummy;
  //Get pointer to vector data
  VecGetArray(x,&xx);
  //clear old matrix
  MatZeroEntries(*B);
  MatZeroEntries(ps->JTmp);
  //build matrix here
  ps->form_jacobian_pn_3E(xx,B,&ps->JTmp);
  *flag = SAME_NONZERO_PATTERN;
  //Restore vector
  VecRestoreArray(x,&xx);
  //Assemble matrix
  MatAssemblyBegin(*B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*B,MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);
  return 0;
}


/* ----------------------------------------------------------------------------
 * error_norm_pn_3E:  This function compute X and RHS error norms.
 */
void EBM_Solver_L3E:: error_norm_pn_3E(PetscScalar *x,PetscScalar *f)
{
  int offset=0;

  // do clear
  potential_norm   = 0;
  electron_norm    = 0;
  hole_norm        = 0;
  temperature_norm = 0;
  elec_temperature_norm = 0;
  hole_temperature_norm = 0;
  possion_norm        = 0;
  elec_continuty_norm = 0;
  hole_continuty_norm = 0;
  heat_equation_norm  = 0;
  elec_energy_equation_norm = 0;
  hole_energy_equation_norm = 0;
  electrode_norm      = 0;

  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        potential_norm   += x[offset+0]*x[offset+0];
        electron_norm    += x[offset+1]*x[offset+1];
        hole_norm        += x[offset+2]*x[offset+2];
        temperature_norm += x[offset+3]*x[offset+3];
        elec_temperature_norm += x[offset+4]/x[offset+1]*x[offset+4]/x[offset+1];
        hole_temperature_norm += x[offset+5]/x[offset+2]*x[offset+5]/x[offset+2];

        possion_norm        += f[offset+0]*f[offset+0];
        elec_continuty_norm += f[offset+1]*f[offset+1];
        hole_continuty_norm += f[offset+2]*f[offset+2];
        heat_equation_norm  += f[offset+3]*f[offset+3];
        elec_energy_equation_norm  += f[offset+4]*f[offset+4];
        hole_energy_equation_norm  += f[offset+5]*f[offset+5];
        offset += 6;
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
        potential_norm   += x[offset]*x[offset];
        temperature_norm += x[offset+1]*x[offset+1];
        possion_norm        += f[offset]*f[offset];
        heat_equation_norm  += f[offset+1]*f[offset+1];
        offset += 2;
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
        potential_norm   += x[offset]*x[offset];
        temperature_norm += x[offset+1]*x[offset+1];
        possion_norm        += f[offset]*f[offset];
        heat_equation_norm  += f[offset+1]*f[offset+1];
        offset += 2;
      }
    }
  }

  potential_norm = sqrt(potential_norm);
  electron_norm  = sqrt(electron_norm);
  hole_norm      = sqrt(hole_norm);
  temperature_norm = sqrt(temperature_norm);
  elec_temperature_norm = sqrt(elec_temperature_norm);
  hole_temperature_norm = sqrt(hole_temperature_norm);
  possion_norm        = sqrt(possion_norm);
  elec_continuty_norm = sqrt(elec_continuty_norm);
  hole_continuty_norm = sqrt(hole_continuty_norm);
  electrode_norm      = sqrt(electrode_norm);
  heat_equation_norm  = sqrt(heat_equation_norm);
  elec_energy_equation_norm = sqrt(elec_energy_equation_norm);
  hole_energy_equation_norm = sqrt(hole_energy_equation_norm);
}


PetscErrorCode LineSearch_ConvergenceTest_L3E(SNES snes,PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm,
    SNESConvergedReason *reason,void *dummy)
{
  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)dummy;

  *reason = SNES_CONVERGED_ITERATING;
  if (!it)
  {
    snes->ttol = fnorm*snes->rtol;

    gss_log.string_buf()<<" "<<"its\t"<<"| Eq(V) | "<<"| Eq(n) | "<<"| Eq(p) | "<<"| Eq(T) | "<<"|Eq(Tn)|  "<<"|Eq(Tp)|  "<<"|delta x|\n"
    <<"-----------------------------------------------------------------------------\n";
    gss_log.record();
  }

  gss_log.string_buf().precision(2);
  gss_log.string_buf()<<"  "<<it<<"\t"
  <<ps->possion_norm<<"  "
  <<ps->elec_continuty_norm<<"  "
  <<ps->hole_continuty_norm<<"  "
  <<ps->heat_equation_norm<<"  "
  <<ps->elec_energy_equation_norm<<"  "
  <<ps->hole_energy_equation_norm<<"  "
  <<pnorm<<"\n";
  gss_log.record();
  gss_log.string_buf().precision(6);

  if (fnorm != fnorm)
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  else if (ps->possion_norm < ps->possion_abs_toler &&
           ps->elec_continuty_norm < ps->elec_continuty_abs_toler &&
           ps->hole_continuty_norm < ps->hole_continuty_abs_toler &&
           ps->electrode_norm      < ps->electrode_abs_toler &&
           ps->heat_equation_norm  < ps->heat_equation_abs_toler &&
           ps->elec_energy_equation_norm <ps->elec_energy_abs_toler &&
           ps->hole_energy_equation_norm <ps->hole_energy_abs_toler)
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
             ps->electrode_norm      < ps->toler_relax*ps->electrode_abs_toler &&
             ps->heat_equation_norm  < ps->toler_relax*ps->heat_equation_abs_toler &&
             ps->elec_energy_equation_norm <ps->toler_relax*ps->elec_energy_abs_toler &&
             ps->hole_energy_equation_norm <ps->toler_relax*ps->hole_energy_abs_toler)
    {
      *reason = SNES_CONVERGED_PNORM_RELATIVE;
    }
  }

  //cout<<pnorm<<" "<<(pnorm != pnorm)<<endl;
  if(*reason>0)
  {
    gss_log.string_buf()<<"----------------------------------------------------------------------\n"
    <<"      "<<SNESConvergedReasons[*reason]<<"    *residual norm = "<<fnorm<<"\n\n\n";
    gss_log.record();
  }

  return(0);
}


PetscErrorCode TrustRegion_ConvergenceTest_L3E(SNES snes,PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm,
    SNESConvergedReason *reason,void *dummy)
{
  SNES_TR *neP = (SNES_TR *)snes->data;
  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)dummy;

  if (!ps->its)
  {
    gss_log.string_buf()<<" "<<"its\t"<<"| Eq(V) | "<<"| Eq(n) | "<<"| Eq(p) | "<<"| Eq(T) | "<<"|Eq(Tn)|  "<<"|Eq(Tp)|  "<<"|delta x|\n"
    <<"-----------------------------------------------------------------------------\n";
    gss_log.record();
  }

  gss_log.string_buf().precision(2);
  gss_log.string_buf()<<"  "<<ps->its++<<"\t"
  <<ps->possion_norm<<"  "
  <<ps->elec_continuty_norm<<"  "
  <<ps->hole_continuty_norm<<"  "
  <<ps->heat_equation_norm<<"  "
  <<ps->elec_energy_equation_norm<<"  "
  <<ps->hole_energy_equation_norm<<"  "
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
  else if (neP->delta < xnorm * snes->deltatol )
  {
    if( ps->possion_norm        < ps->toler_relax*ps->possion_abs_toler &&
        ps->elec_continuty_norm < ps->toler_relax*ps->elec_continuty_abs_toler &&
        ps->hole_continuty_norm < ps->toler_relax*ps->hole_continuty_abs_toler &&
        ps->electrode_norm      < ps->toler_relax*ps->electrode_abs_toler &&
        ps->heat_equation_norm  < ps->toler_relax*ps->heat_equation_abs_toler &&
        ps->elec_energy_equation_norm <ps->toler_relax*ps->elec_energy_abs_toler &&
        ps->hole_energy_equation_norm <ps->toler_relax*ps->hole_energy_abs_toler)
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
           ps->electrode_norm      < ps->electrode_abs_toler &&
           ps->heat_equation_norm  < ps->heat_equation_abs_toler &&
           ps->elec_energy_equation_norm <ps->elec_energy_abs_toler &&
           ps->hole_energy_equation_norm <ps->hole_energy_abs_toler)
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
             ps->electrode_norm      < ps->toler_relax*ps->electrode_abs_toler &&
             ps->heat_equation_norm  < ps->toler_relax*ps->heat_equation_abs_toler &&
             ps->elec_energy_equation_norm <ps->toler_relax*ps->elec_energy_abs_toler &&
             ps->hole_energy_equation_norm <ps->toler_relax*ps->hole_energy_abs_toler)
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


PetscErrorCode LimitorByBankRose_L3E(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  int it;
  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)checkctx;
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

  //prevent negative carrier density and temperature underflow
  int offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        if(fabs(yy[offset])>1.0)                  {ww[offset] = xx[offset] - dsign(yy[offset])*1.0;flag=1;}
        if (ww[offset+3] < 0.7*ps->LatticeTemp )  {ww[offset+3]= 0.7*ps->LatticeTemp;flag=1;}

        pzonedata->mt->mapping(&pzonedata->pzone->danode[i],&pzonedata->aux[i],ps->ODE_formula.clock);
        PetscScalar nie=pzonedata->mt->band->nie(ww[offset+3]);
        PetscScalar Tn=fabs(ww[offset+4]/ww[offset+1]);
        PetscScalar Tp=fabs(ww[offset+5]/ww[offset+2]);

        if (ww[offset+1] < 0 )
        {
          //ww[offset+1]=1*pow(ps->scale_unit.s_centimeter,-3);
          //ww[offset+1]=nie*nie/ww[offset+2];
	  ww[offset+1]=fabs(ww[offset+1]);
          ww[offset+4]=Tn*ww[offset+1];
          flag=1;
        }
        if (ww[offset+2] < 0 )
        {
          //ww[offset+2]=1*pow(ps->scale_unit.s_centimeter,-3);
	  ww[offset+2]=fabs(ww[offset+2]);
          //ww[offset+2]=nie*nie/ww[offset+1];
          ww[offset+5]=Tp*ww[offset+2];
          flag=1;
        }

        if (ww[offset+4] < 0.9*ww[offset+1]*ww[offset+3] )  {ww[offset+4]= 0.9*ww[offset+1]*ww[offset+3];flag=1;}
        if (ww[offset+5] < 0.9*ww[offset+2]*ww[offset+3] )  {ww[offset+5]= 0.9*ww[offset+2]*ww[offset+3];flag=1;}

        offset += 6;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        if (ww[offset+1]<0.7*ps->LatticeTemp ) {ww[offset+1]=0.7*ps->LatticeTemp;flag=1;}
        offset += 2;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        if (ww[offset+1]<0.7*ps->LatticeTemp ) {ww[offset+1]=0.7*ps->LatticeTemp;flag=1;}
        offset += 2;
      }
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


PetscErrorCode LimitorByPotential_L3E(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;
  VecGetArray(x,&xx);
  VecGetArray(y,&yy);
  VecGetArray(w,&ww);
  PetscScalar dV_max=0;

  //search for dV_max;
  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)checkctx;
  int    offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        if(fabs(yy[offset])>dV_max) {dV_max=fabs(yy[offset]); }
        offset += 6;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      offset += 2*pzonedata->pzone->davcell.size();
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      offset += 2*pzonedata->pzone->davcell.size();
    }
  }

  if(dV_max<1e-6) dV_max=1e-6;

  PetscScalar Vt=0.026*ps->LatticeTemp;
  PetscScalar f = log(1+dV_max/Vt)/(dV_max/Vt);

  //prevent negative carrier density and temperature underflow
  offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        
        ww[offset]   = xx[offset] - f*yy[offset];
/*
        ww[offset+1] = xx[offset+1] - f*yy[offset+1];
        ww[offset+2] = xx[offset+2] - f*yy[offset+2];
        ww[offset+3] = xx[offset+3] - f*yy[offset+3];
        ww[offset+4] = xx[offset+4] - f*yy[offset+4];
        ww[offset+5] = xx[offset+5] - f*yy[offset+5];
*/
        if (ww[offset+3] < 0.7*ps->LatticeTemp )  ww[offset+3]= 0.7*ps->LatticeTemp;
        
	pzonedata->mt->mapping(&pzonedata->pzone->danode[i],&pzonedata->aux[i],ps->ODE_formula.clock);
        PetscScalar nie=pzonedata->mt->band->nie(ww[offset+3]);
        PetscScalar Tn=fabs(ww[offset+4]/ww[offset+1]);
        PetscScalar Tp=fabs(ww[offset+5]/ww[offset+2]);
        
	if (ww[offset+1] < 0 )
        {
          ww[offset+1]=fabs(ww[offset+1]);
          //ww[offset+1]=nie*nie/fabs(ww[offset+2]);
	  //ww[offset+1]=1*pow(ps->scale_unit.s_centimeter,-3);
          ww[offset+4]=Tn*ww[offset+1];
        }
        if (ww[offset+2] < 0 )
        {
          ww[offset+2]=fabs(ww[offset+2]);
          //ww[offset+2]=nie*nie/fabs(ww[offset+1]);
	  //ww[offset+2]=1*pow(ps->scale_unit.s_centimeter,-3);
          ww[offset+5]=Tp*ww[offset+2];
        }
        
        if (ww[offset+4] < 0.9*ww[offset+1]*ww[offset+3] )  ww[offset+4]= 0.9*ww[offset+1]*ww[offset+3];
        if (ww[offset+5] < 0.9*ww[offset+2]*ww[offset+3] )  ww[offset+5]= 0.9*ww[offset+2]*ww[offset+3];

        offset += 6;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        ww[offset] = xx[offset] - f*yy[offset];
        if (ww[offset+1]<0.7*ps->LatticeTemp ) ww[offset+1]=0.7*ps->LatticeTemp;
        offset += 2;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        ww[offset] = xx[offset] - f*yy[offset];
        if (ww[offset+1]<0.7*ps->LatticeTemp ) ww[offset+1]=0.7*ps->LatticeTemp;
        offset += 2;
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


PetscErrorCode LimitorNonNegativeCarrier_L3E(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;
  VecGetArray(x,&xx);
  VecGetArray(y,&yy);
  VecGetArray(w,&ww);
  PetscScalar WorstRatio=0.0;

  EBM_Solver_L3E *ps = (EBM_Solver_L3E*)checkctx;
  int flag=0;
  int offset=0;
  for(int z=0;z<ps->zone_num;z++)
  {
    if(ps->zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        if(fabs(yy[offset])>1.0)                  {ww[offset] = xx[offset] - dsign(yy[offset])*1.0;flag=1;}
        if (ww[offset+3] < 0.7*ps->LatticeTemp )  {ww[offset+3]= 0.7*ps->LatticeTemp;flag=1;}

        pzonedata->mt->mapping(&pzonedata->pzone->danode[i],&pzonedata->aux[i],ps->ODE_formula.clock);
        PetscScalar nie=pzonedata->mt->band->nie(ww[offset+3]);
        PetscScalar Tn=fabs(ww[offset+4]/ww[offset+1]);
        PetscScalar Tp=fabs(ww[offset+5]/ww[offset+2]);

        if (ww[offset+1] < 0 )
        {
          //ww[offset+1]=1*pow(ps->scale_unit.s_centimeter,-3);
          //ww[offset+1]=nie*nie/fabs(ww[offset+2]);
	  ww[offset+1]=fabs(ww[offset+1]);
          ww[offset+4]=Tn*ww[offset+1];
          flag=1;
        }
        if (ww[offset+2] < 0 )
        {
          //ww[offset+2]=1*pow(ps->scale_unit.s_centimeter,-3);
          //ww[offset+2]=nie*nie/fabs(ww[offset+1]);
	  ww[offset+2]=fabs(ww[offset+2]);
          ww[offset+5]=Tp*ww[offset+2];
          flag=1;
        }

        if (ww[offset+4] < 0.9*ww[offset+1]*ww[offset+3] )  {ww[offset+4]= 0.9*ww[offset+1]*ww[offset+3];flag=1;}
        if (ww[offset+5] < 0.9*ww[offset+2]*ww[offset+3] )  {ww[offset+5]= 0.9*ww[offset+2]*ww[offset+3];flag=1;}

        offset += 6;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        //if(fabs(yy[offset])>0.5)               {ww[offset] = xx[offset] - dsign(yy[offset])*0.5;flag=1;}
	if (ww[offset+1]<0.7*ps->LatticeTemp ) {ww[offset+1]=0.7*ps->LatticeTemp;flag=1;}
        offset += 2;
      }
      offset += pzonedata->electrode.size();
    }
    else if(ps->zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(ps->zonedata[z]);
      for(int i=0;i<pzonedata->pzone->davcell.size();i++)
      {
        //if(fabs(yy[offset])>0.5)               {ww[offset] = xx[offset] - dsign(yy[offset])*0.5;flag=1;}
	if (ww[offset+1]<0.7*ps->LatticeTemp ) {ww[offset+1]=0.7*ps->LatticeTemp;flag=1;}
        offset += 2;
      }
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


/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::init_solver:  This function do initial setup for nonlinear solver
 */
int EBM_Solver_L3E::init_solver(SolveDefine &sv)
{
  gss_log.string_buf()<<"EBM solver Level 3E init...";
  //set Tolerances
  relative_toler           = sv.relative_toler;
  toler_relax              = sv.toler_relax;
  possion_abs_toler        = sv.possion_abs_toler;
  elec_continuty_abs_toler = sv.elec_continuty_abs_toler;
  hole_continuty_abs_toler = sv.hole_continuty_abs_toler;
  heat_equation_abs_toler  = sv.heat_equation_abs_toler;
  elec_energy_abs_toler    = sv.elec_energy_abs_toler;
  hole_energy_abs_toler    = sv.hole_energy_abs_toler;
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
      N += 6*pzonedata->pzone->davcell.size();
      N += pzonedata->electrode.size();  //additional equation for electrode
      zofs[i+1] = N;
    }
    else if(zonedata[i]->material_type==Insulator) //Insulator zone
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      N += 2*zone[i].davcell.size();
      N += pzonedata->electrode.size();  //additional equation for electrode
      zofs[i+1] = N;
    }
    else if(zonedata[i]->material_type==Conductor) //Electrode zone
    {
      N += 2*zone[i].davcell.size();
      zofs[i+1] = N;
    }
    else  //other zones, such as PML and vacuum
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
  PetscScalar init_value[6];
  int    index[6];
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
        index[3] = offset+3;
        index[4] = offset+4;
        index[5] = offset+5;
        init_value[0] = pzonedata->fs[i].P;
        init_value[1] = pzonedata->fs[i].n;
        init_value[2] = pzonedata->fs[i].p;
        init_value[3] = pzonedata->fs[i].T;
        init_value[4] = pzonedata->fs[i].n*pzonedata->fs[i].Tn;
        init_value[5] = pzonedata->fs[i].p*pzonedata->fs[i].Tp;
        VecSetValues(x,6,index,init_value,INSERT_VALUES);
        offset += 6;
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
        VecSetValue(x,offset++,pzonedata->fs[i].P,INSERT_VALUES);
        VecSetValue(x,offset++,pzonedata->fs[i].T,INSERT_VALUES);
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
        VecSetValue(x,offset++,pzonedata->fs[i].P,INSERT_VALUES);
        VecSetValue(x,offset++,pzonedata->fs[i].T,INSERT_VALUES);
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
          *p++ = 6*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 6*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 6*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 6*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 6*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 6*(pzonedata->pzone->davcell[j].nb_num+1);
        }
        else //interface node, slightly more than needed.
        {
          *p++ = 6*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 6*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 6*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 6*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 6*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 6*(2*pzonedata->pzone->davcell[j].nb_num+4);
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
      {
        *p++ = 2*(zone[i].davcell[j].nb_num+1);
        *p++ = 2*(zone[i].davcell[j].nb_num+1);
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
        *p++ = bc[pzonedata->electrode[j]].psegment->node_array.size()+1;
    }
    //Electrode zones, poisson equation and heat transfer equation
    if(zonedata[i]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      for(int j=0;j<zone[i].davcell.size();j++)
      {
        *p++ = 2*(zone[i].davcell[j].nb_num+1);
        *p++ = 2*(zone[i].davcell[j].nb_num+1);
      }
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
  SNESSetFunction(snes,r,SNES_form_function_pn_3E,this);

  //if user defined a matrix-free method
  if(mf_flg)
  {
    gss_log.string_buf()<<"Matrix-Free...";
    SNESSetJacobian(snes,J,JPrec,SNESMF_form_jacobian_pn_3E,this);
  }
  else
  {
    SNESSetJacobian(snes,J,J,SNES_form_jacobian_pn_3E,this);
  }


  // set the newton method
  SNESSetType(snes,SNESLS);  //default method
  // the maximum number of iterations
  PetscInt maxit = sv.maxit;
  if(sv.Type==EQUILIBRIUM)  maxit = 10*sv.maxit;

  if(sv.NS==LineSearch || sv.NS==Basic)
  {
    if(sv.NS==LineSearch)
      SNESLineSearchSet(snes,SNESLineSearchCubic,PETSC_NULL);
    if(sv.NS==Basic)
      SNESLineSearchSet(snes,SNESLineSearchNo,PETSC_NULL);

    if(sv.Damping==DampingBankRose)
      SNESLineSearchSetPostCheck(snes,LimitorByBankRose_L3E,this);
    else if(sv.Damping==DampingPotential)
      SNESLineSearchSetPostCheck(snes,LimitorByPotential_L3E,this);
    else
      SNESLineSearchSetPostCheck(snes,LimitorNonNegativeCarrier_L3E,this);

    SNESSetTolerances(snes,1e-12*N,1e-14,1e-9,maxit,1000);
    SNESSetConvergenceTest(snes,LineSearch_ConvergenceTest_L3E,this);
  }
  else if(sv.NS==TrustRegion)
  {
    SNESSetType(snes,SNESTR);
    SNESSetTolerances(snes,1e-12*N,1e-14,1e-9,maxit,1000);
    SNESSetTrustRegionTolerance(snes,1e-20);
    SNESSetConvergenceTest(snes,TrustRegion_ConvergenceTest_L3E,this);
    //must set TR delta0 to a sufficient large value, or TR can't find real solution.
    SNES_TR  *neP = (SNES_TR*)snes->data;
    neP->delta0 = 0.1*N;
  }

  //SNESMonitorSet(snes, SNESMonitorDefault,PETSC_NULL,PETSC_NULL);

  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp,&pc);
  //KSPMonitorSet(ksp,KSPMonitorDefault,PETSC_NULL,PETSC_NULL);
  if(sv.LS=="lu"||sv.LS=="superlu"||sv.LS=="umfpack")
  {
    KSPSetType(ksp,KSPPREONLY);
    PCSetType(pc,PCLU);
    PCFactorSetReuseOrdering(pc,PETSC_TRUE);
    PCFactorSetPivoting(pc,1.0);
    PCFactorReorderForNonzeroDiagonal(pc,1e-20);
    PCFactorSetShiftNonzero(pc,1e-20);
  }
  else
  {
    KSPSetType(ksp,sv.LS.c_str());
    if(sv.LS=="gmres")  KSPGMRESSetRestart(ksp,250);
    KSPSetTolerances(ksp,1e-14*N,1e-20*N,PETSC_DEFAULT,max(35,N/10));
    PCSetType(pc,PCILU);
    PCFactorSetLevels(pc,6);
    PCFactorReorderForNonzeroDiagonal(pc,1e-20);
    PCFactorSetShiftNonzero(pc,1e-20);
  }
  KSPSetFromOptions(ksp);
  SNESSetFromOptions(snes);

  gss_log.string_buf()<<"done\n";
  gss_log.record();
  return 0;
}

/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::do_solve:  This function solve the problem
 */
int EBM_Solver_L3E::do_solve(SolveDefine &sv)
{

  if(sv.Type==EQUILIBRIUM)
    solve_equ(sv);

  else if(sv.Type==STEADYSTATE)
    solve_steadystate(sv);

  else if(sv.Type==DCSWEEP)
    solve_dcsweep(sv);

  else if(sv.Type==TRANSIENT)
    solve_transient(sv);

  else
  {
    gss_log.string_buf()<<"EBML3E not support this solver type!\n";
    gss_log.record();
  }

  return 0;

}



/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::solve_equ:  This function compute equilibrium
 * all the stimulate source(s) are zero. time step set to inf
 */
int EBM_Solver_L3E::solve_equ(SolveDefine &sv)
{
  gss_log.string_buf()<<"Compute equilibrium"<<"\n";
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
 * EBM_Solver_L3E::solve_steadystate:  This function compute steadystate
 * set all the stimulate source(s) at there value of t=0. time step set to inf
 */
int EBM_Solver_L3E::solve_steadystate(SolveDefine &sv)
{

  gss_log.string_buf()<<"Compute steady-state"<<"\n";
  gss_log.record();
  bc.Update_Vapp(0);
  bc.Update_Iapp(0);
  ODE_formula.TimeDependent = false;
  ODE_formula.dt = 1e100;
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
 * EBM_Solver_L3E::solve_dcsweep:  This function compute IV curve
 * time step set to inf
 */
int EBM_Solver_L3E::solve_dcsweep(SolveDefine &sv)
{
  PetscScalar I=0;
  FILE *fiv;
  PetscScalar current_scale_mA = scale_unit.s_mA;
  PetscScalar voltage_scale_V =scale_unit.s_volt;

  bc.Update_Vapp(0);
  bc.Update_Iapp(0);

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
    PetscScalar VStep = sv.VStep;
    for(int i=0;i<sv.Electrode_VScan_Name.size();i++)
      bc.Set_electrode_type(sv.Electrode_VScan_Name[i].c_str(),VoltageBC);
    int diverged_retry=0;
    probe_open(DCSWEEP_VSCAN);
    do
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
        if(fabs(VStep) < fabs(sv.VStep))
          VStep *= 1.1;
        V+=VStep;
        if(V*sv.VStep > sv.VStop*sv.VStep && V*sv.VStep < (sv.VStop + VStep - 1e-10*VStep)*sv.VStep)
          V=sv.VStop;
      }
      else // oh, diverged... reduce step and try again
      {
        if(++diverged_retry>=8)
        {
          gss_log.string_buf()<<"------> Too many failed steps, give up tring.\n\n\n";
          gss_log.record();
          break;
        }
        diverged_recovery();
        VStep /= 2.0;
        V-=VStep;
        gss_log.string_buf()<<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        gss_log.record();
        continue;
      }

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
    }
    while(V*sv.VStep < (sv.VStop + 0.5*VStep)*sv.VStep);
  }

  if(sv.Electrode_IScan!=-1)
  {
    bc.Set_electrode_type(sv.Electrode_IScan_Name.c_str(),CurrentBC);
    PetscScalar I = sv.IStart;
    PetscScalar IStep = sv.IStep;
    int diverged_retry=0;
    probe_open(DCSWEEP_ISCAN);
    do
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
        if(fabs(IStep) < fabs(sv.IStep))
          IStep *= 1.1;
        I+=IStep;
        if(I*sv.IStep > sv.IStop*sv.IStep && I*sv.IStep < (sv.IStop + IStep - 1e-10*IStep)*sv.IStep)
          I=sv.IStop;
      }
      else // oh, diverged... reduce step and try again
      {
        if(++diverged_retry>=8) //failed 8 times, stop tring
        {
          gss_log.string_buf()<<"------> Too many failed steps, give up tring.\n\n\n";
          gss_log.record();
          break;
        }
        diverged_recovery();
        IStep /= 2.0;
        I-=IStep;
        gss_log.string_buf()<<"------> nonlinear solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        gss_log.record();
        continue;
      }

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
    while(I*sv.IStep < (sv.IStop+0.5*IStep)*sv.IStep);
  }

  if(!sv.IVFile.empty())        fclose(fiv);
  probe_close();
  return 0;

}


void EBM_Solver_L3E::LET_norm_estimat(PetscScalar & r)
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
        ll[offset+3]=ll[offset+3]/(eps_r*xx[offset+3]+eps_a);
        ll[offset+4]=ll[offset+4]/(eps_r*xx[offset+4]+eps_a);
        ll[offset+5]=ll[offset+5]/(eps_r*xx[offset+5]+eps_a);
        offset += 6;
        N+=5;
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
        ll[offset+0]=0;
        ll[offset+1]=ll[offset+1]/(eps_r*xx[offset+1]+eps_a);
        offset += 2;
        N+=1;
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
        ll[offset+0]=0;
        ll[offset+1]=ll[offset+1]/(eps_r*xx[offset+1]+eps_a);
        offset += 2;
        N+=1;
      }
    }
  }
  VecRestoreArray(x,&xx);
  VecRestoreArray(LTE,&ll);
  VecNorm(LTE,NORM_2,&r);

  r/=sqrt((PetscScalar)N);
}


/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::solve_transient:  This function does transient simulation.
 */
int EBM_Solver_L3E::solve_transient(SolveDefine &sv)
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


  }
  while(clock<sv.TStop+0.5*ODE_formula.dt);

  //close record file
  if(!sv.IVFile.empty())
    fclose(fiv);
  probe_close();

  return 0;

}



/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::solver_pre_compute:  This function precompute some parameters before
 * each newton iteration.
 */
void EBM_Solver_L3E::solver_pre_compute()
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
        pzonedata->aux[i].affinity = pzonedata->mt->basic->Affinity(pfs->T);
        pzonedata->aux[i].density =  pzonedata->mt->basic->Density(pfs->T);
        pzonedata->aux[i].Eg     =   pzonedata->mt->band->Eg(pfs->T);
        pzonedata->aux[i].Nc     =   pzonedata->mt->band->Nc(pfs->T);
        pzonedata->aux[i].Nv     =   pzonedata->mt->band->Nv(pfs->T);
      }
    }
  }
}


/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::solution_update:  This function restore solution data from SNES
 */
void EBM_Solver_L3E::solution_update()
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
        pzonedata->fs[i].T_last = pzonedata->fs[i].T;
        pzonedata->fs[i].Tn_last = pzonedata->fs[i].Tn;
        pzonedata->fs[i].Tp_last = pzonedata->fs[i].Tp;
        pzonedata->fs[i].P = xx[offset+0];
        pzonedata->fs[i].n = fabs(xx[offset+1]);
        pzonedata->fs[i].p = fabs(xx[offset+2]);
        pzonedata->fs[i].T = xx[offset+3];
        pzonedata->fs[i].Tn = xx[offset+4]/pzonedata->fs[i].n;
        pzonedata->fs[i].Tp = xx[offset+5]/pzonedata->fs[i].p;

        pzonedata->aux[i].affinity = pzonedata->mt->basic->Affinity(xx[offset+3]);
        pzonedata->aux[i].Eg = pzonedata->mt->band->Eg(xx[offset+3]);
        pzonedata->aux[i].Nc = pzonedata->mt->band->Nc(xx[offset+3]);
        pzonedata->aux[i].Nv = pzonedata->mt->band->Nv(xx[offset+3]);
        pzonedata->mt->mapping(&pzonedata->pzone->danode[i],&pzonedata->aux[i],0);
        PetscScalar nie = pzonedata->mt->band->nie(pzonedata->fs[i].T);
        pzonedata->aux[i].phi_intrinsic = pzonedata->fs[i].P + pzonedata->aux[i].affinity +
                                          kb*pzonedata->fs[i].T/e*log(pzonedata->aux[i].Nc/nie);
        pzonedata->aux[i].phin = pzonedata->aux[i].phi_intrinsic - log(fabs(pzonedata->fs[i].n)/nie)*kb*pzonedata->fs[i].T/e;
        pzonedata->aux[i].phip = pzonedata->aux[i].phi_intrinsic + log(fabs(pzonedata->fs[i].p)/nie)*kb*pzonedata->fs[i].T/e;

        offset += 6;
      }
      pzonedata->F3E_efield_update(xx,zofs,bc,zonedata);
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar P = xx[offset];
        PetscScalar Pn = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        PetscScalar C =  bc.Get_pointer(pzonedata->electrode[i])->Get_C();
        PetscScalar I =  bc.Get_pointer(pzonedata->electrode[i])->Get_Current_new();
        bc.Get_pointer(pzonedata->electrode[i])->Set_Cap_Current(C*(P-Pn)/ODE_formula.dt);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current(I);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Potential(P);
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      pzonedata->F3_efield_update(xx,zofs,bc,zonedata);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        pzonedata->fs[i].P_last = pzonedata->fs[i].P;
        pzonedata->fs[i].T_last = pzonedata->fs[i].T;
        pzonedata->fs[i].P = xx[offset+0];
        pzonedata->fs[i].T = xx[offset+1];
        offset += 2;
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar P = xx[offset];
        PetscScalar Pn = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        PetscScalar C =  bc.Get_pointer(pzonedata->electrode[i])->Get_C();
        PetscScalar I =  bc.Get_pointer(pzonedata->electrode[i])->Get_Current_new();
        bc.Get_pointer(pzonedata->electrode[i])->Set_Cap_Current(C*(P-Pn)/ODE_formula.dt);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current(I);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Potential(P);
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        pzonedata->fs[i].P_last = pzonedata->fs[i].P;
        pzonedata->fs[i].T_last = pzonedata->fs[i].T;
        pzonedata->fs[i].P = xx[offset+0];
        pzonedata->fs[i].T = xx[offset+1];
        offset += 2;
      }
    }
  }
  ODE_formula.dt_last = ODE_formula.dt;

  VecRestoreArray(x,&xx);
}

/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::diverged_recovery:  This function recovery latest solution data
 * if SNES diverged.
 */
void EBM_Solver_L3E::diverged_recovery()
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
        xx[offset+3] = pzonedata->fs[i].T;
        xx[offset+4] = pzonedata->fs[i].n*pzonedata->fs[i].Tn;
        xx[offset+5] = pzonedata->fs[i].p*pzonedata->fs[i].Tp;
        offset += 6;
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
        xx[offset+0] = pzonedata->fs[i].P;
        xx[offset+1] = pzonedata->fs[i].T;
        offset += 2;
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
        xx[offset+0] = pzonedata->fs[i].P;
        xx[offset+1] = pzonedata->fs[i].T;
        offset += 2;
      }
    }
  }
  VecRestoreArray(x,&xx);
}



/* ----------------------------------------------------------------------------
 * EBM_Solver_L3E::destroy_solver:  This function do destroy the nonlinear solver
 */
int EBM_Solver_L3E::destroy_solver(SolveDefine &sv)
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

