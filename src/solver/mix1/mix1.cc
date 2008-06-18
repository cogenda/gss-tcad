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
/* the interface data structure with ng-spice */

#include "mathfunc.h"
#include "mix1.h"
#include "private/kspimpl.h"
#include "private/snesimpl.h"
#include "src/snes/impls/tr/tr.h"
#include "log.h"
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

/* ----------------------------------------------------------------------------
 * error_norm_pn_Mix1:  This function compute X and RHS error norms. 
 */
void DDM_Mix_Solver_L1E:: error_norm_pn_Mix1(PetscScalar *x,PetscScalar *f)
{
  // do clear
  potential_norm = 0;
  electron_norm  = 0;
  hole_norm      = 0;

  possion_norm        = 0;
  elec_continuty_norm = 0;
  hole_continuty_norm = 0;

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
}


/* ----------------------------------------------------------------------------
 * Convergence Test function for line search method.
 */
PetscErrorCode LineSearch_ConvergenceTest_Mix1(SNES snes,PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm,
    SNESConvergedReason *reason,void *dummy)
{
  DDM_Mix_Solver_L1E *ps = (DDM_Mix_Solver_L1E*)dummy;

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
  
  if (fnorm != fnorm)
  {
    *reason = SNES_DIVERGED_FNORM_NAN;
  }
  else if (ps->possion_norm < ps->possion_abs_toler &&
           ps->elec_continuty_norm < ps->elec_continuty_abs_toler &&
           ps->hole_continuty_norm < ps->hole_continuty_abs_toler )
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
             ps->hole_continuty_norm < ps->toler_relax*ps->hole_continuty_abs_toler )
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


/* ----------------------------------------------------------------------------
 * Convergence Test function for trust region method.
 */
PetscErrorCode TrustRegion_ConvergenceTest_Mix1(SNES snes,PetscInt it,PetscReal xnorm, PetscReal pnorm, PetscReal fnorm,
    SNESConvergedReason *reason,void *dummy)
{
  SNES_TR *neP = (SNES_TR *)snes->data;
  DDM_Mix_Solver_L1E *ps = (DDM_Mix_Solver_L1E*)dummy;

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
        ps->hole_continuty_norm < ps->toler_relax*ps->hole_continuty_abs_toler )
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
           ps->hole_continuty_norm < ps->hole_continuty_abs_toler )
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
             ps->hole_continuty_norm < ps->toler_relax*ps->hole_continuty_abs_toler )
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



/* ----------------------------------------------------------------------------
 * form_function_pn:  This function setup DDM equation F(x)=0
 */
void DDM_Mix_Solver_L1E::form_function_pn_Mix1(PetscScalar *x,PetscScalar *f)
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
      // compute electrode current.
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==OhmicContact)
          pzonedata->F1E_mix_om_electrode_current(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==SchottkyContact)
          pzonedata->F1E_mix_stk_electrode_current(i,x,f,ODE_formula,zofs,bc,DeviceDepth);  
        if(bc[pzonedata->electrode[i]].BCType==InsulatorContact)
          pzonedata->F1E_mix_ins_electrode_current(i,x,f,ODE_formula,zofs,bc,DeviceDepth);    
      }
      // process cell variables and boundaries.
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index||bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F1E_ddm_inner(i,x,f,ODE_formula,zofs);
        //process the electrode bc for mix simulation.  
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
          pzonedata->F1E_mix_ddm_ombc(i,x,f,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
          pzonedata->F1E_mix_ddm_stkbc(i,x,f,ODE_formula,zofs,bc);  
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorContact)
          pzonedata->F1E_mix_ddm_insulator_gate(i,x,f,ODE_formula,zofs,bc); 
        // the routine for interface BC do not change, share the DDML1E routine.     
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
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index||bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->F1_ddm_inner(i,x,f,ODE_formula,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==GateContact)
          pzonedata->F1_mix_ddm_gatebc(i,x,f,ODE_formula,zofs,bc);
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
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
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
        //compute gate electrode current
        if(bc[pzonedata->electrode[i]].BCType==GateContact)
          pzonedata->F1_mix_gate_electrode_current(i,x,f,ODE_formula,zofs,bc,DeviceDepth);
        //addtional equ for charge electrode
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
void DDM_Mix_Solver_L1E::form_jacobian_pn_Mix1(PetscScalar *x, Mat *jac, Mat *jtmp)
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

      // process cell variables and boundaries.
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        if(!pcell->bc_index||bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->J1E_ddm_inner(i,x,jac,&JTmp,ODE_formula,zofs);
        //process the electrode bc for mix simulation.  
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
          pzonedata->J1E_mix_ddm_ombc(i,x,jac,&JTmp,ODE_formula,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
          pzonedata->J1E_mix_ddm_stkbc(i,x,jac,&JTmp,ODE_formula,zofs,bc);  
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorContact)
          pzonedata->J1E_mix_ddm_insulator_gate(i,x,jac,&JTmp,ODE_formula,zofs,bc); 
        // the routine for interface BC do not change, share the DDML1E routine.  
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
          pzonedata->J1_mix_ddm_gatebc(i,x,jac,ODE_formula,zofs,bc);
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
        //addtional equ for charge electrode
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
 * SNES_form_function_pn_Mix1:  wrap function for petsc nonlinear solver
 */
PetscErrorCode SNES_form_function_pn_Mix1(SNES snes, Vec x,Vec f,void *dummy)
{
  PetscScalar    *xx,*ff;
  DDM_Mix_Solver_L1E *ps = (DDM_Mix_Solver_L1E*)dummy;
  VecZeroEntries(f);
  //Get pointers to vector data.
  VecGetArray(x,&xx);
  VecGetArray(f,&ff);

  ps->form_function_pn_Mix1(xx,ff);
  ps->error_norm_pn_Mix1(xx,ff);

  //Restore vectors
  VecRestoreArray(x,&xx);
  VecRestoreArray(f,&ff);
  VecNorm(f,NORM_2,&ps->norm);
  return 0;
}


/* ----------------------------------------------------------------------------
 * SNES_form_jacobian_pn_Mix1:  wrap function for petsc nonlinear solver
 */
PetscErrorCode SNES_form_jacobian_pn_Mix1(SNES snes, Vec x,Mat *jac,Mat *B,MatStructure *flag,void *dummy)
{
  PetscScalar    *xx;
  DDM_Mix_Solver_L1E *ps = (DDM_Mix_Solver_L1E*)dummy;
  //Get pointer to vector data
  VecGetArray(x,&xx);
  //clear old matrix
  MatZeroEntries(*jac);
  MatZeroEntries(ps->JTmp);
  //build matrix here
  ps->form_jacobian_pn_Mix1(xx,jac,&ps->JTmp);
  *flag = SAME_NONZERO_PATTERN;
  //Restore vector
  VecRestoreArray(x,&xx);
  //Assemble matrix
  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  //MatView(*jac,PETSC_VIEWER_DRAW_WORLD);
  return 0;
}

PetscErrorCode LimitorByBankRose_Mix1(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  int it;
  DDM_Mix_Solver_L1E *ps = (DDM_Mix_Solver_L1E*)checkctx;
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


PetscErrorCode LimitorByPotential_Mix1(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;
  VecGetArray(x,&xx);
  VecGetArray(y,&yy);
  VecGetArray(w,&ww);
  PetscScalar dV_max=0;
  
  //search for dV_max;
  DDM_Mix_Solver_L1E *ps = (DDM_Mix_Solver_L1E*)checkctx;
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


PetscErrorCode LimitorNonNegativeCarrier_Mix1(SNES snes, Vec x,Vec y,Vec w,void *checkctx, PetscTruth *changed_y,PetscTruth *changed_w)
{
  PetscScalar    *xx;
  PetscScalar    *yy;
  PetscScalar    *ww;
  VecGetArray(x,&xx);
  VecGetArray(y,&yy);
  VecGetArray(w,&ww);
  
  DDM_Mix_Solver_L1E *ps = (DDM_Mix_Solver_L1E*)checkctx;
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


/* ----------------------------------------------------------------------------
 * DDM_Mix_Solver_L1E::init_solver:  This function do initial setup for nonlinear solver
 */
int DDM_Mix_Solver_L1E::init_solver(SolveDefine &sv)
{
  gss_log.string_buf()<<"DDM Level 1E Device/Circuit Mix Solver init...\n";
  gss_log.record();
  psv = &sv;
  //connect to ngspice
  if(NDEVServer(sv.port,listener,client))
  {
    gss_log.string_buf()<<"I can't connect to ngspice, mix simulation stopped.\n";
    gss_log.record();
    return 1;
  }
  //get ngspice terminal information
  recv(client,&(Deviceinfo),sizeof(sDeviceinfo),MSG_FLAG);
  gss_log.string_buf()<<"Ngspice demands device "<<Deviceinfo.NDEVname<<" to be calculated.\n";
  gss_log.record();
  for(int i=0;i<Deviceinfo.term;i++)
  {
    recv(client,&(PINinfos[i]),sizeof(sPINinfo),MSG_FLAG);
    if(!bc.Get_electrode_pointer_nocase(PINinfos[i].name))
    {
      gss_log.string_buf()<<"I can't link terminal name "<<PINinfos[i].name<<" to GSS eletrode boundary.\n";
      gss_log.record();
      return 1;
    }
    gss_log.string_buf()<<"   Terminal "<<i<<": "<<PINinfos[i].name<<"\n";
    gss_log.record();
  }

  //set Tolerances
  relative_toler           = sv.relative_toler;
  toler_relax              = sv.toler_relax;
  possion_abs_toler        = sv.possion_abs_toler;
  elec_continuty_abs_toler = sv.elec_continuty_abs_toler;
  hole_continuty_abs_toler = sv.hole_continuty_abs_toler;
  
  // compute the scale of problem. 
  zofs.resize(zone_num+1);
  zofs[0] = 0;
  for(int i=0;i<zone_num;i++)
  {
    if(zonedata[i]->material_type==Semiconductor) //semiconductor zone
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      N += 3*pzonedata->pzone->davcell.size();
      //NO extra electrode equation needs.
      zofs[i+1] = N;
    }
    else if(zonedata[i]->material_type==Insulator) //Insulator zone
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      N += zone[i].davcell.size();
      N += pzonedata->electrode.size();  //additional equation for charged electrode
      zofs[i+1] = N;
    }
     else if(zonedata[i]->material_type==Conductor) //Electrode zone
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      N += zone[i].davcell.size();
      zofs[i+1] = N;
    }
    else  //other zones, such as PML and vacuum 
    {
      zofs[i+1] = N;
    }
  }

  VecCreateSeq(PETSC_COMM_SELF,N,&x);
  VecDuplicate(x,&r);
  for(int i=0;i<Deviceinfo.term;i++)
  {
    VecDuplicate(x,&PINconds[i].pdI_pdw);
    VecDuplicate(x,&PINconds[i].pdF_pdV);
    VecDuplicate(x,&PINconds[i].pdw_pdV);
  }

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
        VecSetValue(x,offset,0,INSERT_VALUES);
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
  delete [] nnz;

  SNESCreate(PETSC_COMM_WORLD,&snes);
  SNESSetFunction(snes,r,SNES_form_function_pn_Mix1,this);
  SNESSetJacobian(snes,J,J,SNES_form_jacobian_pn_Mix1,this);

  // set the newton method
  SNESSetType(snes,SNESLS);  //default method

  // the maximum number of iterations
  PetscInt maxit = sv.maxit;

  SNESSetType(snes,SNESLS);  //default method
  if(sv.NS==LineSearch || sv.NS==Basic)
  {
    if(sv.NS==LineSearch)
      SNESLineSearchSet(snes,SNESLineSearchCubic,PETSC_NULL);
    if(sv.NS==Basic)
      SNESLineSearchSet(snes,SNESLineSearchNo,PETSC_NULL);
    
    if(sv.Damping==DampingBankRose)
       SNESLineSearchSetPostCheck(snes,LimitorByBankRose_Mix1,this);
    else if(sv.Damping==DampingPotential)
       SNESLineSearchSetPostCheck(snes,LimitorByPotential_Mix1,this);
    else
       SNESLineSearchSetPostCheck(snes,LimitorNonNegativeCarrier_Mix1,this);
    
    SNESSetTolerances(snes,1e-12*N,1e-14,1e-9,maxit,1000);
    SNESSetConvergenceTest(snes,LineSearch_ConvergenceTest_Mix1,this);
  }
  else if(sv.NS==TrustRegion)
  {
    SNESSetType(snes,SNESTR);
    SNESSetTolerances(snes,1e-12*N,1e-14,1e-9,maxit,1000);
    SNESSetTrustRegionTolerance(snes,1e-30);
    SNESSetConvergenceTest(snes,TrustRegion_ConvergenceTest_Mix1,this);
    //must set TR delta0 to a sufficient large value, or TR can't find real solution.
    SNES_TR  *neP = (SNES_TR*)snes->data;
    neP->delta0 = N;
  }
  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp,&pc);

  if(sv.LS=="lu"||sv.LS=="superlu"||sv.LS=="umfpack")
  {
    KSPSetType(ksp,KSPPREONLY);
    PCSetType(pc,PCLU);
    PCFactorSetReuseOrdering(pc,PETSC_TRUE);
    PCFactorSetReuseFill(pc,PETSC_TRUE);
    PCFactorSetPivoting(pc,1.0);
    PCFactorReorderForNonzeroDiagonal(pc,1e-14);
    PCFactorSetShiftNonzero(pc,1e-12);
  }
  else
  {
    KSPSetType(ksp,sv.LS.c_str());
    if(sv.LS=="gmres")  KSPGMRESSetRestart(ksp,150);
    KSPSetTolerances(ksp,1e-14*N,1e-25*N,PETSC_DEFAULT,max(35,N/10));
    PCSetType(pc,PCILU);
    PCFactorSetLevels(pc,3);
    PCFactorSetShiftNonzero(pc,1e-14);
    PCFactorReorderForNonzeroDiagonal(pc,1e-12);
  }

  KSPSetFromOptions(ksp);
  SNESSetFromOptions(snes);
  gss_log.string_buf()<<"DDM Level 1E Device/Circuit Mix Solver init ok.\n";
  gss_log.record();
  return 0;
}


/* ----------------------------------------------------------------------------
 * DDM_Mix_Solver_L1E::do_solve:  This function receive ngspice device information.
 */
int DDM_Mix_Solver_L1E::do_solve(SolveDefine &sv)
{
  ODE_formula.clock = -1e100;
  int bytes;
  for(;;)
  {
    if (sigsetjmp(net_buf, 1) == 1) break;
    //receive solver information
    bytes=recv(client,&(CKTInfo),sizeof(sCKTinfo),MSG_FLAG);
    //if ngspice terminated, stop loop. 
    if(bytes<sizeof(sCKTinfo)) 
    {
      gss_log.string_buf()<<"Lost connection to remote ngspice.\n";
      gss_log.record();
      break;
    }
    switch(CKTInfo.DEV_CALL)
    {
     case   NDEV_LOAD                : DEV_LOAD();   break;
     case   NDEV_ACCEPT              : DEV_ACCEPT(); break;
     case   NDEV_CONVERGINCE_TEST    : DEV_CONV_TEST(); break;
     default : break;
    }
  }
  return 0;
}


int  DDM_Mix_Solver_L1E::DEV_LOAD()
{
    PetscScalar    *xx;
    //receive terminal voltage of device
    for(int i=0;i<Deviceinfo.term;i++)
    {
      recv(client,&PINinfos[i],sizeof(sPINinfo),MSG_FLAG);
      bc.Set_Vapp_nocase(PINinfos[i].name,PINinfos[i].V);
    }

    //transient calculation settings
    if(CKTInfo.CKTmode & MODETRAN)
    {
      ODE_formula.TimeDependent = true;
      ODE_formula.BDF_Type = BDF2;
      if(CKTInfo.CKTmode & MODEINITTRAN) //the first step of transient simulation
      {
        gss_log.string_buf()<<"----------------------------------------------------------------------\n";
        gss_log.string_buf()<<"NGSPICE Transient Mode Start   Time = "<<CKTInfo.time<<"s  dt = "<<CKTInfo.dt<<"s\n";
        for(int i=0;i<Deviceinfo.term;i++)
        {
          gss_log.string_buf()<<"Set V("<<PINinfos[i].name<<") = "<<PINinfos[i].V<<"V\n";
          gss_log.record();
        }
        gss_log.string_buf()<<"----------------------------------------------------------------------\n\n";
        gss_log.record();
        ODE_formula.BDF2_restart = true;
      }
      else if(CKTInfo.CKTmode & MODEINITPRED) //the time matching flag
      {
        gss_log.string_buf()<<"----------------------------------------------------------------------\n";
        gss_log.string_buf()<<"NGSPICE Transient Mode Update  Time = "<<CKTInfo.time<<"s  dt = "<<CKTInfo.dt<<"s\n";
        for(int i=0;i<Deviceinfo.term;i++)
        {
          gss_log.string_buf()<<"Set V("<<PINinfos[i].name<<") = "<<PINinfos[i].V<<"V\n";
          gss_log.record();
        }
        if(ODE_formula.clock > CKTInfo.time*scale_unit.s_second)
        {
          time_back_recovery();
          ODE_formula.BDF2_restart = true;
          gss_log.string_buf()<<"NGSPICE back trace...\n";
	  gss_log.record();
        }
        else
        {
          ODE_formula.BDF2_restart = false;
        }
	gss_log.string_buf()<<"----------------------------------------------------------------------\n\n";
        gss_log.record();
      }
      ODE_formula.clock = CKTInfo.time*scale_unit.s_second;
      // call the real computation routine
      tran_solve();
    }
    //DC calculation settings
    else if(CKTInfo.CKTmode & MODEDC)
    {
      gss_log.string_buf()<<"----------------------------------------------------------------------\n";
      gss_log.string_buf()<<"NGSPICE DC Mode\n";
      for(int i=0;i<Deviceinfo.term;i++)
      {
          gss_log.string_buf()<<"Set V("<<PINinfos[i].name<<") = "<<PINinfos[i].V<<"V\n";
          gss_log.record();
      }
      gss_log.string_buf()<<"----------------------------------------------------------------------\n\n";
      gss_log.record();
      ODE_formula.TimeDependent = false;
      ODE_formula.dt = 1e100;
      // call the real computation routine
      dc_solve();
    }
    
    //get pdI/pdw, pdI/pdV and pdF/pdV for each electrode
    VecGetArray(x,&xx);
    for(int i=0;i<Deviceinfo.term;i++)
      for(int z=0;z<zone_num;z++)
    {
        if(zonedata[z]->material_type==Semiconductor)
        {
          SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
          for(int j=0;j<pzonedata->electrode.size();j++)
             if(bc.is_electrode_label_nocase(pzonedata->electrode[j], PINinfos[i].name))
            { 
              if(bc[pzonedata->electrode[j]].BCType==OhmicContact)
                pzonedata->F1E_mix_om_electrode_Load(pzonedata->electrode[j],PINinfos[i].I,PINconds[i].pdI_pdV,
                                               PINconds[i].pdI_pdw,PINconds[i].pdF_pdV,
                                               xx,&JTmp,ODE_formula,zofs,bc,DeviceDepth);
              if(bc[pzonedata->electrode[j]].BCType==SchottkyContact)
                pzonedata->F1E_mix_stk_electrode_Load(pzonedata->electrode[j],PINinfos[i].I,PINconds[i].pdI_pdV,
                                               PINconds[i].pdI_pdw,PINconds[i].pdF_pdV,
                                               xx,&JTmp,ODE_formula,zofs,bc,DeviceDepth);
              if(bc[pzonedata->electrode[j]].BCType==InsulatorContact)
                pzonedata->F1E_mix_ins_electrode_Load(pzonedata->electrode[j],PINinfos[i].I,PINconds[i].pdI_pdV,
                                               PINconds[i].pdI_pdw,PINconds[i].pdF_pdV,
                                               xx,&JTmp,ODE_formula,zofs,bc,DeviceDepth);                              
            }
        }
        //
        if(zonedata[z]->material_type==Insulator)
        {
          ISZone *pzonedata= dynamic_cast< ISZone * >(zonedata[z]);
          for(int j=0;j<pzonedata->electrode.size();j++)
            if(bc.is_electrode_label_nocase(pzonedata->electrode[j], PINinfos[i].name))
           {
             pzonedata->F1_mix_gate_electrode_Load(pzonedata->electrode[j],PINinfos[i].I,PINconds[i].pdI_pdV,
                                               PINconds[i].pdI_pdw,PINconds[i].pdF_pdV,
                                               xx,&J,ODE_formula,zofs,bc,DeviceDepth);  
          
           }  
        }
    }
    VecRestoreArray(x,&xx);

    //save conductance and rhs current source.
    for(int i=0;i<Deviceinfo.term;i++)
    {
      PetscScalar Ig=0;
      for(int j=0;j<Deviceinfo.term;j++)
      {
        KSPSolve(ksp,PINconds[j].pdF_pdV,PINconds[i].pdw_pdV);
        PetscScalar pdI_pdV;
        VecDot(PINconds[i].pdI_pdw,PINconds[i].pdw_pdV,&pdI_pdV);
        PINinfos[i].dI_dV[j] = pdI_pdV;
        PINinfos[i].dI_dV[j] += PINconds[i].pdI_pdV;
        PINinfos[i].dI_dV[j] /= scale_unit.s_A; //scale the element of dI/dV.
        Ig += PINinfos[i].dI_dV[j]*PINinfos[j].V;
      }
      PINinfos[i].I = Ig - PINinfos[i].I/scale_unit.s_A;
    }
    
    //sent conductance matrix and rhs current back to ngspice.
    for(int i=0;i<Deviceinfo.term;i++)
      send(client,&PINinfos[i],sizeof(sPINinfo),0);
   
    return 0;  
}


int  DDM_Mix_Solver_L1E::DEV_ACCEPT()
{
  solution_update();
  return 0;
}


int  DDM_Mix_Solver_L1E::DEV_CONV_TEST()
{
  CKTInfo.convergence_flag = reason;
  send(client,&(CKTInfo),sizeof(sCKTinfo),0);
  return 0;
}


int DDM_Mix_Solver_L1E::tran_solve()
{
  PetscInt rework = 1;
  PetscInt max_rework = 64;
  PetscInt Converged = 0;
  
  for(int i=0;i<zone_num;i++)
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->HighFieldMobility    = psv->HighFieldMobility;
      pzonedata->ImpactIonization     = psv->ImpactIonization;
      pzonedata->IIType               = psv->IIType;
      pzonedata->BandBandTunneling    = psv->BandBandTunneling;
      pzonedata->IncompleteIonization = psv->IncompleteIonization;
      pzonedata->QuantumMechanical    = psv->QuantumMechanical;
      pzonedata->Fermi                = psv->Fermi;
      pzonedata->EJModel              = psv->EJModel;
    }
  
  do
  {
    reason = SNES_CONVERGED_ITERATING;
    for(int w=1;w<=rework;w++)
    {
      ODE_formula.dt = w*CKTInfo.dt*scale_unit.s_second/rework;
      
      for(int i=0;i<Deviceinfo.term;i++)
      {
        PetscScalar V_step = PINinfos[i].V_old+(PINinfos[i].V-PINinfos[i].V_old)*w/rework;
        bc.Set_Vapp_nocase(PINinfos[i].name,V_step);
      }
      
      SNESSolve(snes,PETSC_NULL,x);
      SNESGetConvergedReason(snes,&reason);
      if(reason<0)
      {
        gss_log.string_buf()<<"------> GSS mixed solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        gss_log.record();
        diverged_recovery();
        rework*=2;
        break;
      }
    }
    if(reason>0) {Converged = 1; break;}
    if(rework>max_rework) {Converged = 0; break;}
  }  while(1);
  
  return Converged;
}


int DDM_Mix_Solver_L1E::dc_solve()
{
  PetscInt rework = 1;
  PetscInt max_rework = 64;
  PetscInt Converged = 0;
    
  for(int i=0;i<zone_num;i++)
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->HighFieldMobility    = psv->HighFieldMobility;
      pzonedata->ImpactIonization     = psv->ImpactIonization;
      pzonedata->IIType               = psv->IIType;
      pzonedata->BandBandTunneling    = psv->BandBandTunneling;
      pzonedata->IncompleteIonization = psv->IncompleteIonization;
      pzonedata->QuantumMechanical    = psv->QuantumMechanical;
      pzonedata->Fermi                = psv->Fermi;
      pzonedata->EJModel              = psv->EJModel;
    }
  
  do
  {
    reason = SNES_CONVERGED_ITERATING;
    for(int w=1;w<=rework;w++)
    {
      for(int i=0;i<Deviceinfo.term;i++)
      {
        PetscScalar V_step = PINinfos[i].V_old+(PINinfos[i].V-PINinfos[i].V_old)*w/rework;
        bc.Set_Vapp_nocase(PINinfos[i].name,V_step);
      }
      
      SNESSolve(snes,PETSC_NULL,x);
      SNESGetConvergedReason(snes,&reason);
      if(reason<0)
      {
        gss_log.string_buf()<<"------> GSS mixed solver "<<SNESConvergedReasons[reason]<<", do recovery...\n\n\n";
        gss_log.record();
        diverged_recovery();
        rework*=2;
        break;
      }
    }
    if(reason>0) {Converged = 1; break;}
    if(rework>max_rework) {Converged = 0; break;}
  }  while(1);
  
  return Converged;
}

    
/* ----------------------------------------------------------------------------
 * DDM_Mix_Solver_L1E::solution_update:  This function restore solution data from SNES
 */
void DDM_Mix_Solver_L1E::solution_update()
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
        pzonedata->fs[i].n = xx[offset+1];
        pzonedata->fs[i].p = xx[offset+2];

        pzonedata->mt->mapping(&pzonedata->pzone->danode[i],&pzonedata->aux[i],0);
        PetscScalar nie = pzonedata->mt->band->nie(pzonedata->fs[i].T);
        pzonedata->aux[i].phi_intrinsic = pzonedata->fs[i].P + pzonedata->aux[i].affinity +
                                          kb*pzonedata->fs[i].T/e*log(pzonedata->aux[i].Nc/nie);
        pzonedata->aux[i].phin = pzonedata->aux[i].phi_intrinsic - log(fabs(pzonedata->fs[i].n)/nie)*kb*pzonedata->fs[i].T/e;
        pzonedata->aux[i].phip = pzonedata->aux[i].phi_intrinsic + log(fabs(pzonedata->fs[i].p)/nie)*kb*pzonedata->fs[i].T/e;

        offset += 3;
      }
      pzonedata->F1E_efield_update(xx,zofs,bc,zonedata);
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar I =  bc.Get_pointer(pzonedata->electrode[i])->Get_Current_new();
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current(I);
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
        PetscScalar I =  bc.Get_pointer(pzonedata->electrode[i])->Get_Current_new();
        bc.Get_pointer(pzonedata->electrode[i])->Set_Current(I);
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          bc.Get_pointer(pzonedata->electrode[i])->Set_Potential(xx[offset]);
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
  ODE_formula.dt_last = ODE_formula.dt;

  VecRestoreArray(x,&xx);
}


/* ----------------------------------------------------------------------------
 * DDM_Mix_Solver_L1E::time_back_recovery:  This function recovery latest solution data
 * if NGSPICE time drops back.
 */
void DDM_Mix_Solver_L1E::time_back_recovery()
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
        pzonedata->fs[i].P = pzonedata->fs[i].P_last;
        pzonedata->fs[i].n = pzonedata->fs[i].n_last;
        pzonedata->fs[i].p = pzonedata->fs[i].p_last;
        offset += 3;
      }
    }
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        xx[offset] = pzonedata->fs[i].P;
        pzonedata->fs[i].P = pzonedata->fs[i].P_last;
        offset += 1;
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          xx[offset] = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        else
          xx[offset] = 0;
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        xx[offset] = pzonedata->fs[i].P;
        pzonedata->fs[i].P = pzonedata->fs[i].P_last;
        offset += 1;
      }
    }
  }
  VecRestoreArray(x,&xx);
}



/* ----------------------------------------------------------------------------
 * DDM_Mix_Solver_L1E::diverged_recovery:  This function recovery latest solution data
 * if SNES diverged.
 */
void DDM_Mix_Solver_L1E::diverged_recovery()
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
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          xx[offset] = bc.Get_pointer(pzonedata->electrode[i])->Get_Potential();
        else
          xx[offset] = 0;
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
 * DDM_Mix_Solver_L1E::destroy_solver:  This function do destroy the nonlinear solver
 */
int DDM_Mix_Solver_L1E::destroy_solver(SolveDefine &sv)
{
  // free work space
  N = 0;
  zofs.clear();
  for(int i=0;i<Deviceinfo.term;i++)
  {
    VecDestroy(PINconds[i].pdI_pdw);
    VecDestroy(PINconds[i].pdF_pdV);
    VecDestroy(PINconds[i].pdw_pdV);
  }
  VecDestroy(x);
  VecDestroy(r);
  MatDestroy(J);
  MatDestroy(JTmp);
  SNESDestroy(snes);

  close(client);
  close(listener);
  return 0;
}

