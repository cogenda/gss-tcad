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
/*  Last update: Nov 19, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "mathfunc.h"
#include "ddm_nt1ac.h"
#include "log.h"

/* ----------------------------------------------------------------------------
 * DDM_Solver_L1AC::form_linear_system_ac1:  This function build the linear system.
 */
void DDM_Solver_L1AC::form_linear_system_ac1(PetscScalar omega, PetscScalar *x, Mat *jac, Vec *r, Mat *jtmp)
{
  MatZeroEntries(*jac);
  VecZeroEntries(*r);

  //Compute complex Jacobian matrix and right vector.
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[z]);

      // process electrode.
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==OhmicContact)
          pzonedata->AC1_om_electrode(i,omega,x,jac,r,&JTmp,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==SchottkyContact)
          pzonedata->AC1_stk_electrode(i,omega,x,jac,r,&JTmp,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==InsulatorContact)
          pzonedata->AC1_ins_electrode(i,omega,x,jac,r,&JTmp,zofs,bc,DeviceDepth);
      }

      // process cell variables and boundaries.
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        const VoronoiCell* pcell = zone[z].davcell.GetPointer(i);
        SemiData *pfs = &pzonedata->fs[i];
        if(!pcell->bc_index||bc[pcell->bc_index-1].BCType==NeumannBoundary)
          pzonedata->AC1_ddm_inner(i,omega,x,jac,r,&JTmp,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==OhmicContact)
          pzonedata->AC1_ddm_ombc(i,omega,x,jac,r,&JTmp,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==SchottkyContact)
          pzonedata->AC1_ddm_stkbc(i,omega,x,jac,r,&JTmp,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorContact)
          pzonedata->AC1_ddm_insulator_gate(i,omega,x,jac,r,&JTmp,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_interface(i,omega,x,jac,r,&JTmp,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HomoInterface)
        {
          HomoInterfaceBC *pbc;
          pbc = dynamic_cast<HomoInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_homojunction(i,omega,x,jac,r,&JTmp,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==HeteroInterface)
        {
          HeteroInterfaceBC *pbc;
          pbc = dynamic_cast<HeteroInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_heterojunction(i,omega,x,jac,r,&JTmp,zofs,bc,pz,n_node);
        }
        else
          pzonedata->AC1_ddm_inner(i,omega,x,jac,r,&JTmp,zofs);
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
          pzonedata->AC1_ddm_inner(i,omega,x,r,jac,zofs);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==GateContact)
          pzonedata->AC1_ddm_gatebc(i,omega,x,r,jac,zofs,bc);
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==ChargedContact)
          pzonedata->AC1_ddm_chargebc(i,omega,x,r,jac,zofs,bc);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
        {
          ElectrodeInsulatorBC *pbc;
          pbc = dynamic_cast<ElectrodeInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ElZone * pz = dynamic_cast<ElZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_electrode_insulator_interface(i,omega,x,r,jac,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_semiconductor_insulator_interface(i,omega,x,r,jac,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Insulator_Insulator)
        {
          InsulatorInsulatorBC *pbc;
          pbc = dynamic_cast<InsulatorInsulatorBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_insulator_insulator_interface(i,omega,x,r,jac,zofs,bc,pz,n_node);
        }
        else
          pzonedata->AC1_ddm_inner(i,omega,x,r,jac,zofs);
      }
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        if(bc[pzonedata->electrode[i]].BCType==GateContact)
          pzonedata->AC1_gate_electrode(i,omega,x,r,jac,zofs,bc,DeviceDepth);
        if(bc[pzonedata->electrode[i]].BCType==ChargedContact)
          pzonedata->AC1_charge_electrode(i,omega,x,r,jac,zofs,bc);
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
          pzonedata->AC1_ddm_inner(i,omega,x,r,jac,zofs);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==IF_Electrode_Insulator)
          pzonedata->AC1_ddm_inner(i,omega,x,r,jac,zofs);
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==OhmicContact)
        {
          OhmicBC *pbc = dynamic_cast<OhmicBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_om_contact(i,omega,x,r,jac,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==SchottkyContact)
        {
          SchottkyBC *pbc = dynamic_cast<SchottkyBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          SMCZone * pz = dynamic_cast<SMCZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_stk_contact(i,omega,x,r,jac,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==GateContact)
        {
          GateBC *pbc = dynamic_cast<GateBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_gate_contact(i,omega,x,r,jac,zofs,bc,pz,n_node);
        }
        else if(pcell->bc_index && bc[pcell->bc_index-1].BCType==ChargedContact)
        {
          ChargedContactBC *pbc = dynamic_cast<ChargedContactBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          pzonedata->AC1_ddm_charge_contact(i,omega,x,r,jac,zofs,bc,pz,n_node);
        }
        else
          pzonedata->AC1_ddm_inner(i,omega,x,r,jac,zofs);
      }
    }

  }
  MatAssemblyBegin(*jac,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*jac,MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(*r);
  VecAssemblyEnd(*r);
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1AC::Get_Jacobian:  This function get a real Jacobian matrix.
 */
void DDM_Solver_L1AC::Get_Jacobian()
{
  // compute the scale of JTmp
  zofs.resize(zone_num+1);
  zofs[0] = 0;
  N = 0;
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
      N += pzonedata->pzone->davcell.size();
      N += pzonedata->electrode.size();  //additional equation for electrode
      zofs[i+1] = N;
    }
    else if(zonedata[i]->material_type==Conductor) //Electrode zone
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      N += pzonedata->pzone->davcell.size();
      zofs[i+1] = N;
    }
    else  //other zones, such as PML and vacuum
    {
      zofs[i+1] = N;
    }
  }
  // create DC solution array v
  x = new PetscScalar[N];
  // fill initial value
  int    offset=0;
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        x[offset+0] = pzonedata->fs[i].P;
        x[offset+1] = pzonedata->fs[i].n;
        x[offset+2] = pzonedata->fs[i].p;
        offset += 3;
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        PetscScalar P = bc.Get_pointer(pzonedata->electrode[j])->Get_Potential();
        x[offset+0] = P;
        offset += 1;
      }
    }

    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        x[offset+0] = pzonedata->fs[i].P;
        offset += 1;
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        PetscScalar P = bc.Get_pointer(pzonedata->electrode[j])->Get_Potential();
        x[offset+0] = P;
        offset += 1;
      }
    }
    if(zonedata[z]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        x[offset+0] = pzonedata->fs[i].P;
        offset += 1;
      }
    }
  }
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
      {
        *p++ = bc[pzonedata->electrode[j]].psegment->node_array.size()+1;
      }

    }
    //Insulator zones, poisson equation only
    if(zonedata[i]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      for(int j=0;j<zone[i].davcell.size();j++)
      {
        *p++ = zone[i].davcell[j].nb_num+1;
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        *p++ = bc[pzonedata->electrode[j]].psegment->node_array.size()+1;
      }
    }
    //Electrode zones, poisson equation only
    if(zonedata[i]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      for(int j=0;j<zone[i].davcell.size();j++)
        *p++ = 1*(zone[i].davcell[j].nb_num+1);
    }
  }
  MatCreateSeqAIJ(PETSC_COMM_SELF,N,N,0,nnz,&JTmp);

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
        pzonedata->J1E_Tri_ddm(ptri,x,&JTmp,zofs);
      }
    }
  }
  MatAssemblyBegin(JTmp,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(JTmp,MAT_FINAL_ASSEMBLY);
  delete [] nnz;

}

/* ----------------------------------------------------------------------------
 * DDM_Solver_L1AC::init_solver:  This function do initial setup for AC linear solver
 */
int DDM_Solver_L1AC::init_solver(SolveDefine &sv)
{
  gss_log.string_buf()<<"DDM solver Level 1 AC Sweep init...";

  //distable advaced features.
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

  Get_Jacobian();

  // scale of complex Jac
  VecCreateSeq(PETSC_COMM_SELF,2*N,&v);
  VecDuplicate(v,&r);

  // Create Jacobian matrix data structure
  // pre-alloc approximate memory
  int *nnz = new int[2*N];
  int *p = nnz;
  //fill half of the array nnz
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
          *p++ = 2*3*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 2*3*(pzonedata->pzone->davcell[j].nb_num+1);
          *p++ = 2*3*(pzonedata->pzone->davcell[j].nb_num+1);
        }
        else //interface node, slightly more than needed.
        {
          *p++ = 2*3*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 2*3*(2*pzonedata->pzone->davcell[j].nb_num+4);
          *p++ = 2*3*(2*pzonedata->pzone->davcell[j].nb_num+4);
        }
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        *p++ = 2*(bc[pzonedata->electrode[j]].psegment->node_array.size()+1);
      }

    }
    //Insulator zones, poisson equation only
    if(zonedata[i]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      for(int j=0;j<zone[i].davcell.size();j++)
      {
        *p++ = 2*(zone[i].davcell[j].nb_num+1);
      }
      for(int j=0;j<pzonedata->electrode.size();j++)
      {
        *p++ = 2*(bc[pzonedata->electrode[j]].psegment->node_array.size()+1);
      }
    }
    //Electrode zones, poisson equation only
    if(zonedata[i]->material_type==Conductor)
    {
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      for(int j=0;j<zone[i].davcell.size();j++)
        *p++ = 2*(zone[i].davcell[j].nb_num+1);
    }

  }
  //fill another half of nnz array
  for(int i=0;i<N;i++)
    nnz[N+i]=nnz[i];

  MatCreate(PETSC_COMM_SELF,&J);
  MatSetSizes(J,2*N,2*N,2*N,2*N);
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

  //MatCreateSeqAIJ(PETSC_COMM_SELF,2*N,2*N,0,nnz,&J);  
  delete [] nnz;

  KSPCreate(PETSC_COMM_SELF,&ksp);
  KSPSetOperators(ksp,J,J,SAME_NONZERO_PATTERN);
  KSPGetPC(ksp,&pc);
  //KSPSetMonitor(ksp,KSPDefaultMonitor,PETSC_NULL,PETSC_NULL);
  if(sv.LS=="lu"||sv.LS=="superlu"||sv.LS=="umfpack")
  {
    KSPSetType(ksp,KSPGMRES);
    KSPSetTolerances(ksp,1e-20*N,1e-18*N,PETSC_DEFAULT,N/2);
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
    KSPSetTolerances(ksp,1e-20*N,1e-18*N,PETSC_DEFAULT,max(35,N/2));
    PCSetType(pc,PCILU);
    PCFactorSetLevels(pc,3);
    PCFactorSetShiftNonzero(pc,1e-14);
    PCFactorReorderForNonzeroDiagonal(pc,1e-12);
  }
  KSPSetFromOptions(ksp);
  gss_log.string_buf()<<"done\n";
  gss_log.record();
  return 0;
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1AC::do_solve:  This function solve the problem
 */
int DDM_Solver_L1AC::do_solve(SolveDefine &sv)
{
  if(sv.Type==ACSWEEP)
    solve_ac(sv);
  else
  {
    gss_log.string_buf()<<"DDML1AC not support this solver type!\n";
    gss_log.record();
  }
  return 0;
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1AC::solve_ac:  This function compute device response to
 * small AC signal AC.
 */
int DDM_Solver_L1AC::solve_ac(SolveDefine &sv)
{
  FILE *fiv;
  PetscScalar current_scale_mA = scale_unit.s_mA;
  PetscScalar voltage_scale_V =scale_unit.s_volt;
  KSPConvergedReason reason;
  PetscInt   its;
  PetscReal  rnorm;
  bc.Update_Vapp(0);
  bc.Update_Iapp(0);

  //set which electrodes are required to record IV
  for(int i=0;i<zone_num;i++)
  {
    if(zonedata[i]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->ImpactIonization = sv.ImpactIonization;
      pzonedata->IncompleteIonization = sv.IncompleteIonization;
      pzonedata->QuantumMechanical = sv.QuantumMechanical;
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
        {
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
        }
    }
    if(zonedata[i]->material_type == Insulator)
    {
      ISZone * pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      for(int k=0;k<sv.Electrode_Record_Name.size();k++)
        for(int j=0;j<pzonedata->electrode.size();j++)
        {
          if(bc.is_electrode_label(pzonedata->electrode[j], sv.Electrode_Record_Name[k].c_str()))
            sv.Electrode_Record.push_back(pzonedata->electrode[j]);
        }
    }
  }

  // output DC Scan information
  gss_log.string_buf()<<"AC scan from "<<sv.FStart*scale_unit.s_second/1e6
  <<" MHz step size "<<sv.FMultiple
  <<" to "<<sv.FStop*scale_unit.s_second/1e6<<" MHz"<<"\n";
  gss_log.record();

  // prepare IV file. If no file given, output to screen.
  if(!sv.IVFile.empty())
  {
    fiv=fopen(sv.IVFile.c_str(),"w");
    fprintf(fiv,"#");
    fprintf(fiv,"Freq[MHz]  ");
    for(int j=0;j<sv.Electrode_Record.size();j++)
      fprintf(fiv,"Vp(%s)[V]  I(%s)[mA]   ",
              sv.Electrode_Record_Name[j].c_str(),
              sv.Electrode_Record_Name[j].c_str()
             );
    fprintf(fiv,"\n");
  }
  else
    fiv = stdout;

  // begin
  PetscScalar omega = 2*PI*sv.FStart;
  if(omega<=0)
  {
    gss_log.string_buf()<<"The start frequency shouldn't be zero or negtive.\n";
    gss_log.record();
    return 0;
  }

  bc.Set_electrode_type(sv.Electrode_ACScan_Name.c_str(),VoltageBC);
  bc.Set_Vac(sv.Electrode_ACScan_Name.c_str(),sv.VAC);

  while( omega < 2*PI*sv.FStop)
  {
    gss_log.string_buf()<<"AC Scan: f("<<sv.Electrode_ACScan_Name<<") = "
    <<omega/(2*PI)*scale_unit.s_second/1e6<<" MHz "<<"\n";
    gss_log.record();
    form_linear_system_ac1(omega,x,&J,&r,&JTmp);
    KSPSetUp(ksp);
    KSPSolve(ksp,r,v);
    KSPGetConvergedReason(ksp,&reason);
    KSPGetIterationNumber(ksp,&its);
    KSPGetResidualNorm(ksp,&rnorm);
    gss_log.string_buf()<<"------> residual norm = "<<rnorm<<" its = "<<its<<" with "<<SNESConvergedReasons[reason]<<"\n\n";
    gss_log.record();
    compute_electrode_current(omega);
    fprintf(fiv,"%e\t",double(omega/(2*PI)*scale_unit.s_second/1e6));
    for(int j=0;j<sv.Electrode_Record.size();j++)
    {
      int zone_index = bc[sv.Electrode_Record[j]].psegment->zone_index;
      complex<PetscScalar> I = bc.Get_pointer(sv.Electrode_Record[j])->Get_Iac()/current_scale_mA;
      complex<PetscScalar> P = bc.Get_pointer(sv.Electrode_Record[j])->Get_Pac()/voltage_scale_V;
      fprintf(fiv,"%e %e %e\t %e %e %e\t",
              double(P.real()),double(P.imag()),double(abs(P)),
              double(I.real()),double(I.imag()),double(abs(I)));
    }
    fprintf(fiv,"\n");
    fflush(fiv);
    omega*=sv.FMultiple;
  }

  if(!sv.IVFile.empty())        fclose(fiv);
  return 0;
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1AC::compute_electrode_current:  This function compute electrode current after
 * AC sweep.
 */
void DDM_Solver_L1AC::compute_electrode_current(PetscScalar omega)
{
  //------------------------------------------------------------
  PetscScalar    *vv;
  VecGetArray(v,&vv);

  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      int equ_num = 3;
      int size = pzonedata->pzone->davcell.size();
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar IacRe=0,IacIm=0;
        PetscScalar PRe = vv[zofs[z]+size*equ_num+i];
        PetscScalar PIm = vv[N+zofs[z]+size*equ_num+i];

        if(bc[pzonedata->electrode[i]].BCType==OhmicContact)
          pzonedata->AC1_om_electrode_current(i,omega,x,&v,&JTmp, zofs, bc, DeviceDepth, IacRe, IacIm);
        if(bc[pzonedata->electrode[i]].BCType==SchottkyContact)
          pzonedata->AC1_stk_electrode_current(i,omega,x,&v,&JTmp, zofs, bc, DeviceDepth, IacRe, IacIm);
        if(bc[pzonedata->electrode[i]].BCType==InsulatorContact)
          pzonedata->AC1_ins_electrode_current(i,omega,x,&v,&JTmp, zofs, bc, DeviceDepth, IacRe, IacIm);
        complex<PetscScalar> I(IacRe,IacIm);
        complex<PetscScalar> P(PRe,PIm);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Iac(I);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Pac(P);

      }
    }
    if(zonedata[z]->material_type==Insulator)
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
      int equ_num = 1;
      int size = pzonedata->pzone->davcell.size();
      for(int i=0;i<pzonedata->electrode.size();i++)
      {
        PetscScalar IacRe=0,IacIm=0;
        PetscScalar PRe = vv[zofs[z]+size*equ_num+i];
        PetscScalar PIm = vv[N+zofs[z]+size*equ_num+i];
        if(bc[pzonedata->electrode[i]].BCType==GateContact)
          pzonedata->AC1_gate_electrode_current(i,omega,x,&v,zofs, bc, DeviceDepth, IacRe, IacIm);
        complex<PetscScalar> I(IacRe,IacIm);
        complex<PetscScalar> P(PRe,PIm);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Iac(I);
        bc.Get_pointer(pzonedata->electrode[i])->Set_Pac(P);
      }
    }
  }

  VecRestoreArray(v,&vv);
}


/* ----------------------------------------------------------------------------
 * DDM_Solver_L1AC::destroy_solver:  This function do destroy the nonlinear solver
 */
int DDM_Solver_L1AC::destroy_solver(SolveDefine &sv)
{
  // free work space
  N = 0;
  zofs.clear();
  delete [] x;
  VecDestroy(v);
  VecDestroy(r);
  MatDestroy(J);
  MatDestroy(JTmp);
  KSPDestroy(ksp);
  return 0;
}

