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
/*  Last update: Nov 27, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "bc.h"
#include "geom.h"
#include "material.h"
#include "log.h"

// ----------------------------------------------------------------------------
// constructors of each boundary class
// ----------------------------------------------------------------------------
NeumannBC::NeumannBC()
{
  heat_transfer = 0;
}

OhmicBC::OhmicBC()
{
  R=C=L=0;
  Vapp = 0;
  Vac = 0;
  heat_transfer = 0;
  electrode_type = VoltageBC;
  current       = 0;
  cap_current   = 0;
  potential     = 0;
  inner_connect = -1;
  pinterface    = NULL;
}

SchottkyBC::SchottkyBC()
{
  R=C=L=0;
  Vapp = 0;
  Vac = 0;
  heat_transfer = 0;
  electrode_type = VoltageBC;
  current       = 0;
  cap_current   = 0;
  potential     = 0;
  WorkFunction  = 0;
  pinterface    = NULL;
}

GateBC::GateBC()
{
  R=C=L=0;
  Vapp = 0;
  Vac = 0;
  current       = 0;
  cap_current   = 0;
  potential     = 0;
  electrode_type = VoltageBC;
  heat_transfer = 0;
  WorkFunction  = 0;
  pinterface    = NULL;
}

InsulatorContactBC::InsulatorContactBC()
{
  R=C=L=0;
  Vapp = 0;
  current       = 0;
  cap_current   = 0;
  potential     = 0;
  electrode_type = VoltageBC;
  heat_transfer = 0;
  WorkFunction  = 0;
  Thick = 1;
}

ChargedContactBC::ChargedContactBC()
{
  potential     = 0;
  QF = 0;
}

InsulatorInterfaceBC::InsulatorInterfaceBC()
{}

HeteroInterfaceBC::HeteroInterfaceBC()
{
  QF = 0;
}

AbsorbingBC::AbsorbingBC() //by zhangxih ----06-09-27
{
  heat_transfer = 0;
}

SourceBC::SourceBC() //by zhangxih ----06-09-28
{
  heat_transfer = 0;
}

// ----------------------------------------------------------------------------

DABC::DABC()
{
  bc_num = 0;
}

DABC::~DABC()
{
  clear();
}
/* ----------------------------------------------------------------------------
 * DABC::size():  This function return the number of bcs in dabc
 */
int DABC::size()
{
  return bc_num;
}

/* ----------------------------------------------------------------------------
 * DABC::clear():  This function delete all the bc in dabc
 */
void DABC::clear()
{
  for(int i=0;i<bc_point_array.size();i++)
    if(bc_point_array[i])
    {
      delete bc_point_array[i];
      bc_point_array[i]=NULL;
    }
  bc_point_array.clear();
  bc_num = 0;
}

/* ----------------------------------------------------------------------------
 * DABC::Get_pointer(i):  This function return the pointer of the ith bc
 */
BaseBC * DABC::Get_pointer(int bc_index)
{
  return        bc_point_array[bc_index];
}

/* ----------------------------------------------------------------------------
 * DABC::operator[]:  This function return the reference of the ith bc
 */
const BaseBC & DABC::operator[](int bc_index)
{
  return        *bc_point_array[bc_index];
}

/* ----------------------------------------------------------------------------
 * DABC::Get_bc_index():  This function find the bc belongs to zone
 */
int DABC::Get_bc_index(int zone_index, const char * bc_name)
{
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,bc_name) && zone_index == bc_point_array[i]->psegment->zone_index)
    {
      return i;
    }
  return -1;
}

/* ----------------------------------------------------------------------------
 * DABC::Get_bc_index():  This function find the bc belongs to zone
 */
int DABC::Get_bc_index(const char *zone_name, const char * bc_name)
{
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,bc_name) && 
       !strcmp(bc_point_array[i]->pzone->zonename,zone_name))
    {
      return i;
    }
  return -1;
}

/* ----------------------------------------------------------------------------
 * DABC::Get_bc_index():  This function find the index of bc by its name
 * Warning, the two bcs for interface has the same name, but belogs to different zone
 */
int DABC::Get_bc_index(const char * bc_name)
{
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,bc_name))
    {
      return i;
    }
  return -1;
}

/* ----------------------------------------------------------------------------
 * DABC::Get_bc_index():  This function find the index of bc by its name (case insensitive)
 */
int DABC::Get_bc_index_nocase(const char * bc_name)
{
  for(int i=0;i<bc_point_array.size();i++)
  {
    if(!strcasecmp(bc_point_array[i]->psegment->label,bc_name))
    {
      return i;
    }
  }
  return -1;
}


// ----------------------------------------------------------------------------
// DABC:: electrode functions
// search electrode, attach source(s) to electrode, clear electrode source(s)...
// One can use segment label or electrode region label
// to specify the electrode
//
int DABC::is_electrode(int bc_index)
{
  return (bc_point_array[bc_index]->BCType & ContactMask);
}

int DABC::is_electrode(const char * bc_label)
{
  int flag=0;
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,bc_label) ||
        !strcmp(bc_point_array[i]->pzone->zonename,bc_label))
    {
      flag=(bc_point_array[i]->BCType & ContactMask);
      if (flag) break;
    }
  return flag;
}

BaseBC * DABC::Get_electrode_pointer(const char * bc_label)
{
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,bc_label) ||
        !strcmp(bc_point_array[i]->pzone->zonename,bc_label))
    {
      if(bc_point_array[i]->BCType & ContactMask) return bc_point_array[i];
    }
  return NULL;
}

BaseBC * DABC::Get_electrode_pointer_nocase(const char * bc_label)
{
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcasecmp(bc_point_array[i]->psegment->label,bc_label) ||
        !strcasecmp(bc_point_array[i]->pzone->zonename,bc_label))
    {
      if(bc_point_array[i]->BCType & ContactMask) return bc_point_array[i];
    }
  return NULL;
}

const char * DABC::format_electrode_name(const char * electrode)
{
  const char * formated_name = electrode;
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,electrode))
    {
      formated_name =electrode;
    }
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->pzone->zonename,electrode) && is_electrode(i))
    {
      formated_name = bc_point_array[i]->psegment->label;
    }
  return formated_name;
}

const char * DABC::format_electrode_name_nocase(const char * electrode)
{
  const char * formated_name = electrode;
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcasecmp(bc_point_array[i]->psegment->label,electrode))
    {
      formated_name =electrode;
    }
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcasecmp(bc_point_array[i]->pzone->zonename,electrode) && is_electrode(i))
    {
      formated_name = bc_point_array[i]->psegment->label;
    }
  return formated_name;
}

int  DABC::is_electrode_label(int bc_index, const char * bc_name)
{
  if(!strcmp(bc_point_array[bc_index]->psegment->label,format_electrode_name(bc_name)))
    return 1;
  else
    return 0;
}

int  DABC::is_electrode_label_nocase(int bc_index, const char * bc_name)
{
  if(!strcasecmp(bc_point_array[bc_index]->psegment->label,format_electrode_name_nocase(bc_name)))
    return 1;
  else
    return 0;
}

void DABC::Set_electrode_type(const char *bc_label,int type)
{
  const char * formated_name = format_electrode_name(bc_label);
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,formated_name))
      bc_point_array[i]->Set_electrode_type(type);
}

void DABC::Attach_Vapp(const char *bc_label, vector<VSource *> & vsrc)
{
  const char * formated_name = format_electrode_name(bc_label);
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,formated_name))
    {
      if(bc_point_array[i]->BCType==OhmicContact)
      {
        OhmicBC *pbc = dynamic_cast<OhmicBC *>(bc_point_array[i]);
        pbc->electrode_type = VoltageBC;
        pbc->vsrc = vsrc;
      }
      else if(bc_point_array[i]->BCType==SchottkyContact)
      {
        SchottkyBC *pbc = dynamic_cast<SchottkyBC *>(bc_point_array[i]);
        pbc->electrode_type = VoltageBC;
        pbc->vsrc = vsrc;
      }
      else if(bc_point_array[i]->BCType==GateContact)
      {
        GateBC *pbc = dynamic_cast<GateBC *>(bc_point_array[i]);
        pbc->vsrc = vsrc;
      }
      else if(bc_point_array[i]->BCType==InsulatorContact)
      {
        InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC *>(bc_point_array[i]);
        pbc->vsrc = vsrc;
      }
      else
      {
        sprintf(log_buf,"\nWarning: can't set voltage source(s) to BC %s, ignored.\n",bc_label);
        GSS_LOG();
      }
    }
  Update_Vapp(0);  
}

void DABC::Attach_Iapp(const char *bc_label, vector<ISource *> & isrc)
{
  const char * formated_name = format_electrode_name(bc_label);
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,formated_name))
    {
      if(bc_point_array[i]->BCType==OhmicContact)
      {
        OhmicBC *pbc = dynamic_cast<OhmicBC *>(bc_point_array[i]);
        pbc->electrode_type = CurrentBC;
        pbc->isrc = isrc;
      }
      else if(bc_point_array[i]->BCType==SchottkyContact)
      {
        SchottkyBC *pbc = dynamic_cast<SchottkyBC *>(bc_point_array[i]);
        pbc->electrode_type = CurrentBC;
        pbc->isrc = isrc;
      }
      else
      {
        sprintf(log_buf,"\nWarning: can't set current source(s) to BC %s, ignored.\n",bc_label);
        GSS_LOG();
      }
    }
  Update_Iapp(0);  
}

void DABC::Set_Vapp(const char *bc_label, PetscScalar V)
{
  const char * formated_name = format_electrode_name(bc_label);
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,formated_name))
    {
      bc_point_array[i]->Set_Vapp(V);
    }
}

void DABC::Set_Vapp_nocase(const char *bc_label, PetscScalar V)
{
  const char * formated_name = format_electrode_name_nocase(bc_label);
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcasecmp(bc_point_array[i]->psegment->label,formated_name) )
    {
      bc_point_array[i]->Set_Vapp(V);
    }
}

void DABC::Set_Vac(const char *bc_label, PetscScalar V)
{
  const char * formated_name = format_electrode_name(bc_label);
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,formated_name))
    {
      bc_point_array[i]->Set_Vac(V);
    }
}

void DABC::Set_Iapp(const char *bc_label, PetscScalar I)
{
  const char * formated_name = format_electrode_name(bc_label);
  for(int i=0;i<bc_point_array.size();i++)
    if(!strcmp(bc_point_array[i]->psegment->label,formated_name))
      bc_point_array[i]->Set_Iapp(I);
}

void DABC::Update_Vapp(PetscScalar t)
{
  for(int i=0;i<bc_point_array.size();i++)
  {
    if(bc_point_array[i]->BCType==OhmicContact)
    {
      OhmicBC *pbc = dynamic_cast<OhmicBC *>(Get_pointer(i));
      pbc->Vapp = 0;
      for(int j=0;j<pbc->vsrc.size();j++)
        pbc->Vapp+=pbc->vsrc[j]->vapp(t);
    }
    else if(bc_point_array[i]->BCType==SchottkyContact)
    {
      SchottkyBC *pbc = dynamic_cast<SchottkyBC *>(Get_pointer(i));
      pbc->Vapp = 0;
      for(int j=0;j<pbc->vsrc.size();j++)
        pbc->Vapp+=pbc->vsrc[j]->vapp(t);
    }
    else if(bc_point_array[i]->BCType==GateContact)
    {
      GateBC *pbc = dynamic_cast<GateBC *>(Get_pointer(i));
      pbc->Vapp = 0;
      for(int j=0;j<pbc->vsrc.size();j++)
        pbc->Vapp+=pbc->vsrc[j]->vapp(t);
    }
    else if(bc_point_array[i]->BCType==InsulatorContact)
    {
      InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC *>(Get_pointer(i));
      pbc->Vapp = 0;
      for(int j=0;j<pbc->vsrc.size();j++)
        pbc->Vapp+=pbc->vsrc[j]->vapp(t);
    }
  }

}

void DABC::Clear_Vapp()
{
  for(int i=0;i<bc_point_array.size();i++)
  {
    if(bc_point_array[i]->BCType==OhmicContact)
    {
      OhmicBC *pbc = dynamic_cast<OhmicBC *>(Get_pointer(i));
      pbc->Vapp = 0;
    }
    else if(bc_point_array[i]->BCType==SchottkyContact)
    {
      SchottkyBC *pbc = dynamic_cast<SchottkyBC *>(Get_pointer(i));
      pbc->Vapp = 0;
    }
    else if(bc_point_array[i]->BCType==GateContact)
    {
      GateBC *pbc = dynamic_cast<GateBC *>(Get_pointer(i));
      pbc->Vapp = 0;
    }
    else if(bc_point_array[i]->BCType==InsulatorContact)
    {
      InsulatorContactBC *pbc = dynamic_cast<InsulatorContactBC *>(Get_pointer(i));
      pbc->Vapp = 0;
    }
  }
}

void DABC::Update_Iapp(PetscScalar t)
{
  for(int i=0;i<bc_point_array.size();i++)
  {
    if(bc_point_array[i]->BCType==OhmicContact)
    {
      OhmicBC *pbc = dynamic_cast<OhmicBC *>(Get_pointer(i));
      pbc->Iapp = 0;
      for(int j=0;j<pbc->isrc.size();j++)
        pbc->Iapp+=pbc->isrc[j]->iapp(t);
    }
    else if(bc_point_array[i]->BCType==SchottkyContact)
    {
      SchottkyBC *pbc = dynamic_cast<SchottkyBC *>(Get_pointer(i));
      pbc->Iapp = 0;
      for(int j=0;j<pbc->isrc.size();j++)
        pbc->Iapp+=pbc->isrc[j]->iapp(t);
    }
  }
}

void DABC::Clear_Iapp()
{
  for(int i=0;i<bc_point_array.size();i++)
  {
    if(bc_point_array[i]->BCType==OhmicContact)
    {
      OhmicBC *pbc = dynamic_cast<OhmicBC *>(Get_pointer(i));
      pbc->Iapp = 0;
    }
    else if(bc_point_array[i]->BCType==SchottkyContact)
    {
      SchottkyBC *pbc = dynamic_cast<SchottkyBC *>(Get_pointer(i));
      pbc->Iapp = 0;
    }
  }
}



int InterfaceType(char *mat_name1,char *mat_name2)
{
  int mat1 = IsSemiconductor(mat_name1);
  int mat2 = IsSemiconductor(mat_name2);
  if(mat1>0 && mat2>0 && mat1 != mat2) return HeteroInterface;
  if(mat1>0 && mat1 == mat2) return HomoInterface;
  if(mat1*mat2==0 && mat1+mat2>0) //one material is semiconductor
  {
    if(IsInsulator(mat_name1)+IsInsulator(mat_name2)>0) return IF_Insulator_Semiconductor;
    if(IsVacuum(mat_name1)+IsVacuum(mat_name2)>0) return IF_Semiconductor_Vacuum;
    if(IsElectrode(mat_name1)+IsElectrode(mat_name2)>0) return IF_Electrode_Semiconductor;
  }

  mat1 = IsInsulator(mat_name1);
  mat2 = IsInsulator(mat_name2);
  if(mat1>0 && mat2>0) return IF_Insulator_Insulator;
  if(mat1*mat2==0 && mat1+mat2>0) //one material is Insulator
  {
    if(IsVacuum(mat_name1)+IsVacuum(mat_name2)>0) return IF_Insulator_Vacuum;
    if(IsElectrode(mat_name1)+IsElectrode(mat_name2)>0) return IF_Electrode_Insulator;
  }

  mat1 = IsElectrode(mat_name1);
  mat2 = IsElectrode(mat_name2);
  if(mat1>0 && mat2>0) return IF_Electrode_Electrode;
  if(mat1*mat2==0 && mat1+mat2>0) //one material is Electrode
  {
    if(IsVacuum(mat_name1)+IsVacuum(mat_name2)>0) return IF_Electrode_Vacuum;
  }
  return 0;
}

/* ----------------------------------------------------------------------------
 * DABC::InitBC():  This function set up each bc by segment information from mesh
 *  and boundary information from command file.
 */
int DABC::InitBC(int nbc, Segment** asegment,int nzone,ZONE* zone, ZoneInterface &azintface, CmdBuf *pcmdbuf,
                 double Tl,PhysicalUnitScale *unit_scale)
{

  bc_num = nbc;
  bc_point_array.resize(bc_num);
  scale = unit_scale;
  lattice_temperature = Tl;
  psegment_array = asegment;
  pzintface = &azintface;

  for(int i=0; i<bc_num; i++)
    bc_point_array[i] = NULL;
  sprintf(log_buf,"\n--------------------------------------------------\n");GSS_LOG();
  sprintf(log_buf,"\ninit each boundary condition...\n");GSS_LOG();
  //set bc parameters by command file
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
  {
    if(pcmdbuf->is_current_cmd("BOUNDARY"))   // It's a boundary card
    {
      list<Cmd>::iterator    pcmd = pcmdbuf->get_current_cmd();
      if(pcmd->is_arg_exist("type"))
      {
        if(pcmd->is_arg_value("type","NeumannBoundary"))
        {if(SetBCNeumannBoundary(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","OhmicContact"))
        {if(SetBCOhmicContact(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","InsulatorInterface"))
        {if(SetBCInsulatorInterface(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","Heterojunction"))
        {if(SetBCHeterojunction(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","SchottkyContact"))
        {if(SetBCSchottkyContact(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","GateContact"))
        {if(SetBCGateContact(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","InsulatorContact"))
        {if(SetBCInsulatorContact(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","ChargedContact"))
        {if(SetBCChargedContact(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","AbsorbingBoundary"))  //```
        {if(SetBCAbsorbingBoundary(pcmd)) return 1;}
        else if(pcmd->is_arg_value("type","SourceBoundary"))
        {if(SetBCSourceBoundary(pcmd)) return 1;}
        else
        {
          sprintf(log_buf,"\nline %d BOUNDARY: no such BC Type!\n",pcmd->lineno);GSS_LOG();
          return 1;
        }
      }
    }

    if(pcmdbuf->is_current_cmd("CONTACT"))   // It's a contact card
    {
      list<Cmd>::iterator    pcmd = pcmdbuf->get_current_cmd();
      if(pcmd->is_arg_exist("type"))
      {
        if(pcmd->is_arg_value("type","OhmicContact"))
        {if(SetElectrodeOhmicContact(pcmd,nzone,zone)) return 1;}
        if(pcmd->is_arg_value("type","GateContact"))
        {if(SetElectrodeGateContact(pcmd,nzone,zone)) return 1;}
        if(pcmd->is_arg_value("type","SchottkyContact"))
        {if(SetElectrodeSchottkyContact(pcmd,nzone,zone)) return 1;}
        if(pcmd->is_arg_value("type","FloatMetal"))
        {if(SetFloatMetal(pcmd,nzone,zone)) return 1;}
      }
    }
  }

  //set Electrode, Insulator, Homo and Hetero Interface
  for(int i=0; i<bc_num; i++)
    if(bc_point_array[i] == NULL && psegment_array[i]->interface!=-1)
    {
      Interface *pinterface = &pzintface->interface[psegment_array[i]->interface];
      //int j;
      //for(j=0; j<bc_num; j++)
      //  if(i!=j && !strcmp(psegment_array[j]->label,psegment_array[i]->label)) 	break;
      int interface=InterfaceType(pinterface->pzone1->zonelabel,pinterface->pzone2->zonelabel);

      if(interface==IF_Electrode_Semiconductor)
      {
        bc_point_array[i] = new OhmicBC;
        bc_point_array[i]->BCType = OhmicContact;
        bc_point_array[i]->psegment = psegment_array[i];
        OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc_point_array[i]);
        pbc->T_external = lattice_temperature;
        pbc->heat_transfer = 1e3*scale->s_joule/scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
        pbc->R = 0*scale->s_volt/scale->s_A;
        pbc->C = 0*scale->s_coulomb/scale->s_volt;
        pbc->L = 0*scale->s_volt*scale->s_second/scale->s_A;
	pbc->pinterface = pinterface;
        sprintf(log_buf,"\nBC[%d] label \"%s\" init as OhmicContact...done",i+1,psegment_array[i]->label);
        GSS_LOG();
      }
      if(interface==IF_Insulator_Semiconductor)
      {
        bc_point_array[i] = new InsulatorInterfaceBC;
        bc_point_array[i]->BCType = InsulatorInterface;
        bc_point_array[i]->psegment = psegment_array[i];
        InsulatorInterfaceBC * pbc = dynamic_cast <InsulatorInterfaceBC * > (bc_point_array[i]);
        pbc->pinterface = pinterface;
        pbc->QF = 1e10/pow(scale->s_centimeter,2);
        sprintf(log_buf,"\nBC[%d] label \"%s\" init as insulator interface...done",i+1,psegment_array[i]->label);
        GSS_LOG();
      }
      //Some Ohmic/Schottky Contact electrode have IF to Insulator. Set as ElectrodeInsulatorBC.
      if(interface==IF_Electrode_Insulator)
      {
        bc_point_array[i] = new ElectrodeInsulatorBC;
        bc_point_array[i]->BCType = IF_Electrode_Insulator;
        bc_point_array[i]->psegment = psegment_array[i];
        ElectrodeInsulatorBC *pbc = dynamic_cast <ElectrodeInsulatorBC * > (bc_point_array[i]);
	pbc->pinterface = pinterface;
	sprintf(log_buf,"\nBC[%d] label \"%s\" init as ElectrodeInsulatorBC...done",i+1,psegment_array[i]->label);
        GSS_LOG();
      }
      if(interface==HeteroInterface)
      {
        bc_point_array[i] = new HeteroInterfaceBC;
        bc_point_array[i]->BCType = HeteroInterface;
        bc_point_array[i]->psegment = psegment_array[i];
        HeteroInterfaceBC *pbc = dynamic_cast <HeteroInterfaceBC * > (bc_point_array[i]);
        pbc->pinterface = pinterface;
        pbc->QF = 0;
        sprintf(log_buf,"\nBC[%d] label \"%s\" init as heterojuction...done",i+1,psegment_array[i]->label);
        GSS_LOG();
      }
      if(interface==HomoInterface)
      {
        bc_point_array[i] = new HomoInterfaceBC;
        bc_point_array[i]->BCType = HomoInterface;
        bc_point_array[i]->psegment = psegment_array[i];
        HomoInterfaceBC *pbc = dynamic_cast <HomoInterfaceBC * > (bc_point_array[i]);
        pbc->pinterface = pinterface;
        sprintf(log_buf,"\nBC[%d] label \"%s\" init as homojuction...done",i+1,psegment_array[i]->label);
        GSS_LOG();
      }
      if(interface==IF_Insulator_Insulator)
      {
        bc_point_array[i] = new InsulatorInsulatorBC;
        bc_point_array[i]->BCType = IF_Insulator_Insulator;
        bc_point_array[i]->psegment = psegment_array[i];
        InsulatorInsulatorBC *pbc = dynamic_cast <InsulatorInsulatorBC * > (bc_point_array[i]);
        pbc->pinterface = pinterface;
        sprintf(log_buf,"\nBC[%d] label \"%s\" init as InsulatorInsulatorBC...done",i+1,psegment_array[i]->label);
        GSS_LOG();
      }
      if(interface==IF_Electrode_Electrode)
      {
        bc_point_array[i] = new ElectrodeElectrodeBC;
        bc_point_array[i]->BCType = IF_Electrode_Electrode;
        bc_point_array[i]->psegment = psegment_array[i];
        ElectrodeElectrodeBC *pbc = dynamic_cast <ElectrodeElectrodeBC * > (bc_point_array[i]);
        pbc->pinterface = pinterface;
        sprintf(log_buf,"\nBC[%d] label \"%s\" init as ElectrodeElectrodeBC...done",i+1,psegment_array[i]->label);
        GSS_LOG();
      }
    }
    
  //set remain bc as NeumannBoundary
  for(int i=0; i<bc_num; i++)
    if(bc_point_array[i] == NULL)
    {
      bc_point_array[i] = new NeumannBC;
      bc_point_array[i]->BCType = NeumannBoundary;
      bc_point_array[i]->psegment = psegment_array[i];
      NeumannBC *pbc = dynamic_cast <NeumannBC * > (bc_point_array[i]);
      pbc->T_external = lattice_temperature;
      if(IsElectrode(zone[psegment_array[i]->zone_index].zonelabel))
        pbc->heat_transfer = 1e3*scale->s_joule/scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
      else
        pbc->heat_transfer = 0;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as NeumannBoundary...done",i+1,psegment_array[i]->label);
      GSS_LOG();
    }

  //set the pzone pointer
  for(int i=0; i<bc_num; i++)
    bc_point_array[i]->pzone = & zone[bc_point_array[i]->psegment->zone_index];

  //ok
  sprintf(log_buf,"\nboundary condition init ok\n\n");GSS_LOG();

  return 0;
}

// ----------------------------------------------------------------------------
int DABC::SetBCNeumannBoundary(list<Cmd>::iterator pcmd)
{
  int flag = 0;
  list<Arg>::iterator  parg;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(4,"type","id","ext.temp","heat.transfer"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as NeumannBoundary...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new NeumannBC;
      bc_point_array[i]->BCType = NeumannBoundary;
      bc_point_array[i]->psegment = psegment_array[i];
      NeumannBC *pbc = dynamic_cast <NeumannBC * > (bc_point_array[i]);
      pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
      pbc->heat_transfer = pcmd->get_number("heat.transfer",0,0.0)*scale->s_joule
                           /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}

int DABC::SetBCOhmicContact(list<Cmd>::iterator pcmd)
{
  int flag = 0;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary\n!",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(8,"type","id","ext.temp","heat.transfer","res","cap","ind","connectto"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as OhmicContact...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new OhmicBC;
      bc_point_array[i]->BCType = OhmicContact;
      bc_point_array[i]->psegment = psegment_array[i];
      OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc_point_array[i]);
      pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
      pbc->heat_transfer = pcmd->get_number("heat.transfer",0,1e3)*scale->s_joule
                           /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
      pbc->R = pcmd->get_number("res",0,0.0)*scale->s_volt/scale->s_A;
      pbc->C = pcmd->get_number("cap",0,0.0)*scale->s_coulomb/scale->s_volt;
      pbc->L = pcmd->get_number("ind",0,0.0)*scale->s_volt*scale->s_second/scale->s_A;
      pbc->connect_elec=pcmd->get_string("connectto",0,"");
      if(psegment_array[i]->interface!=-1)
        pbc->pinterface= & pzintface->interface[psegment_array[i]->interface];
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


int DABC::SetBCSchottkyContact(list<Cmd>::iterator pcmd)
{
  int flag = 0;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(9,"type","id","ext.temp","heat.transfer","res","cap","ind","workfunction","connectto"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as Schottky...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new SchottkyBC;
      bc_point_array[i]->BCType = SchottkyContact;
      bc_point_array[i]->psegment = psegment_array[i];
      SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc_point_array[i]);
      pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
      pbc->heat_transfer = pcmd->get_number("heat.transfer",0,1e3)*scale->s_joule
                           /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
      pbc->R = pcmd->get_number("res",0,0.0)*scale->s_volt/scale->s_A;
      pbc->C = pcmd->get_number("cap",0,0.0)*scale->s_coulomb/scale->s_volt;
      pbc->L = pcmd->get_number("ind",0,0.0)*scale->s_volt*scale->s_second/scale->s_A;
      pbc->WorkFunction = pcmd->get_number("workfunction",0,4.7)*scale->s_volt;
      pbc->connect_elec=pcmd->get_string("connectto",0,"");
      if(psegment_array[i]->interface!=-1)
        pbc->pinterface= & pzintface->interface[psegment_array[i]->interface];
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


int DABC::SetBCGateContact(list<Cmd>::iterator pcmd)
{
  int flag = 0;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary\n!",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(8,"type","id","ext.temp","heat.transfer","res","cap","ind","workfunction"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as GateContact...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new GateBC;
      bc_point_array[i]->BCType = GateContact;
      bc_point_array[i]->psegment = psegment_array[i];
      GateBC *pbc = dynamic_cast <GateBC * > (bc_point_array[i]);
      pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
      pbc->heat_transfer = pcmd->get_number("heat.transfer",0,1e3)*scale->s_joule
                           /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
      pbc->R = pcmd->get_number("res",0,0.0)*scale->s_volt/scale->s_A;
      pbc->C = pcmd->get_number("cap",0,0.0)*scale->s_coulomb/scale->s_volt;
      pbc->L = pcmd->get_number("ind",0,0.0)*scale->s_volt*scale->s_second/scale->s_A;
      pbc->WorkFunction=pcmd->get_number("workfunction",0,4.17)*scale->s_volt;
      if(psegment_array[i]->interface!=-1)
        pbc->pinterface= & pzintface->interface[psegment_array[i]->interface];
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


int DABC::SetBCInsulatorInterface(list<Cmd>::iterator pcmd)
{
  int flag = 0;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(3,"type","id","qf"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as InsulatorInterface...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new InsulatorInterfaceBC;
      bc_point_array[i]->BCType = InsulatorInterface;
      bc_point_array[i]->psegment = psegment_array[i];
      InsulatorInterfaceBC * pbc = dynamic_cast <InsulatorInterfaceBC * > (bc_point_array[i]);
      pbc->pinterface = &pzintface->interface[psegment_array[i]->interface];
      pbc->QF=pcmd->get_number("qf",0,1e10)/pow(scale->s_centimeter,2);
    }
  }
  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}

int DABC::SetBCHeterojunction(list<Cmd>::iterator pcmd)
{
  int flag = 0;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(3,"type","id","qf"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as heterojunction...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new HeteroInterfaceBC;
      bc_point_array[i]->BCType = HeteroInterface;
      bc_point_array[i]->psegment = psegment_array[i];
      HeteroInterfaceBC * pbc = dynamic_cast <HeteroInterfaceBC * > (bc_point_array[i]);
      pbc->pinterface = & pzintface->interface[psegment_array[i]->interface];
      pbc->QF=pcmd->get_number("qf",0,0.0)/pow(scale->s_centimeter,2);
    }
  }
  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


int DABC::SetBCInsulatorContact(list<Cmd>::iterator pcmd)
{
  int flag = 0;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(11,"type","id","ext.temp","heat.transfer","thickness","eps","res","cap","ind","qf","workfunction"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as InsulatorContact...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new InsulatorContactBC;
      bc_point_array[i]->BCType = InsulatorContact;
      bc_point_array[i]->psegment = psegment_array[i];
      InsulatorContactBC *pbc = dynamic_cast <InsulatorContactBC * > (bc_point_array[i]);
      pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
      pbc->heat_transfer = pcmd->get_number("heat.transfer",0,1e3)*scale->s_joule
                           /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
      pbc->Thick = pcmd->get_number("thickness",0,2e-7)*scale->s_centimeter;
      pbc->eps = pcmd->get_number("eps",0,3.9);
      pbc->R = pcmd->get_number("res",0,0.0)*scale->s_volt/scale->s_A;
      pbc->C = pcmd->get_number("cap",0,0.0)*scale->s_coulomb/scale->s_volt;
      pbc->L = pcmd->get_number("ind",0,0.0)*scale->s_volt*scale->s_second/scale->s_A;
      pbc->WorkFunction = pcmd->get_number("workfunction",0,4.7)*scale->s_volt;
      pbc->QF = pcmd->get_number("qf",0,1e10)/pow(scale->s_centimeter,2);
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


int DABC::SetBCChargedContact(list<Cmd>::iterator pcmd)
{
  int flag = 0;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(3,"type","id","qf"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as ChargedContact...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new ChargedContactBC;
      bc_point_array[i]->BCType = ChargedContact;
      bc_point_array[i]->psegment = psegment_array[i];
      ChargedContactBC *pbc = dynamic_cast <ChargedContactBC * > (bc_point_array[i]);
      pbc->pinterface = & pzintface->interface[psegment_array[i]->interface];
      pbc->QF = pcmd->get_number("qf",0,0.0)*scale->s_coulomb/scale->s_micron;
      pbc->pinterface= & pzintface->interface[psegment_array[i]->interface];
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


int DABC::SetBCAbsorbingBoundary(list<Cmd>::iterator pcmd)//by zhangxih ----06-09-27
{
  int flag = 0;
  list<Arg>::iterator  parg;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(4,"type","id","ext.temp","heat.transfer"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as AbsorbingBoundary...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new AbsorbingBC;
      bc_point_array[i]->BCType = AbsorbingBoundary;
      bc_point_array[i]->psegment = psegment_array[i];
      AbsorbingBC *pbc = dynamic_cast <AbsorbingBC * > (bc_point_array[i]);
      pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
      pbc->heat_transfer = pcmd->get_number("heat.transfer",0,0.0)*scale->s_joule
                           /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


int DABC::SetBCSourceBoundary(list<Cmd>::iterator pcmd)//by zhangxih ----06-09-27
{
  int flag = 0;
  list<Arg>::iterator  parg;

  char *Identifier;
  if(pcmd->is_arg_exist("id"))
    Identifier=pcmd->get_string("id",0,"");
  else
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: You must give ID for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(4,"type","id","ext.temp","heat.transfer"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  for(int i=0;i<bc_num;i++)
  {
    if(!strcmp(psegment_array[i]->label,Identifier))
    {
      flag = 1;
      sprintf(log_buf,"\nBC[%d] label \"%s\" init as SourceBoundary...done",i+1,Identifier);
      GSS_LOG();
      bc_point_array[i] = new SourceBC;
      bc_point_array[i]->BCType = SourceBoundary;
      bc_point_array[i]->psegment = psegment_array[i];
      SourceBC *pbc = dynamic_cast <SourceBC * > (bc_point_array[i]);
      pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
      pbc->heat_transfer = pcmd->get_number("heat.transfer",0,0.0)*scale->s_joule
                           /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
    }
  }

  if(!flag)
  {
    sprintf(log_buf,"\nerror: BC label \"%s\" not matched\n",Identifier);
    GSS_LOG();
    return 1;
  }
  return 0;
}


//------------------------------------------------------------------------------------

int DABC::SetElectrodeOhmicContact(list<Cmd>::iterator pcmd,int nzone,ZONE* zone)
{
  ZONE* pzone=NULL;
  char *Identifier;
  Interface *pinterface;

  Identifier=pcmd->get_string("id",0,"");
  for(int i=0;i<nzone;i++)
    if(!strcmp(zone[i].zonename,Identifier) && IsElectrode(zone[i].zonelabel))
      pzone=&zone[i];
  if(!pzone)
  {
    sprintf(log_buf,"\nLine %d CONTACT: You must give correct region name for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(8,"type","id","ext.temp","heat.transfer","res","cap","ind","connectto"))
  {
    sprintf(log_buf,"\nLine %d CONTACT: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  for(int i=0;i<bc_num;i++)
    for(int j=0;j<pzone->dasegment.size();j++)
    {
      if(!strcmp(psegment_array[i]->label,pzone->dasegment[j].label))
      {
        if(pzone->dasegment[j].interface!=-1) // this segment is an interface
        {
          pinterface = & pzintface->interface[pzone->dasegment[j].interface];
          int interface=InterfaceType(pinterface->pzone1->zonelabel,pinterface->pzone2->zonelabel);
          if(interface==IF_Electrode_Semiconductor) //the segment is an interface to semiconductor, set it as ohmic contact
          {
            sprintf(log_buf,"\nBC[%d] label \"%s\" init as OhmicContact...done",i+1,pzone->dasegment[j].label);
            GSS_LOG();
            bc_point_array[i] = new OhmicBC;
            bc_point_array[i]->BCType = OhmicContact;
            bc_point_array[i]->psegment = psegment_array[i];
            OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc_point_array[i]);
            pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
            pbc->heat_transfer = pcmd->get_number("heat.transfer",0,1e3)*scale->s_joule
                                 /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
            pbc->R = pcmd->get_number("res",0,0.0)*scale->s_volt/scale->s_A;
            pbc->C = pcmd->get_number("cap",0,0.0)*scale->s_coulomb/scale->s_volt;
            pbc->L = pcmd->get_number("ind",0,0.0)*scale->s_volt*scale->s_second/scale->s_A;
            pbc->connect_elec=pcmd->get_string("connectto",0,"");
            pbc->pinterface = pinterface;
          }
          if(interface==IF_Electrode_Insulator)
          {
            sprintf(log_buf,"\nBC[%d] label \"%s\" init as ElectrodeInsulatorBC...done",i+1,pzone->dasegment[j].label);
            GSS_LOG();
            bc_point_array[i] = new ElectrodeInsulatorBC;
            bc_point_array[i]->BCType = IF_Electrode_Insulator;
            bc_point_array[i]->psegment = psegment_array[i];
	    ElectrodeInsulatorBC *pbc = dynamic_cast <ElectrodeInsulatorBC * > (bc_point_array[i]);
	    pbc->pinterface = pinterface;
          }
        }
        // the remain segment will be set as Newmann later, no process is needed here.
      }
    }
  return 0;
}

int DABC::SetElectrodeSchottkyContact(list<Cmd>::iterator pcmd,int nzone,ZONE* zone)
{
  ZONE* pzone=NULL;
  char *Identifier;
  Interface *pinterface;

  Identifier=pcmd->get_string("id",0,"");
  for(int i=0;i<nzone;i++)
    if(!strcmp(zone[i].zonename,Identifier) && IsElectrode(zone[i].zonelabel))
      pzone=&zone[i];
  if(!pzone)
  {
    sprintf(log_buf,"\nLine %d CONTACT: You must give correct region name for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(9,"type","id","ext.temp","heat.transfer","res","cap","ind","workfunction","connectto"))
  {
    sprintf(log_buf,"\nLine %d CONTACT: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  for(int i=0;i<bc_num;i++)
    for(int j=0;j<pzone->dasegment.size();j++)
    {
      if(!strcmp(psegment_array[i]->label,pzone->dasegment[j].label))
      {
        if(pzone->dasegment[j].interface!=-1) // this segment is an interface
        {
          pinterface = & pzintface->interface[pzone->dasegment[j].interface];
          int interface=InterfaceType(pinterface->pzone1->zonelabel,pinterface->pzone2->zonelabel);
          if(interface==IF_Electrode_Semiconductor) //the segment is an interface to semiconductor, set it as ohmic contact
          {
            sprintf(log_buf,"\nBC[%d] label \"%s\" init as SchottkyContact...done",i+1,pzone->dasegment[j].label);
            GSS_LOG();
            bc_point_array[i] = new SchottkyBC;
            bc_point_array[i]->BCType = SchottkyContact;
            bc_point_array[i]->psegment = psegment_array[i];
            SchottkyBC *pbc = dynamic_cast <SchottkyBC * > (bc_point_array[i]);
            pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
            pbc->heat_transfer = pcmd->get_number("heat.transfer",0,1e3)*scale->s_joule
                                 /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
            pbc->R = pcmd->get_number("res",0,0.0)*scale->s_volt/scale->s_A;
            pbc->C = pcmd->get_number("cap",0,0.0)*scale->s_coulomb/scale->s_volt;
            pbc->L = pcmd->get_number("ind",0,0.0)*scale->s_volt*scale->s_second/scale->s_A;
            pbc->WorkFunction = pcmd->get_number("workfunction",0,4.7)*scale->s_volt;
            pbc->connect_elec=pcmd->get_string("connectto",0,"");
            pbc->pinterface = pinterface;
          }
          if(interface==IF_Electrode_Insulator)
          {
            sprintf(log_buf,"\nBC[%d] label \"%s\" init as ElectrodeInsulatorBC...done",i+1,pzone->dasegment[j].label);
            GSS_LOG();
            bc_point_array[i] = new ElectrodeInsulatorBC;
            bc_point_array[i]->BCType = IF_Electrode_Insulator;
            bc_point_array[i]->psegment = psegment_array[i];
	    ElectrodeInsulatorBC *pbc = dynamic_cast <ElectrodeInsulatorBC * > (bc_point_array[i]);
	    pbc->pinterface = pinterface;
          }
        }
        // the remain segment will be set as Newmann later, no process is needed here.
      }
    }
  return 0;
}


int DABC::SetElectrodeGateContact(list<Cmd>::iterator pcmd,int nzone,ZONE* zone)
{
  ZONE* pzone=NULL;
  char *Identifier;
  Interface *pinterface;

  Identifier=pcmd->get_string("id",0,"");
  for(int i=0;i<nzone;i++)
    if(!strcmp(zone[i].zonename,Identifier) && IsElectrode(zone[i].zonelabel))
      pzone=&zone[i];
  if(!pzone)
  {
    sprintf(log_buf,"\nLine %d CONTACT: You must give correct region name for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(9,"type","id","ext.temp","heat.transfer","res","cap","ind","workfunction","connectto"))
  {
    sprintf(log_buf,"\nLine %d CONTACT: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  for(int i=0;i<bc_num;i++)
    for(int j=0;j<pzone->dasegment.size();j++)
    {
      if(!strcmp(psegment_array[i]->label,pzone->dasegment[j].label))
      {
        if(pzone->dasegment[j].interface!=-1) // this segment is an interface
        {
          pinterface = & pzintface->interface[pzone->dasegment[j].interface];
          int interface=InterfaceType(pinterface->pzone1->zonelabel,pinterface->pzone2->zonelabel);

          if(interface==IF_Electrode_Insulator)
          {
            sprintf(log_buf,"\nBC[%d] label \"%s\" init as GateContact...done",i+1,pzone->dasegment[j].label);
            GSS_LOG();
            bc_point_array[i] = new GateBC;
            bc_point_array[i]->BCType = GateContact;
            bc_point_array[i]->psegment = psegment_array[i];
            GateBC *pbc = dynamic_cast <GateBC * > (bc_point_array[i]);
            pbc->T_external = pcmd->get_number("ext.temp",0,lattice_temperature/scale->s_kelvin)*scale->s_kelvin;
            pbc->heat_transfer = pcmd->get_number("heat.transfer",0,1e3)*scale->s_joule
                                 /scale->s_second/pow(scale->s_centimeter,2)/scale->s_kelvin;
            pbc->R = pcmd->get_number("res",0,0.0)*scale->s_volt/scale->s_A;
            pbc->C = pcmd->get_number("cap",0,0.0)*scale->s_coulomb/scale->s_volt;
            pbc->L = pcmd->get_number("ind",0,0.0)*scale->s_volt*scale->s_second/scale->s_A;
            pbc->WorkFunction=pcmd->get_number("workfunction",0,4.17)*scale->s_volt;
	    pbc->pinterface = pinterface;
          }
        }
        // the remain segment will be set as Newmann later, no process is needed here.
      }
    }
  return 0;
}


int DABC::SetFloatMetal(list<Cmd>::iterator pcmd,int nzone,ZONE* zone)
{
  ZONE* pzone=NULL;
  char *Identifier;
  Interface *pinterface;

  Identifier=pcmd->get_string("id",0,"");
  for(int i=0;i<nzone;i++)
    if(!strcmp(zone[i].zonename,Identifier) && IsElectrode(zone[i].zonelabel))
      pzone=&zone[i];
  if(!pzone)
  {
    sprintf(log_buf,"\nLine %d CONTACT: You must give correct region name for this boundary!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(!pcmd->allowed_args(3,"type","id","qf"))
  {
    sprintf(log_buf,"\nLine %d BOUNDARY: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  for(int i=0;i<bc_num;i++)
    for(int j=0;j<pzone->dasegment.size();j++)
    {
      if(!strcmp(psegment_array[i]->label,pzone->dasegment[j].label))
      {
        if(pzone->dasegment[j].interface!=-1) // this segment is an interface
        {
          pinterface = & pzintface->interface[pzone->dasegment[j].interface];
          int interface=InterfaceType(pinterface->pzone1->zonelabel,pinterface->pzone2->zonelabel);

          if(interface==IF_Electrode_Insulator)
          {
            sprintf(log_buf,"\nBC[%d] label \"%s\" init as ChargedContact...done",i+1,Identifier);
            GSS_LOG();
            bc_point_array[i] = new ChargedContactBC;
            bc_point_array[i]->BCType = ChargedContact;
            bc_point_array[i]->psegment = psegment_array[i];
            ChargedContactBC *pbc = dynamic_cast <ChargedContactBC * > (bc_point_array[i]);
            pbc->pinterface = & pzintface->interface[psegment_array[i]->interface];
            pbc->QF = pcmd->get_number("qf",0,0.0)*scale->s_coulomb/scale->s_micron;
          }
        }
        // the remain segment will be set as Newmann later, no process is needed here.
      }
    }
  return 0;
}


/* ----------------------------------------------------------------------------
 * DABC::EmitTo  This function return the likely bc which tunneling electron may inject to
 */
int  DABC::EmitTo(int zone_index, PetscScalar source_x, PetscScalar source_y, PetscScalar ex, PetscScalar ey)
{
  for( int i=0;i<size(); i++)
    if(bc_point_array[i]->pzone->zone_index==zone_index &&
        (bc_point_array[i]->BCType==InsulatorInterface ||
         bc_point_array[i]->BCType==IF_Insulator_Semiconductor ||
         bc_point_array[i]->BCType==GateContact ))
    {
      for (int j=0;j<bc_point_array[i]->psegment->edge_array.size();j++)
      {
        Edge e = bc_point_array[i]->psegment->edge_array[j];
        PetscScalar d1x=bc_point_array[i]->pzone->danode[e.p1].x;
        PetscScalar d1y=bc_point_array[i]->pzone->danode[e.p1].y;
        PetscScalar d2x=bc_point_array[i]->pzone->danode[e.p2].x;
        PetscScalar d2y=bc_point_array[i]->pzone->danode[e.p2].y;
        if(RaySegmentIntersectTest(source_x,source_y,ex,ey,d1x,d1y,d2x,d2y)==INTERESECTING)
          return i;
      }
    }
  return -1;
}

