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
/*  Last update: March 13, 2006                                              */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "bsolver.h"
#include "log.h"
#include "typedef.h"
#include "adolc.h"
unsigned int adtl::AutoDScalar::numdir = 9;

int BSolver::build_zonedata()
{
  zonedata.resize(zone_num);
  // determine zone material by zone label. the label is compatible with medici
  for(int i=0; i<zone_num; i++)
  { //it is an Insulator zone
    if(IsInsulator(zone[i].zonelabel))
    {
      zonedata[i] = new ISZone;
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[i]);
      pzonedata->pzone = &zone[i];
      pzonedata->zone_index = zone[i].zone_index;
      pzonedata->pzone = &zone[i];
      pzonedata->node_num = zone[i].davcell.size();
      pzonedata->tri_num  = zone[i].datri.size();
      pzonedata->material_type = Insulator;
      pzonedata->material = IsInsulator(zone[i].zonelabel);
      try {pzonedata->mt = new MatInsulator(zone[i].zonename,
                          FormatInsulatorString(zone[i].zonelabel),&scale_unit);}
      catch(int){return 1;}
      pzonedata->fs = new ISData[pzonedata->node_num];
      pzonedata->aux = new ISAuxData[pzonedata->node_num];
      pzonedata->electrode.clear();
      for(int j=0;j<pzonedata->pzone->dasegment.size();j++)
      {
        int bc_index = pzonedata->pzone->dasegment[j].bc_index-1;
        if(bc[bc_index].BCType==GateContact ||
           bc[bc_index].BCType==ChargedContact )
          pzonedata->electrode.push_back(bc_index);
      }
      //pzonedata->report();
    }
    //if this is a semiconductor zone
    else if(IsSemiconductor(zone[i].zonelabel))
    {
      zonedata[i] = new SMCZone;
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[i]);
      pzonedata->pzone = &zone[i];
      pzonedata->zone_index = zone[i].zone_index;
      pzonedata->material = IsSemiconductor(zone[i].zonelabel);
      pzonedata->material_type = Semiconductor;
      try{pzonedata->mt = new MatSemiconductor(zone[i].zonename,
                          FormatSemiconductorString(zone[i].zonelabel),PMIS_define,&scale_unit);}
      catch(int){return 1;}
      pzonedata->node_num = zone[i].davcell.size();
      pzonedata->tri_num  = zone[i].datri.size();
      pzonedata->fs  = new SemiData[pzonedata->node_num];
      pzonedata->aux = new SemiAuxData[pzonedata->node_num];
      pzonedata->electrode.clear();
      for(int j=0;j<pzonedata->pzone->dasegment.size();j++)
      {
        int bc_index = pzonedata->pzone->dasegment[j].bc_index-1;
        if(bc[bc_index].BCType==OhmicContact||
            bc[bc_index].BCType==SchottkyContact||
            bc[bc_index].BCType==InsulatorContact)
          pzonedata->electrode.push_back(bc_index);
      }
      //pzonedata->report();
    }
    //oh it is an electrode,medici can generate this useless region
    else if(IsElectrode(zone[i].zonelabel))
    {
      zonedata[i] = new ElZone;
      ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[i]);
      pzonedata->pzone = &zone[i];
      pzonedata->zone_index = zone[i].zone_index;
      pzonedata->node_num = zone[i].davcell.size();
      pzonedata->tri_num  = zone[i].datri.size();
      pzonedata->material_type = Conductor;
      pzonedata->material = Conductor;
      try {pzonedata->mt = new MatConductor(zone[i].zonename, 
                          FormatElectrodeString(zone[i].zonelabel),&scale_unit);}
      catch(int){return 1;}
      pzonedata->fs = new ELData[pzonedata->node_num];
      pzonedata->aux = new ELAuxData[pzonedata->node_num];
      pzonedata->electrode.clear();
      for(int j=0;j<pzonedata->pzone->dasegment.size();j++)
      {
        int bc_index = pzonedata->pzone->dasegment[j].bc_index-1;
        if(bc[bc_index].BCType==OhmicContact||
            bc[bc_index].BCType==SchottkyContact||
            bc[bc_index].BCType==GateContact)
          pzonedata->electrode.push_back(bc_index);
      }
      if(pzonedata->electrode.size()>1)
      {
        gss_log.string_buf()<<"Electrode region "<<zone[i].zonename<<" has more than one Contact, not supported yet.\n";
        gss_log.record();
        return 1;
      }
      //pzonedata->report();
    }
    //Data for Vacuum zone
    else if(IsVacuum(zone[i].zonelabel))
    {
      zonedata[i] = new VacuumZone;
      VacuumZone *pzonedata = dynamic_cast< VacuumZone * >(zonedata[i]);
      pzonedata->pzone = &zone[i];
      pzonedata->zone_index = zone[i].zone_index;
      pzonedata->node_num = zone[i].davcell.size();
      pzonedata->tri_num  = zone[i].datri.size();
      pzonedata->material_type = Vacuum;
      pzonedata->material = Vacuum;
      try {pzonedata->mt = new MatVacuum(zone[i].zonename,
                          FormatVacuumString(zone[i].zonelabel),&scale_unit);}
      catch(int){return 1;}
      pzonedata->fs = new VacuumData[pzonedata->node_num];
      pzonedata->aux = new VacuumAuxData[pzonedata->node_num];
    }
    //PML data
    else if(IsPML(zone[i].zonelabel))
    {
      zonedata[i] = new PMLZone;
      PMLZone *pzonedata = dynamic_cast< PMLZone * >(zonedata[i]);
      pzonedata->pzone = &zone[i];
      pzonedata->zone_index = zone[i].zone_index;
      pzonedata->node_num = zone[i].davcell.size();
      pzonedata->tri_num  = zone[i].datri.size();
      pzonedata->material_type = PML;
      pzonedata->material = PML;
      try {pzonedata->mt = new MatPML(zone[i].zonename,
                          FormatPMLString(zone[i].zonelabel),&scale_unit);}
      catch(int){return 1;}
      pzonedata->fs = new PMLData[pzonedata->node_num];
      pzonedata->aux = new PMLAuxData[pzonedata->node_num];
    }
    else
    {
      gss_log.string_buf()<<"error, no such region type "<<zone[i].zonelabel<<" specified!\n";
      gss_log.record();
      return 1;
    }
  }

  // process inter_connect here
  for(int i=0; i<bc.size(); i++)
  {
    if(bc.Get_pointer(i)->BCType==OhmicContact)
    {
      OhmicBC *pbc = dynamic_cast <OhmicBC * > (bc.Get_pointer(i));
      if(pbc->connect_elec!="")
      {
        int flag=0;
        for(int j=0;j<bc.size();j++)
        {
          if(pbc->connect_elec==bc.Get_pointer(j)->psegment->label&&
              IsSemiconductor(zone[bc.Get_pointer(j)->psegment->zone_index].zonelabel))
          {
            flag=1;
            pbc->inner_connect=j;
            pbc->connect_zone = bc.Get_pointer(j)->psegment->zone_index;
            if(bc.Get_pointer(pbc->inner_connect)->BCType==OhmicContact)
            {
              OhmicBC *om_bc = dynamic_cast <OhmicBC * >(bc.Get_pointer(pbc->inner_connect));
              om_bc->inner_connect=i;
              om_bc->connect_zone=pbc->psegment->zone_index;
              om_bc->connect_elec=pbc->psegment->label;
              om_bc->pzonedata=zonedata[pbc->psegment->zone_index];
            }
            else
            {
              sprintf(log_buf,"\nCheck Inner Connect failed, BC %s should connect to an Ohmic/Schottky Electrode!\n",
                      bc.Get_pointer(pbc->inner_connect)->psegment->label);
              GSS_LOG();
              return 1;
            }
          }
        }
        if(!flag)
        {
          sprintf(log_buf,"\nSorry, I can't find a suitable inner conect for %s!\n",pbc->connect_elec.c_str());
          GSS_LOG();
          return 1;
        }
      }
    }
  }

  return 0;
}

int  BSolver::setup_init_data()
{
  for(int i=0; i<zone_num; i++)
    if(zonedata[i]->Init(&zone[i],LatticeTemp,&scale_unit)) return 1;
  return 0;
}


int  BSolver::setup_doping()
{
  for(int z=0; z<zone_num; z++)
    if(zonedata[z]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
      for(int i=0;i<pzonedata->node_num;i++)
      {
        double x = pzonedata->pzone->danode[i].x;
        double y = pzonedata->pzone->danode[i].y;
        for(int j=0;j<doping_func.size();j++)
        {
          double doping = doping_func[j]->profile(x,y);
          if(doping > 0)
            pzonedata->aux[i].Nd += doping;
          else
            pzonedata->aux[i].Na += fabs(doping);
        }
      }
    }
  return 0;
}


int  BSolver::import_doping_from_cgns(char* cgnsfile)
{
  for(int z=0; z<zone_num; z++)
  {
    if(zonedata[z]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
      if(pzonedata->import_doping(cgnsfile,&scale_unit)) return 1;
    }
  }
  return 0;
}


int  BSolver::import_mole_from_cgns(char* cgnsfile)
{
  for(int z=0; z<zone_num; z++)
  {
    if(IsSingleCompSemiconductor(zone[z].zonelabel))
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
      if(pzonedata->import_mole(cgnsfile,1)) return 1;
    }
  }
  return 0;
}

int  BSolver::setup_bc()
{
  //the bc init here
  Segment **segment = new Segment*[bc_num];
  for(int i=1; i<=bc_num;i++)
    segment[i-1] = Get_segment(i);
  if(bc.InitBC(bc_num,segment,zone_num,zone,zinterface,pcmdbuf,LatticeTemp,&scale_unit)) return 1;
  delete [] segment;

  //for bc marks are assigned to edge, some nodes may have two bc marks
  //assign this type of nodes to correct bc type
  for(int i=0;i<bc.size();i++)
  {
      if(bc[i].psegment->interface==-1)
      {
        int zone_index = bc[i].psegment->zone_index;
        for(int j=0;j<bc[i].psegment->node_array.size();j++)
        {
          int node_index = bc[i].psegment->node_array[j];
          int bc_index=zone[zone_index].danode[node_index].bc_index;
          if(bc[bc_index-1].BCType < bc[i].BCType)
          {
            zone[zone_index].danode[node_index].bc_index = i+1;
            zone[zone_index].davcell.GetPointer(node_index)->bc_index = i+1;
          }  
        }
      }
      else
      {
        int inf = bc[i].psegment->interface;
        int z1 = zinterface[inf].zone1;
        int z2 = zinterface[inf].zone2;
        if(z1==bc[i].psegment->zone_index)
          for(int j=0;j<zinterface[inf].index_array1.size();j++)
          {
            int node_index = zinterface[inf].index_array1[j];
            int bc_index=zone[z1].danode[node_index].bc_index;
            if(bc[bc_index-1].BCType < bc[i].BCType)
            {
              zone[z1].danode[node_index].bc_index = i+1;
              zone[z1].davcell.GetPointer(node_index)->bc_index = i+1;
            }  
          }
        if(z2==bc[i].psegment->zone_index)
          for(int j=0;j<zinterface[inf].index_array2.size();j++)
          {
            int node_index = zinterface[inf].index_array2[j];
            int bc_index=zone[z2].danode[node_index].bc_index;
            if(bc[bc_index-1].BCType < bc[i].BCType)
            {
              zone[z2].danode[node_index].bc_index = i+1;
              zone[z2].davcell.GetPointer(node_index)->bc_index = i+1;
            }
          }
      }
  }
  
  return 0;
}


void BSolver::reorder()
{

  for(int z=0; z<zone_num; z++)
  {
    int * old_to_new = reorder_zone(z);
    if(zonedata[z]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
      SemiData* fs_new = new SemiData[pzonedata->node_num];
      SemiAuxData* aux_new = new SemiAuxData[pzonedata->node_num];
      for(int i=0;i<pzonedata->node_num;i++)
      {
        fs_new[old_to_new[i]] = pzonedata->fs[i];
        aux_new[old_to_new[i]] = pzonedata->aux[i];
      }
      delete [] pzonedata->fs;
      delete [] pzonedata->aux;
      pzonedata->fs = fs_new;
      pzonedata->aux = aux_new;
    }
    if(zonedata[z]->material_type == Insulator)
    {
      ISZone *pzonedata = dynamic_cast<ISZone *> (zonedata[z]);
      ISData* fs_new = new ISData[zonedata[z]->node_num];
      ISAuxData* aux_new = new ISAuxData[pzonedata->node_num];
      for(int i=0;i<pzonedata->node_num;i++)
      {
        fs_new[old_to_new[i]] = pzonedata->fs[i];
        aux_new[old_to_new[i]] = pzonedata->aux[i];
      }
      delete [] pzonedata->fs;
      delete [] pzonedata->aux;
      pzonedata->fs = fs_new;
      pzonedata->aux = aux_new;
    }
    if(zonedata[z]->material_type == Conductor)
    {
      ElZone *pzonedata = dynamic_cast<ElZone *> (zonedata[z]);
      ELData* fs_new = new ELData[zonedata[z]->node_num];
      ELAuxData* aux_new = new ELAuxData[pzonedata->node_num];
      for(int i=0;i<pzonedata->node_num;i++)
      {
        fs_new[old_to_new[i]] = pzonedata->fs[i];
        aux_new[old_to_new[i]] = pzonedata->aux[i];
      }
      delete [] pzonedata->fs;
      delete [] pzonedata->aux;
      pzonedata->fs = fs_new;
      pzonedata->aux = aux_new;
    }
    if(zonedata[z]->material_type == Vacuum)
    {
      VacuumZone *pzonedata = dynamic_cast<VacuumZone *> (zonedata[z]);
      VacuumData* fs_new = new VacuumData[zonedata[z]->node_num];
      VacuumAuxData* aux_new = new VacuumAuxData[pzonedata->node_num];
      for(int i=0;i<pzonedata->node_num;i++)
      {
        fs_new[old_to_new[i]] = pzonedata->fs[i];
        aux_new[old_to_new[i]] = pzonedata->aux[i];
      }
      delete [] pzonedata->fs;
      delete [] pzonedata->aux;
      pzonedata->fs = fs_new;
      pzonedata->aux = aux_new;
    }
    if(zonedata[z]->material_type == PML)
    {
      PMLZone *pzonedata = dynamic_cast<PMLZone *> (zonedata[z]);
      PMLData* fs_new = new PMLData[zonedata[z]->node_num];
      PMLAuxData* aux_new = new PMLAuxData[pzonedata->node_num];
      for(int i=0;i<pzonedata->node_num;i++)
      {
        fs_new[old_to_new[i]] = pzonedata->fs[i];
        aux_new[old_to_new[i]] = pzonedata->aux[i];
      }
      delete [] pzonedata->fs;
      delete [] pzonedata->aux;
      pzonedata->fs = fs_new;
      pzonedata->aux = aux_new;
    }
    delete [] old_to_new;
  }
  build_zone();
}


int BSolver::build_least_squares()
{
  double a=0,b=0,c=0,w=0;
  for(int z=0; z<zone_num; z++)
    if(zonedata[z]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
      for(int i=0; i<zone[z].davcell.size(); i++)
      {
        a=0,b=0,c=0,w=0;
        VoronoiCell *pcell = zone[z].davcell.GetPointer(i);
        for(int k=0;k<pcell->nb_num;k++)
        {
          const VoronoiCell *ncell = zone[z].davcell.GetPointer(pcell->nb_array[k]);
          double dx = ncell->x - pcell->x;
          double dy = ncell->y - pcell->y;
          w=1.0/sqrt(dx*dx+dy*dy);
          a+=w*w*dx*dx;
          b+=w*w*dx*dy;
          c+=w*w*dy*dy;
        }
        if(pcell->bc_index&&bc[pcell->bc_index-1].BCType==InsulatorInterface)
        {
          InsulatorInterfaceBC *pbc;
          pbc = dynamic_cast<InsulatorInterfaceBC*>(bc.Get_pointer(pcell->bc_index-1));
          int n_zone = pbc->pinterface->Find_neighbor_zone_index(z);
          int n_node = pbc->pinterface->Find_neighbor_node_index(z,i);
          ISZone * pz = dynamic_cast<ISZone *>(zonedata[n_zone]);
          const VoronoiCell* dcell = pz->pzone->davcell.GetPointer(n_node);
          for(int k=1;k<dcell->nb_num-1;k++)
          {
            const VoronoiCell *ncell = pz->pzone->davcell.GetPointer(dcell->nb_array[k]);
            double dx = ncell->x - dcell->x;
            double dy = ncell->y - dcell->y;
            w=1.0/sqrt(dx*dx+dy*dy);
            a+=w*w*dx*dx;
            b+=w*w*dx*dy;
            c+=w*w*dy*dy;
          }
        }
        pcell->sa = a/(a*c-b*b);
        pcell->sb = b/(a*c-b*b);
        pcell->sc = c/(a*c-b*b);
      }
    }
  return 0;
}

int  VacuumZone::Init(ZONE* zone, double Tl,PhysicalUnitScale *scale)
{
  for(int i=0;i<node_num; i++)
  {
    aux[i].eps = mt->eps0;
    aux[i].mu  = mt->mu0;
  }
  return 0;
}

int  PMLZone::Init(ZONE* zone, double Tl,PhysicalUnitScale *scale)
{
  for(int i=0;i<node_num; i++)
  {
    aux[i].eps = mt->eps0;
    aux[i].mu  = mt->mu0;
  }
  return 0;
}

int  ISZone::Init(ZONE* zone, double Tl,PhysicalUnitScale *scale)
{
  for(int i=0;i<node_num; i++)
  {
    mt->mapping(&pzone->danode[i],&aux[i],0);
    aux[i].Ex = aux[i].Ey = 0;
    aux[i].affinity = mt->basic->Affinity(Tl);
    aux[i].density =  mt->basic->Density(Tl);
    aux[i].eps = mt->eps0*mt->basic->Permittivity();
    aux[i].mu  = mt->mu0*mt->basic->Permeability();
    fs[i].T =  Tl;
    fs[i].P =  -aux[i].affinity;
  }
  return 0;
}

int  ElZone::Init(ZONE* zone, double Tl,PhysicalUnitScale *scale)
{
  for(int i=0;i<node_num; i++)
  {
    aux[i].affinity = mt->basic->Affinity(Tl);
    aux[i].density =  mt->basic->Density(Tl);
    aux[i].eps = mt->eps0*mt->basic->Permittivity();
    aux[i].mu  = mt->mu0*mt->basic->Permeability();
    fs[i].T =  Tl;
    fs[i].P =  -aux[i].affinity;
  }
  return 0;
}


int  SMCZone::Init(ZONE* zone, double Tl,PhysicalUnitScale *scale)
{
  PetscScalar e  =  mt->e;
  PetscScalar kb =  mt->kb;
  for(int i=0;i<node_num; i++)
  {
    fs[i].T =  Tl;
    fs[i].Tn=fs[i].Tp=fs[i].T;
    mt->mapping(&pzone->danode[i],&aux[i],0);
    double nie = mt->band->nie(fs[i].T);
    double Nc  = mt->band->Nc(fs[i].T);
    if( aux[i].Net_doping() > 0 )   //n-type
    {
      fs[i].n=(aux[i].Net_doping()+sqrt(pow(aux[i].Net_doping(),2)+4*nie*nie))/2;
      fs[i].p=nie*nie/fs[i].n;
    }
    else                       //p-type
    {
      fs[i].p=(-aux[i].Net_doping()+sqrt(pow(-aux[i].Net_doping(),2)+4*nie*nie))/2;
      fs[i].n=nie*nie/fs[i].p;
    }
    aux[i].affinity = mt->basic->Affinity(fs[i].T);
    aux[i].density =  mt->basic->Density(fs[i].T);
    aux[i].eps = mt->eps0*mt->basic->Permittivity();
    aux[i].mu  = mt->mu0*mt->basic->Permeability();
    aux[i].Eg     =   mt->band->Eg(fs[i].T)-mt->band->EgNarrow(fs[i].T);
    aux[i].Nc     =   mt->band->Nc(fs[i].T);
    aux[i].Nv     =   mt->band->Nv(fs[i].T);
    fs[i].P       =   kb*fs[i].T/e*asinh(aux[i].Net_doping()/(2*nie)) - kb*fs[i].T/e*log(Nc/nie)
                      - aux[i].affinity - mt->band->EgNarrowToEc(fs[i].T);
    fs[i].Eqc = -(e*fs[i].P + aux[i].affinity + mt->band->EgNarrowToEc(fs[i].T) + kb*fs[i].T*log(aux[i].Nc));
    fs[i].Eqv = -(e*fs[i].P + aux[i].affinity - mt->band->EgNarrowToEv(fs[i].T) - kb*fs[i].T*log(aux[i].Nv)+ aux[i].Eg);
    aux[i].mun = aux[i].mup = 0;
  }
  return 0;
}

void SMCZone::report()
{
  gss_log.string_buf()<<"\nThis is a semiconductor Zone.\n";
  gss_log.record();
}

void ISZone::report()
{
  gss_log.string_buf()<<"\nThis is a insulator Zone.\n";
  gss_log.record();
}

void ElZone::report()
{
  gss_log.string_buf()<<"\nThis is a electrode.\n";
  gss_log.record();
}

void VacuumZone::report()
{
  gss_log.string_buf()<<"\nThis is vacuum area.\n";
  gss_log.record();
}

void PMLZone::report()
{
  gss_log.string_buf()<<"\nThis is PML.\n";
  gss_log.record();
}



