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
/*  Last update: Feb 15 2006                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "matdefine.h"
#include "material.h"
#include "log.h"
#include <stdlib.h>
#include <math.h>
#include <dlfcn.h>



int init_PMI_Semiconductor(CmdBuf *pcmdbuf, PMISDefine & PMIS_define)
{
  //process PMIS (Physical Model Interface of Semiconductor) define here
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
  {
    if(pcmdbuf->is_current_cmd("PMIS"))   //it is a set structure
    {
      list<Cmd>::iterator    pcmd = pcmdbuf->get_current_cmd();
      if(!pcmd->allowed_args(3,"region","mobility","ii.model"))
      {
        sprintf(log_buf,"line %d PMIS: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
        GSS_LOG();
        return 1;
      }
      PMIS_define.Material.push_back(pcmd->get_string("material",0,""));
      PMIS_define.Region.push_back(pcmd->get_string("region",0,""));
      PMIS_define.Mobility.push_back(pcmd->get_string("mobility",0,"Default"));
      PMIS_define.IIModel.push_back(pcmd->get_string("ii.model",0,"Default"));
      PMIS_define.OpticalModel.push_back(pcmd->get_string("optical.model",0,"Default"));
    }
  }
  return 0;
}


PhysicalConst::PhysicalConst(PhysicalUnitScale *unit_scale)
{
  kb   = 1.3806503e-23*unit_scale->s_joule/unit_scale->s_kelvin;
  e    = 1.602176462e-19*unit_scale->s_coulomb;
  me   = 9.10938188e-31*unit_scale->s_kg;
  eps0 = 8.854187818e-12*unit_scale->s_coulomb/unit_scale->s_volt/unit_scale->s_meter;
  mu0  = 12.56637061e-7*pow(unit_scale->s_second,2)/unit_scale->s_coulomb*unit_scale->s_volt/unit_scale->s_meter;
  h    = 6.62606876e-34*unit_scale->s_joule*unit_scale->s_second;
  hbar = 1.054571596e-34*unit_scale->s_joule*unit_scale->s_second;
}


MatSemiconductor::MatSemiconductor(const char * region_name, const char *material,
                                   PMISDefine & PMIS_define, PhysicalUnitScale *unit_scale)
    :PhysicalConst(unit_scale)
{
  pscale = unit_scale;
  char filename[128];
  char model_fun_name[128];
  int  region_index = -1;

  PMIS_Environment env(pscale,&pnode,&plocal_semi_data,&time);
  for(int i=0;i<PMIS_define.Region.size();i++)
  {
    if(PMIS_define.Region[i] == region_name) region_index = i;
  }
  sprintf(filename,"%s/lib/lib%s.so",getenv("GSS_DIR"),material);
  dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL) {fprintf(stderr,"Open material file lib%s.so error\n%s\n",material,dlerror()); throw ERROR;}
  PMIS_BasicParameter*(*wbasic)    (const PMIS_Environment& env);
  PMIS_BandStructure* (*wband)     (const PMIS_Environment& env);
  PMIS_Mobility*      (*wmob)      (const PMIS_Environment& env);
  PMIS_Avalanche*     (*wgen)      (const PMIS_Environment& env);
  PMIS_Thermal*       (*wthermal)  (const PMIS_Environment& env);
  PMIS_Optical*       (*woptical)  (const PMIS_Environment& env);
  
  //init AD indepedent variable set routine
  set_ad_num = (void* (*) (const unsigned int))dlsym(dp,"set_ad_number");
  if(dlerror()) {fprintf(stderr,"Open PMIS AD_SET_VARIABLE function error!\n");throw ERROR;}
  
  // init basic parameters for the material
  sprintf(model_fun_name,"PMIS_%s_BasicParameter_Default",material);
  wbasic = (PMIS_BasicParameter* (*) (const PMIS_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIS %s BasicParameter function error!\n",material);throw ERROR;}

  // init band structure model
  sprintf(model_fun_name,"PMIS_%s_BandStructure_Default",material);
  wband =  (PMIS_BandStructure* (*) (const PMIS_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIS %s BandStructure function error!\n",material);throw ERROR;}

  // init mobility model
  if(region_index != -1)
    sprintf(model_fun_name,"PMIS_%s_Mob_%s",material,PMIS_define.Mobility[region_index].c_str());
  else
    sprintf(model_fun_name,"PMIS_%s_Mob_%s",material,"Default");
  wmob  =  (PMIS_Mobility* (*) (const PMIS_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIS Mobility function %s error!\n",model_fun_name);throw ERROR;}

  // init Avalanche generation model
  if(region_index != -1)
    sprintf(model_fun_name,"PMIS_%s_Avalanche_%s",material,PMIS_define.IIModel[region_index].c_str());
  else
    sprintf(model_fun_name,"PMIS_%s_Avalanche_Default",material);
  wgen  =  (PMIS_Avalanche* (*) (const PMIS_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIS %s Avalanche function error!\n",material);throw ERROR;}

  // init Thermal model for lattice temperature equation
  sprintf(model_fun_name,"PMIS_%s_Thermal_Default",material);
  wthermal  = (PMIS_Thermal* (*) (const PMIS_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIS %s Thermal function error!\n",material);throw ERROR;}

  // init optical data
  if(region_index != -1)
    sprintf(model_fun_name,"PMIS_%s_Optical_%s",material,PMIS_define.OpticalModel[region_index].c_str());
  else
    sprintf(model_fun_name,"PMIS_%s_Optical_Default",material);
  woptical  = (PMIS_Optical* (*) (const PMIS_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIS %s Optical function error!\n",material);throw ERROR;}
  
  basic = wbasic(env);
  band  = wband(env);
  mob   = wmob(env);
  gen   = wgen(env);
  thermal  = wthermal(env);
  optical  = woptical(env);
}

void MatSemiconductor::mapping(Node* node, SemiAuxData *data, PetscScalar clock)
{
  plocal_semi_data = data;
  pnode            = node;
  time             = clock;
}

MatSemiconductor::~MatSemiconductor()
{
  delete basic;
  delete band;
  delete mob;
  delete gen;
  delete thermal;
  delete optical;
  dlclose(dp);
}




//------------------------------------------------------------------------------------------------------

MatInsulator::MatInsulator(const char * region_name, const char *material,
                           PhysicalUnitScale *unit_scale)
    :PhysicalConst(unit_scale)
{
  pscale = unit_scale;
  char filename[128];
  char model_fun_name[128];
  char *error;
  PMII_Environment env(pscale,&pnode,&plocal_is_data,&time);
  sprintf(filename,"%s/lib/lib%s.so",getenv("GSS_DIR"),material);
  dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL) {fprintf(stderr,"open material file lib%s.so error\n",material); throw ERROR;}
  PMII_BasicParameter*(*wbasic)    (const PMII_Environment& env);
  PMII_Thermal*       (*wthermal)  (const PMII_Environment& env);
  PMII_Optical*       (*woptical)  (const PMII_Environment& env);
  // init basic parameters for the material
  sprintf(model_fun_name,"PMII_%s_BasicParameter_Default",material);
  wbasic = (PMII_BasicParameter* (*) (const PMII_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMII %s BasicParameter function error!\n",material);throw ERROR;}

  // init Thermal model for lattice temperature equation
  sprintf(model_fun_name,"PMII_%s_Thermal_Default",material);
  wthermal  = (PMII_Thermal* (*) (const PMII_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMII %s Thermal function error!\n",material);throw ERROR;}

  // init optical data
  sprintf(model_fun_name,"PMII_%s_Optical_Default",material);
  woptical  = (PMII_Optical* (*) (const PMII_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMII %s Optical function error!\n",material);throw ERROR;}
  
  basic    = wbasic(env);
  thermal  = wthermal(env);
  optical  = woptical(env);
}

MatInsulator::~MatInsulator()
{
  delete basic;
  delete thermal;
  delete optical;
  dlclose(dp);
}

void MatInsulator::mapping(Node* node, ISAuxData *data, PetscScalar clock)
{
  plocal_is_data   = data;
  pnode            = node;
  time             = clock;
}

//--------------------------------------------------------------------------------------------------------------
MatConductor::MatConductor(const char * region_name, const char *material,
                           PhysicalUnitScale *unit_scale):
    PhysicalConst(unit_scale)
{
  pscale = unit_scale;
  char filename[128];
  char model_fun_name[128];
  char *error;
  PMIC_Environment env(pscale,&pnode,&plocal_el_data,&time);
  sprintf(filename,"%s/lib/lib%s.so",getenv("GSS_DIR"),material);
  dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL) {fprintf(stderr,"open material file lib%s.so error\n",material); throw ERROR;}
  PMIC_BasicParameter*(*wbasic)    (const PMIC_Environment& env);
  PMIC_Thermal*       (*wthermal)  (const PMIC_Environment& env);
  PMIC_Optical*       (*woptical)  (const PMIC_Environment& env);
  // init basic parameters for the material
  sprintf(model_fun_name,"PMIC_%s_BasicParameter_Default",material);
  wbasic = (PMIC_BasicParameter* (*) (const PMIC_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIC %s BasicParameter function error!\n",material);throw ERROR;}

  // init Thermal model for lattice temperature equation
  sprintf(model_fun_name,"PMIC_%s_Thermal_Default",material);
  wthermal  = (PMIC_Thermal* (*) (const PMIC_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIC %s Thermal function error!\n",material);throw ERROR;}

  // init optical data
  sprintf(model_fun_name,"PMIC_%s_Optical_Default",material);
  woptical  = (PMIC_Optical* (*) (const PMIC_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIC %s Optical function error!\n",material);throw ERROR;}
  
  basic    = wbasic(env);
  thermal  = wthermal(env);
  optical  = woptical(env);
}

MatConductor::~MatConductor()
{
  delete basic;
  delete thermal;
  delete optical;
  dlclose(dp);
}

void MatConductor::mapping(Node* node, ELAuxData *data, PetscScalar clock)
{
  plocal_el_data   = data;
  pnode            = node;
  time             = clock;
}

//----------------------------------------------------------------------------------------------------------

MatVacuum::MatVacuum(const char * region_name, const char *material,
                     PhysicalUnitScale *unit_scale)
    :PhysicalConst(unit_scale)
{
  pscale = unit_scale;
  char filename[128];
  char model_fun_name[128];
  char *error;
  PMIV_Environment env(pscale,&pnode,&plocal_vac_data,&time);
  sprintf(filename,"%s/lib/lib%s.so",getenv("GSS_DIR"),material);
  dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL) {fprintf(stderr,"open material file lib%s.so error\n",material); throw ERROR;}
  PMIV_BasicParameter*(*wbasic)    (const PMIV_Environment& env);
  PMIV_Thermal*       (*wthermal)  (const PMIV_Environment& env);
  // init basic parameters for the material
  sprintf(model_fun_name,"PMIV_%s_BasicParameter_Default",material);
  wbasic = (PMIV_BasicParameter* (*) (const PMIV_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIV %s BasicParameter function error!\n",material);throw ERROR;}

  // init Thermal model for lattice temperature equation
  sprintf(model_fun_name,"PMIV_%s_Thermal_Default",material);
  wthermal  = (PMIV_Thermal* (*) (const PMIV_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIV %s Thermal function error!\n",material);throw ERROR;}

  basic    = wbasic(env);
  thermal  = wthermal(env);
}

MatVacuum::~MatVacuum()
{
  delete basic;
  delete thermal;
  dlclose(dp);
}

void MatVacuum::mapping(Node* node, VacuumAuxData *data, PetscScalar clock)
{
  plocal_vac_data   = data;
  pnode            = node;
  time             = clock;
}

//------------------------------------------------------------------------------------------------------------

MatPML::MatPML(const char * region_name, const char *material, PhysicalUnitScale *unit_scale)
    :PhysicalConst(unit_scale)
{
  pscale = unit_scale;
  char filename[128];
  char model_fun_name[128];
  char *error;
  PMIP_Environment env(pscale,&pnode,&plocal_pml_data,&time);
  sprintf(filename,"%s/lib/lib%s.so",getenv("GSS_DIR"),material);
  dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL) {fprintf(stderr,"open material file lib%s.so error\n",material); throw ERROR;}
  PMIP_BasicParameter*(*wbasic)    (const PMIP_Environment& env);
  PMIP_Thermal*       (*wthermal)  (const PMIP_Environment& env);
  // init basic parameters for the material
  sprintf(model_fun_name,"PMIP_%s_BasicParameter_Default",material);
  wbasic = (PMIP_BasicParameter* (*) (const PMIP_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIP %s BasicParameter function error!\n",material);throw ERROR;}

  // init Thermal model for lattice temperature equation
  sprintf(model_fun_name,"PMIP_%s_Thermal_Default",material);
  wthermal  = (PMIP_Thermal* (*) (const PMIP_Environment& env))dlsym(dp,model_fun_name);
  if(dlerror()) {fprintf(stderr,"Open PMIP %s Thermal function error!\n",material);throw ERROR;}

  basic    = wbasic(env);
  thermal  = wthermal(env);
}

MatPML::~MatPML()
{
  delete basic;
  delete thermal;
  dlclose(dp);
}

void MatPML::mapping(Node* node, PMLAuxData *data, PetscScalar clock)
{
  plocal_pml_data   = data;
  pnode            = node;
  time             = clock;
}


