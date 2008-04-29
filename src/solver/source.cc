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
/*  Last update: May 11, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "cmdbuf.h"
#include "bsolver.h"
#include "log.h"
#include <string.h>


VSource* SetVDC(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  VSource* vsrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw VSOURCE_NO_ID;
  if(!pcmd->allowed_args(4,"type","id","vconst","tdelay"))
        throw VSOURCE_UNKNOW_PARAMETER;
  double v = pcmd->get_number("vconst",0,0)*scale->s_volt;
  double Tdelay = pcmd->get_number("tdelay",0,0)*scale->s_second;
  vsrc=new VDC(Tdelay,v);
  strcpy(vsrc->label,ID);
  
  gss_log.string_buf()<<"\nVDC    : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" "
                      <<"v="<<v/scale->s_volt<<"\n";
  gss_log.record();
  return vsrc;
}


VSource* SetVSIN(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  VSource* vsrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw VSOURCE_NO_ID;
  if(!pcmd->allowed_args(7,"type","id","vconst","vamp","tdelay","freq","alpha"))
        throw VSOURCE_UNKNOW_PARAMETER;
  double v0 = pcmd->get_number("vconst",0,0.0)*scale->s_volt;
  double vamp = pcmd->get_number("vamp",0,0.0)*scale->s_volt;
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;
  double freq = pcmd->get_number("freq",0,0.0)*1.0/scale->s_second;
  double alpha= pcmd->get_number("alpha",0,0.0)*1.0/scale->s_second;


  vsrc=new VSIN(Tdelay,v0,vamp,freq,alpha);
  strcpy(vsrc->label,ID);

  gss_log.string_buf()<<"\nVSIN   : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" "
		      <<"v0="<<v0/scale->s_volt<<" "
                      <<"vamp="<<vamp/scale->s_volt<<" "
                      <<"freq="<<freq*scale->s_second<<"\n"
                      <<"alpha="<<alpha*scale->s_second<<"\n";
  gss_log.record();

  return vsrc;
}


VSource* SetVEXP(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  VSource* vsrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw VSOURCE_NO_ID;
  if(!pcmd->allowed_args(8,"type","id","tdelay","trc","tfd","tfc","vlo","vhi"))
        throw VSOURCE_UNKNOW_PARAMETER;
  double Tdelay = pcmd->get_number("tdelay",0,0)*scale->s_second;
  double trc = pcmd->get_number("trc",0,0.0)*scale->s_second;
  double tfd = pcmd->get_number("tfd",0,0.0)*scale->s_second;
  double tfc = pcmd->get_number("tfc",0,0.0)*scale->s_second;
  double vlo = pcmd->get_number("vlo",0,0.0)*scale->s_volt;
  double vhi = pcmd->get_number("vhi",0,0.0)*scale->s_volt;

  vsrc=new VEXP(Tdelay,vlo,vhi,trc,tfd,tfc);
  strcpy(vsrc->label,ID);

  gss_log.string_buf()<<"\nVEXP   : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" " 
                      <<"trc="<<trc/scale->s_second<<" " 
                      <<"tfd="<<tfd/scale->s_second<<" " 
                      <<"tfc="<<tfc/scale->s_second<<"\n\t\t" 
                      <<"vlo="<<vlo/scale->s_volt<<" " 
                      <<"vhi="<<vhi/scale->s_volt<<"\n";
  gss_log.record();

  return vsrc;
}

VSource* SetVPULSE(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  VSource* vsrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw VSOURCE_NO_ID;
  if(!pcmd->allowed_args(9,"type","id","tdelay","tr","tf","pw","pr","vlo","vhi"))
        throw VSOURCE_UNKNOW_PARAMETER;
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;
  double tr = pcmd->get_number("tr",0,1e-12)*scale->s_second;
  double tf = pcmd->get_number("tf",0,1e-12)*scale->s_second;
  double pw = pcmd->get_number("pw",0,0.0)*scale->s_second;
  double pr = pcmd->get_number("pr",0,0.0)*scale->s_second;
  double vlo = pcmd->get_number("vlo",0,0.0)*scale->s_volt;
  double vhi = pcmd->get_number("vhi",0,0.0)*scale->s_volt;

  vsrc=new VPULSE(Tdelay,vlo,vhi,tr,tf,pw,pr);
  strcpy(vsrc->label,ID);

  gss_log.string_buf()<<"\nVPULSE : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" " 
                      <<"tr="<<tr/scale->s_second<<" " 
                      <<"tf="<<tf/scale->s_second<<" " 
                      <<"pw="<<pw/scale->s_second<<" " 
                      <<"pr="<<pr/scale->s_second<<"\n\t\t" 
                      <<"vlo="<<vlo/scale->s_volt<<" "
                      <<"vhi="<<vhi/scale->s_volt<<"\n";
  gss_log.record();

  return vsrc;
}

VSource* SetVSHELL(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  VSource* vsrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw VSOURCE_NO_ID;
  if(!pcmd->allowed_args(4,"type","id","dll","func"))
        throw VSOURCE_UNKNOW_PARAMETER;
  char * filename = pcmd->get_string("dll",0,"");
  char * funcname = pcmd->get_string("func",0,"");
  void * dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL)
  {
    gss_log.string_buf()<<"Open dynamic libaray "<<filename<<" error!\n";
    gss_log.record();
    throw VSOURCE_UNKNOW_PARAMETER;
  }
  void *fp = dlsym(dp,funcname);
  if(fp==NULL)
  {
        gss_log.string_buf()<<"The dynamic library "<<filename<<" don't have function "<<funcname<<"\n";
        gss_log.record();
        throw VSOURCE_UNKNOW_PARAMETER;
  }
  vsrc= new VSHELL(dp,fp,scale->s_second,scale->s_volt);
  strcpy(vsrc->label,ID);
  gss_log.string_buf()<<"\nVSHELL : "<<ID<<" load from "<<filename<<"\n";
  gss_log.record();
  return vsrc;
}


ISource* SetIDC(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  ISource* isrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw ISOURCE_NO_ID;
  if(!pcmd->allowed_args(4,"type","id","iconst","tdelay"))
        throw ISOURCE_UNKNOW_PARAMETER;
  double i = pcmd->get_number("iconst",0,0)*scale->s_mA;
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;

  isrc=new IDC(Tdelay,i);
  strcpy(isrc->label,ID);
  gss_log.string_buf()<<"\nIDC    : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" "
                      <<"i="<<i/scale->s_mA<<"\n";
  gss_log.record();
  return isrc;
}

ISource* SetISIN(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  ISource* isrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw ISOURCE_NO_ID;
  if(!pcmd->allowed_args(5,"type","id","iamp","tdelay","freq"))
        throw ISOURCE_UNKNOW_PARAMETER;
  double iamp = pcmd->get_number("iamp",0,0.0)*scale->s_mA;
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;
  double freq = pcmd->get_number("freq",0,0.0)*1.0/scale->s_second;

  isrc=new ISIN(Tdelay,iamp,freq);
  strcpy(isrc->label,ID);

  gss_log.string_buf()<<"\nISIN   : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" "
                      <<"iamp="<<iamp/scale->s_mA<<" "
                      <<"freq="<<freq*scale->s_second<<"\n";
  gss_log.record();
  return isrc;
}

ISource* SetIEXP(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  ISource* isrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw ISOURCE_NO_ID;
  if(!pcmd->allowed_args(8,"type","id","tdelay","trc","tfd","tfc","ilo","ihi"))
        throw ISOURCE_UNKNOW_PARAMETER;
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;
  double trc = pcmd->get_number("trc",0,0.0)*scale->s_second;
  double tfd = pcmd->get_number("tfd",0,0.0)*scale->s_second;
  double tfc = pcmd->get_number("tfc",0,0.0)*scale->s_second;
  double ilo = pcmd->get_number("ilo",0,0.0)*scale->s_mA;
  double ihi = pcmd->get_number("ihi",0,0.0)*scale->s_mA;
  isrc=new IEXP(Tdelay,ilo,ihi,trc,tfd,tfc);
  strcpy(isrc->label,ID);

  gss_log.string_buf()<<"\nIEXP   : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" " 
                      <<"trc="<<trc/scale->s_second<<" " 
                      <<"tfd="<<tfd/scale->s_second<<" " 
                      <<"tfc="<<tfc/scale->s_second<<"\n\t\t" 
                      <<"ilo="<<ilo/scale->s_mA<<" " 
                      <<"ihi="<<ihi/scale->s_mA<<"\n";
  gss_log.record();
  return isrc;
}

ISource* SetIPULSE(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  ISource* isrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw ISOURCE_NO_ID;
  if(!pcmd->allowed_args(9,"type","id","tdelay","tr","tf","pw","pr","ilo","ihi"))
        throw ISOURCE_UNKNOW_PARAMETER;
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;
  double tr = pcmd->get_number("tr",0,1e-12)*scale->s_second;
  double tf = pcmd->get_number("tf",0,1e-12)*scale->s_second;
  double pw = pcmd->get_number("pw",0,0.0)*scale->s_second;
  double pr = pcmd->get_number("pr",0,0.0)*scale->s_second;
  double ilo = pcmd->get_number("ilo",0,0.0)*scale->s_mA;
  double ihi = pcmd->get_number("ihi",0,0.0)*scale->s_mA;

  isrc=new IPULSE(Tdelay,ilo,ihi,tr,tf,pw,pr);
  strcpy(isrc->label,ID);

  gss_log.string_buf()<<"\nIPULSE : "<<ID<<"\t"
                      <<"td="<<Tdelay/scale->s_second<<" " 
                      <<"tr="<<tr/scale->s_second<<" " 
                      <<"tf="<<tf/scale->s_second<<" " 
                      <<"pw="<<pw/scale->s_second<<" " 
                      <<"pr="<<pr/scale->s_second<<"\n\t\t" 
                      <<"ilo="<<ilo/scale->s_mA<<" "
                      <<"ihi="<<ihi/scale->s_mA<<"\n";
  gss_log.record();
  return isrc;
}


ISource* SetISHELL(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  ISource* isrc;
  char   *ID;
  if(pcmd->is_arg_exist("id"))
    ID = pcmd->get_string("id",0,"");
  else
    throw ISOURCE_NO_ID;
  if(!pcmd->allowed_args(4,"type","id","dll","func"))
        throw ISOURCE_UNKNOW_PARAMETER;
  char * filename = pcmd->get_string("dll",0,"");
  char * funcname = pcmd->get_string("func",0,"");
  void * dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL)
  {
    gss_log.string_buf()<<"Open dynamic libaray "<<filename<<" error!\n";
    gss_log.record();
    throw ISOURCE_UNKNOW_PARAMETER;
  }
  void *fp = dlsym(dp,funcname);
  if(fp==NULL)
  {
        gss_log.string_buf()<<"The dynamic library "<<filename<<" don't have function "<<funcname<<"\n";
        gss_log.record();
        throw ISOURCE_UNKNOW_PARAMETER;
  }
  isrc= new ISHELL(dp,fp,scale->s_second,scale->s_mA);
  strcpy(isrc->label,ID);
  gss_log.string_buf()<<"\nISHELL : "<<ID<<" load from "<<filename<<"\n";
  gss_log.record();
  return isrc;
}


LSource* SetLightUNIFORM(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  LSource* lsrc;  
  
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;
  double power = pcmd->get_number("power",0,1.0);
  lsrc=new LSource_UNIFORM(Tdelay,power);
  gss_log.string_buf()<<"\nUNIFORM Light   :\t"
                      <<"Tdelay="<<Tdelay/scale->s_second<<"  "
                      <<"Power="<<power<<"\n ";                    
  gss_log.record();
  return lsrc;
}

LSource* SetLightGAUSSIAN(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  LSource* lsrc; 
  
  double tpeak = pcmd->get_number("tpeak",0,1e-12)*scale->s_second;
  double tradius = pcmd->get_number("tradius",0,4e-12)*scale->s_second;
  double power = pcmd->get_number("power",0,1.0);
  lsrc=new LSource_GAUSSIAN(tpeak,tradius,power);
  
  gss_log.string_buf()<<"\nGAUSSIAN Light  : \t"
                      <<"Tpeak="<<tpeak/scale->s_second<<" " 
                      <<"Tradius="<<tradius/scale->s_second<<"  "
                      <<"Power="<<power<<"\n ";  
                      
  gss_log.record();
  return lsrc;
}

LSource* SetLightPULSE(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  LSource* lsrc;
    
  double Tdelay = pcmd->get_number("tdelay",0,0.0)*scale->s_second;
  double tr = pcmd->get_number("tr",0,1e-15)*scale->s_second;
  double tf = pcmd->get_number("tf",0,1e-15)*scale->s_second;
  double pw = pcmd->get_number("pw",0,0.0)*scale->s_second;
  double pr = pcmd->get_number("pr",0,0.0)*scale->s_second; 
  double powerlo = pcmd->get_number("powerlo",0,0.0);
  double powerhi = pcmd->get_number("powerhi",0,1.0);
  
  lsrc=new LSource_PULSE(Tdelay,tr,tf,pw,pr,powerlo,powerhi);
  
  gss_log.string_buf()<<"\nPULSE Light  : \t"
                      <<"Tdelay="<<Tdelay/scale->s_second<<" " 
                      <<"tr="<<tr/scale->s_second<<" " 
                      <<"tf="<<tf/scale->s_second<<" " 
                      <<"pw="<<pw/scale->s_second<<" " 
                      <<"pr="<<pr/scale->s_second<<"\n\t\t" 
                      <<"powerlo="<<powerlo<<" "
                      <<"powerhi="<<powerhi<<"\n";
                      
  gss_log.record();
  return lsrc;
}


LSource* SetLightSHELL(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale)
{
  LSource* lsrc;
  if(!pcmd->allowed_args(3,"type","dll","func"))
        throw LSOURCE_UNKNOW_PARAMETER;
  char * filename = pcmd->get_string("dll",0,"");
  char * funcname = pcmd->get_string("func",0,"");
  void * dp = dlopen(filename,RTLD_LAZY);
  if(dp==NULL)
  {
    gss_log.string_buf()<<"Open dynamic libaray "<<filename<<" error!\n";
    gss_log.record();
    throw LSOURCE_UNKNOW_PARAMETER;
  }
  void *fp = dlsym(dp,funcname);
  if(fp==NULL)
  {
        gss_log.string_buf()<<"The dynamic library "<<filename<<" don't have function "<<funcname<<"\n";
        gss_log.record();
        throw LSOURCE_UNKNOW_PARAMETER;
  }
  lsrc= new LSource_SHELL(dp,fp,scale->s_second);
  gss_log.string_buf()<<"\nLSHELL load from "<<filename<<"\n";
  gss_log.record();
  return lsrc;
}


int BSolver::init_esource(list<Cmd> &cmdlist)
{
  list<Cmd>::iterator  pcmd;
  //setup voltage, current source and light source
  gss_log.string_buf()<<"\nSetting each voltage, current source and find light source here...\n";
  gss_log.record();

  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
  {
    if(pcmdbuf->is_current_cmd("VSOURCE"))   // It's a VSOURCE card
    {
      pcmd = pcmdbuf->get_current_cmd();
      if(pcmd->is_arg_exist("type"))
      {
        try
        {
          if(pcmd->is_arg_value("type","VDC"))
            {vsrc.push_back(SetVDC(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","VSIN"))
            {vsrc.push_back(SetVSIN(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","VEXP"))
            {vsrc.push_back(SetVEXP(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","VPULSE"))
            {vsrc.push_back(SetVPULSE(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","VSHELL"))
            {vsrc.push_back(SetVSHELL(pcmd,&scale_unit));}
          else
          {
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" VSOURCE: no such source Type!\n";
               gss_log.record();
               return 1;
          }
        }
        catch (VSOURCE_CARD_ERROR err)
        {
            if(err == VSOURCE_NO_ID)
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" VSOURCE: you should give an ID!\n";
            else if(err==VSOURCE_UNKNOW_PARAMETER)
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" VSOURCE: unrecognized parameter(s)!\n";
            gss_log.record();
          
            return 1;
        }
      }
    }
    else if(pcmdbuf->is_current_cmd("ISOURCE"))  //it is a isource structure
    {
      pcmd = pcmdbuf->get_current_cmd();
      if(pcmd->is_arg_exist("type"))
      {
        try
        {
          if(pcmd->is_arg_value("type","IDC"))
            {isrc.push_back(SetIDC(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","ISIN"))
            {isrc.push_back(SetISIN(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","IEXP"))
            {isrc.push_back(SetIEXP(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","IPULSE"))
            {isrc.push_back(SetIPULSE(pcmd,&scale_unit));}
          else if(pcmd->is_arg_value("type","ISHELL"))
            {isrc.push_back(SetISHELL(pcmd,&scale_unit));}
          else
          {
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" ISOURCE: no such source Type!\n";
               gss_log.record();
               return 1;
          }
        }
        catch (ISOURCE_CARD_ERROR err)
        {
                
            if(err == ISOURCE_NO_ID)
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" ISOURCE: you should give an ID!\n";
            else if(err == ISOURCE_UNKNOW_PARAMETER)
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" ISOURCE: unrecognized parameter(s)!\n";
            gss_log.record();

            return 1;
        }
      }
    }
    else if(pcmdbuf->is_current_cmd("LSOURCE"))  //there should be a light source structure
    {
      pcmd = pcmdbuf->get_current_cmd();
      if(pcmd->is_arg_exist("type"))
      {
        try
        {
          if(pcmd->is_arg_value("type","UNIFORM"))
            {lsrc=SetLightUNIFORM(pcmd,&scale_unit);}
          else if(pcmd->is_arg_value("type","GAUSSIAN"))
            {lsrc=SetLightGAUSSIAN(pcmd,&scale_unit);}
          else if(pcmd->is_arg_value("type","PULSE"))
            {lsrc=SetLightPULSE(pcmd,&scale_unit);}          
          else if(pcmd->is_arg_value("type","LSHELL"))
            {lsrc=SetLightSHELL(pcmd,&scale_unit);}
          else
          {
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" LSOURCE: no such source Type!\n";
               gss_log.record();
               return 1;
          }
        }
        catch (LSOURCE_CARD_ERROR err)
        {
            if(err == LSOURCE_UNKNOW_PARAMETER)
               gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" LSOURCE: unrecognized parameter(s)!\n";
            gss_log.record();
            return 1;
        }              
      }
      else
      {
           gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" LSOURCE: You should give LSOURCE type!\n";
           gss_log.record();
           return 1;
      
      }
    }
  }
  gss_log.string_buf()<<"done.\n";
  gss_log.record();
  return 0;
}

/* ----------------------------------------------------------------------------
 * optgen_update:  This function update optical generation rate
 */
int BSolver::optgen_update(PetscScalar currentT)
{

  double funcT = 0;
  if (lsrc) funcT = lsrc->currentft(currentT);
  for(int z=0;z<zone_num;z++)
  {
    if(zonedata[z]->material_type==Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *>(zonedata[z]);
      for(int i=0;i<zone[z].davcell.size();i++)
      {
        pzonedata->aux[i].RealOptG = pzonedata->aux[i].OptG*funcT;
      }
    }
  }
  return 0;
}
