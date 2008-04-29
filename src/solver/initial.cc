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
/*  Last update: Apr 21, 2006                                                */
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


//simple error hanfle function
inline int read_config_err_handle(char *s)
{
  gss_log.string_buf()<<s;
  gss_log.record();
  return 1;
}

/* ----------------------------------------------------------------------------
 * BSolver::static_config:  This function get global config data from CmdParseBuf
 */
int BSolver::static_config(CmdBuf *cmdbuf)
{
  pcmdbuf = cmdbuf;
  //set default value
  mesh_scale = 1e6;
  Carrier = PN_Both;
  LatticeTemp = 300;
  DeviceDepth = 1.0;

  //search the CmdParseBuf. setup global variables
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();)
  {
    if(pcmdbuf->is_current_cmd("SET"))   //it is a set structure
    {
      list<Cmd>::iterator    pcmd = pcmdbuf->get_current_cmd();
      if(pcmd->is_arg_exist("latticetemp"))
        LatticeTemp = pcmd->get_number("latticetemp",0,300.0);
      if(pcmd->is_arg_exist("z.width"))
        DeviceDepth = pcmd->get_number("z.width",0,1.0);
      if(pcmd->is_arg_exist("dopingscale"))
        mesh_scale  = pow(pcmd->get_number("dopingscale",0,1e18),1.0/3);
      if(pcmd->is_arg_exist("carrier"))
      {
          if(pcmd->is_arg_value("carrier","pn"))        Carrier=PN_Both;
          else if(pcmd->is_arg_value("carrier","p"))    Carrier=P_Type;
          else if(pcmd->is_arg_value("carrier","n"))    Carrier=N_Type;
          else {sprintf(log_buf,"line %d Set: Unrecognized carrier type!\n",pcmd->lineno);GSS_LOG();return 1;}
      }
      pcmdbuf->delete_current_cmd();
    }
    else
      pcmdbuf->goto_next_cmd();
  }

  //scale variables
  scale_unit.SetPhysicalUnitScale(mesh_scale);
  DeviceDepth = DeviceDepth*scale_unit.s_micron;
  LatticeTemp = LatticeTemp*scale_unit.s_kelvin;
  set_mesh_scale(scale_unit.s_centimeter);
  show_config_info();

  //set voltage and current sources
  if(init_esource(pcmdbuf->cmd_list)) return 1;
  //set PMI
  if(init_PMI_Semiconductor(pcmdbuf,PMIS_define)) return 1;
  //set doping function
  if(init_doping_func()) return 1;

  //set initial grid
  if(init_grid(pcmdbuf->cmd_list)) return 1;

 

  gss_log.string_buf()<<"\n";
  gss_log.record();
  return 0;
}


/* ----------------------------------------------------------------------------
 * BSolver::config_info:  This function print global config data to screen
 * mainly for debug
 */
void BSolver::show_config_info()
{

  gss_log.string_buf()<<"\n--------------------------------------------------\n";
  gss_log.record();
  switch(Carrier)
  {
  case N_Type:  gss_log.string_buf()<<"carrier            : n  only\n";gss_log.record();break;
  case P_Type:  gss_log.string_buf()<<"carrier            : p  only\n";gss_log.record();break;
  case PN_Both: gss_log.string_buf()<<"carrier            : pn both\n";gss_log.record();break;
  }

  gss_log.string_buf()<<"device Z-width     : "<<DeviceDepth/scale_unit.s_micron<< "um\n";
  gss_log.record();

  gss_log.string_buf()<<"lattice temprature : "<<LatticeTemp/scale_unit.s_kelvin<< "K\n";
  gss_log.record();
}

/* ----------------------------------------------------------------------------
 * BSolver::clear_data:  This function delete solution data structure
 */
void BSolver::clear_data()
{
  if(zonedata.size()==0) return;
  for(int i=0;i<zonedata.size();i++)
    delete zonedata[i];
  zonedata.clear();  
}

BSolver::~BSolver()
{
  clear_data();
  for(int i=0;i<vsrc.size();i++)
    delete vsrc[i];
  for(int i=0;i<isrc.size();i++)
    delete isrc[i];
  for(int i=0;i<doping_func.size();i++)
    delete doping_func[i];
}


