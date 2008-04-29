/*****************************************************************************/
/*   	        8888888         88888888         88888888                    */
/*  	      8                8                8                            */
/* 	     8                 8                8                            */
/*  	     8                  88888888         88888888                    */
/* 	     8      8888                8                8                   */
/* 	      8       8                 8                8                   */
/* 	        888888         888888888        888888888                    */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "bsolver.h"
#include "log.h"
#include "typedef.h"

void SetDopingUniform(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale, vector<DopingFunc *> & doping_func_array)
{
  double ion  = 0;
  double xmax = pcmd->get_number("x.max",1,"x.right",0.0)*scale->s_micron;
  double xmin = pcmd->get_number("x.min",1,"x.left",0.0)*scale->s_micron;
  double ymax = pcmd->get_number("y.max",1,"y.top",0.0)*scale->s_micron;
  double ymin = pcmd->get_number("y.min",1,"y.bottom",0.0)*scale->s_micron;
  double peak = pcmd->get_number("n.peak",0,0.0)/pow(scale->s_centimeter,3);
  if(xmin-xmax>1e-8) throw PROFILE_XLOCATION;
  if(ymin-ymax>1e-8) throw PROFILE_YLOCATION;
  if(peak<0)    throw PROFILE_NEGTIVE_DOSE;
  // check doping ion
  if(pcmd->is_arg_exist("ion"))
  {
      if(pcmd->is_arg_value("ion","Donor"))
        ion = DONOR;
      else if(pcmd->is_arg_value("ion","Acceptor"))
        ion = ACCEPTOR;
      else
        throw PROFILE_ION_TYPE;
  }
  else
      throw PROFILE_ION_TYPE;

  if(!pcmd->allowed_args(11,"type","x.min","x.left","x.max","x.right","y.max",
                             "y.top","y.min","y.bottom","n.peak","ion"))
  {
	throw PROFILE_UNKNOW_PARAMETER;
  }

  DopingFunc * df = new UniformDopingFunc(xmin,xmax,ymax,ymin,ion,peak);
  doping_func_array.push_back(df);
  gss_log.string_buf()<<"\nDoping Uniform   : "<<peak*pow(scale->s_centimeter,3)<<" cm-3\n";
  gss_log.record();
}


void SetDopingGauss(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale, vector<DopingFunc *> & doping_func_array)
{
  double ion  = 0;
  if(pcmd->is_arg_exist("ion"))
  {
      if(pcmd->is_arg_value("ion","Donor"))
        ion = DONOR;
      else if(pcmd->is_arg_value("ion","Acceptor"))
        ion = ACCEPTOR;
      else
        throw PROFILE_ION_TYPE;
  }
  else
      throw PROFILE_ION_TYPE;
      
  double xmin = pcmd->get_number("x.min",1,"x.left",0.0)*scale->s_micron;
  double xmax = pcmd->get_number("x.max",1,"x.right",xmin)*scale->s_micron;
  double ymax = pcmd->get_number("y.max",1,"y.top",0.0)*scale->s_micron;
  double ymin = pcmd->get_number("y.min",1,"y.bottom",ymax)*scale->s_micron;
  double slice = 0.5*(xmax+xmin);
  if(xmin-xmax>1e-8) throw PROFILE_XLOCATION;
  if(ymin-ymax>1e-8) throw PROFILE_YLOCATION;
  
   
  double peak = pcmd->get_number("n.peak",0,0.0)/pow(scale->s_centimeter,3);

  double YCHAR,XCHAR,YJUNC,dop=0;
  if(pcmd->is_arg_exist("y.char"))
    	YCHAR = pcmd->get_number("y.char",0,0.25)*scale->s_micron;
  else if(pcmd->is_arg_exist("y.junction"))
  {
        //Junction is an absolute location.
	YJUNC = pcmd->get_number("y.junction",0,0)*scale->s_micron;
	//get conc. of background doping
	for(int i=0;i<doping_func_array.size();i++)
		dop += doping_func_array[i]->profile(slice,YJUNC);
	//Can we even find a junction?
	if(dop*ion>0) 	 throw PROFILE_ION_TYPE;
	if(peak<=0.0)    throw PROFILE_NEGTIVE_CONCENTR;
	//Now convert junction depth into char. length
	YCHAR = (ymin-YJUNC)/sqrt(log(fabs(peak/dop)));
  }
  if(pcmd->is_arg_exist("x.char"))
    XCHAR = pcmd->get_number("x.char",0,YCHAR/scale->s_micron)*scale->s_micron;
  else if(pcmd->is_arg_exist("xy.ratio"))
    XCHAR = YCHAR*pcmd->get_number("xy.ratio",0,1.0);
  else   
    XCHAR = YCHAR;
  if(XCHAR<=0 || YCHAR<=0 ) throw PROFILE_DOPING_CHAR_LENGTH;

  if(pcmd->is_arg_exist("dose"))
  {
    double dose = pcmd->get_number("dose",0,0.0)/pow(scale->s_centimeter,2);
    if(dose<0.0)    throw PROFILE_NEGTIVE_DOSE;
    peak=dose/(YCHAR*sqrt(3.14159265359)); 
  }
  
  if(!pcmd->allowed_args(16,"type","x.min","x.left","x.max","x.right","y.max",
                        "y.top","y.min","y.bottom","n.peak","dose","x.char","y.char","xy.ratio","y.junction","ion"))
  {
	throw PROFILE_UNKNOW_PARAMETER;
  }

  DopingFunc * df = new GaussDopingFunc(xmin,xmax,ymax,ymin,ion,peak,XCHAR,YCHAR);
  doping_func_array.push_back(df);
  gss_log.string_buf()<<"\nDoping Gauss     : "<<peak*pow(scale->s_centimeter,3)<<" cm-3\n";
  gss_log.record();
}

void SetDopingErf(list<Cmd>::iterator  pcmd, PhysicalUnitScale *scale, vector<DopingFunc *> & doping_func_array)
{
  double ion  = 0;
  double xmin = pcmd->get_number("x.min",1,"x.left",0.0)*scale->s_micron;
  double xmax = pcmd->get_number("x.max",1,"x.right",xmin)*scale->s_micron;
  double ymax = pcmd->get_number("y.max",1,"y.top",0.0)*scale->s_micron;
  double ymin = pcmd->get_number("y.min",1,"y.bottom",ymax)*scale->s_micron;
  double slice = 0.5*(xmax+xmin);
  double peak = pcmd->get_number("n.peak",0,0.0)/pow(scale->s_centimeter,3);
  if(xmin-xmax>1e-8) throw PROFILE_XLOCATION;
  if(ymin-ymax>1e-8) throw PROFILE_YLOCATION;
  if(peak<0)    throw PROFILE_NEGTIVE_DOSE;

  double YCHAR,XCHAR,YJUNC,dop=0;
  if(pcmd->is_arg_exist("ion"))
  {
      if(pcmd->is_arg_value("ion","Donor"))
        ion = DONOR;
      else if(pcmd->is_arg_value("ion","Acceptor"))
        ion = ACCEPTOR;
      else
        throw PROFILE_ION_TYPE;
  }
  else
      throw PROFILE_ION_TYPE;

  YCHAR = pcmd->get_number("y.char",0,0.25)*scale->s_micron;
  XCHAR = pcmd->get_number("x.char",0,YCHAR/scale->s_micron)*scale->s_micron;
  if(XCHAR<=0 || YCHAR<=0 ) throw PROFILE_DOPING_CHAR_LENGTH;

  if(!pcmd->allowed_args(13,"type","x.min","x.left","x.max","x.right","y.max",
                        "y.top","y.min","y.bottom","n.peak","x.char","y.char","ion"))
  {
	throw PROFILE_UNKNOW_PARAMETER;
  }

  DopingFunc * df = new ErfDopingFunc(xmin,xmax,ymax,ymin,ion,peak,XCHAR,YCHAR);
  doping_func_array.push_back(df);
  gss_log.string_buf()<<"\nDoping Errorfunc : "<<peak*pow(scale->s_centimeter,3)<<" cm-3\n";
  gss_log.record();
}

int BSolver::init_doping_func()
{
  double xmin=0,xmax=0,ymin=0,ymax=0;
  double PEAK=0,YCHAR=0.25,XCHAR=0.25;
  double ion=0;

  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
    if(pcmdbuf->is_current_cmd("PROFILE"))   // It's a doping card
    {
      list<Cmd>::iterator    pcmd = pcmdbuf->get_current_cmd();

      if(pcmd->is_arg_exist("type"))
      {
           try
	   {
	     if(pcmd->is_arg_value("type","Uniform"))
	   	SetDopingUniform(pcmd,&scale_unit,doping_func);
             else if(pcmd->is_arg_value("type","Gauss"))
		SetDopingGauss(pcmd,&scale_unit,doping_func);
             else if(pcmd->is_arg_value("type","ErrorFunc"))
                SetDopingErf(pcmd,&scale_unit,doping_func);
	     else
	     {
	       gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: no such doping Type!\n";
               gss_log.record();
	       return 1;
	     }
	   }
	   catch (PROFILE_CARD_ERROR err)
	   {
	       if(err==PROFILE_ION_TYPE)
	         gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: ION_TYPE error!\n";
	       else if(err==PROFILE_XLOCATION)
	       	 gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: X position of doping region error!\n";
	       else if(err==PROFILE_YLOCATION)
	       	 gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: Y position of doping region error!\n";
	       else if(err==PROFILE_NEGTIVE_DOSE)
	         gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: Must give a positive dose!\n";
	       else if(err==PROFILE_NEGTIVE_CONCENTR)
	         gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: Must give a positive concentration!\n";  
	       else if(err==PROFILE_DOPING_CHAR_LENGTH)
	         gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: char length must be positive!\n";
	       else if(err==PROFILE_UNKNOW_PARAMETER)
	         gss_log.string_buf()<<"line "<<pcmd->get_current_lineno()<<" PROFILE: unrecognized parameter(s)!\n";
	       gss_log.record();
	       
               return 1;
	   }
      }
    }
  return 0;
}

