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
/*  Last update: May 14, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include <string.h>
#include "cmdbuf.h"


bool  Cmd::allowed_args(int argc,...)
{
  va_list arg;
  va_start(arg,argc);

  for(parg=arg_map.begin();parg!=arg_map.end();parg++)
  {
    int flag = 0;
    va_start(arg,argc);
    for(int i=0;i<argc;i++)
      if( parg->second.arg_label==va_arg(arg,char *)) flag = 1;
    if(!flag)
    {
        va_end(arg);
        return false;
    }
    va_end(arg);
  }
  return true;
}

double Cmd::get_number(const char *name, int alias_num, ...)
{
  string arg_name=name;
  string alias[20];
  va_list arg;
  va_start(arg,alias_num);
  for(int i=0;i<alias_num;i++)
       alias[i] = va_arg(arg,char *);
  double value = va_arg(arg,double);
  va_end(arg);

  if((parg=arg_map.find(arg_name))!=arg_map.end())
        return parg->second.arg_value.dval;

  for(int i=0;i<alias_num;i++)
      if((parg=arg_map.find(alias[i]))!=arg_map.end())
        return parg->second.arg_value.dval;

  return value;
}


int Cmd::get_integer(const char *name, int alias_num, ...)
{
  string arg_name=name;
  string alias[20];
  va_list arg;
  va_start(arg,alias_num);
  for(int i=0;i<alias_num;i++)
       alias[i] = va_arg(arg,char *);
  int value = va_arg(arg,int);
  va_end(arg);

  if((parg=arg_map.find(arg_name))!=arg_map.end())
        return parg->second.arg_value.ival;

  for(int i=0;i<alias_num;i++)
      if((parg=arg_map.find(alias[i]))!=arg_map.end())
        return parg->second.arg_value.ival;

  return value;
}


bool Cmd::get_bool(const char *name, int alias_num, ...)
{
  string arg_name=name;
  string alias[20];
  va_list arg;
  va_start(arg,alias_num);
  for(int i=0;i<alias_num;i++)
       alias[i] = va_arg(arg,char *);
  bool value = va_arg(arg,int);
  va_end(arg);

  if((parg=arg_map.find(arg_name))!=arg_map.end())
        return parg->second.arg_value.bval;

  for(int i=0;i<alias_num;i++)
      if((parg=arg_map.find(alias[i]))!=arg_map.end())
        return parg->second.arg_value.bval;

  return value;
}


char * Cmd::get_string(const char *name, int alias_num, ...)
{
  string arg_name=name;
  string alias[20];
  va_list arg;
  va_start(arg,alias_num);
  for(int i=0;i<alias_num;i++)
       alias[i] = va_arg(arg,char *);
  char * value = va_arg(arg,char *);
  va_end(arg);

  if((parg=arg_map.find(arg_name))!=arg_map.end())
        return parg->second.arg_value.sval;

  for(int i=0;i<alias_num;i++)
      if((parg=arg_map.find(alias[i]))!=arg_map.end())
        return parg->second.arg_value.sval;

  return value;
}


bool Cmd::is_arg_exist(const char *name)
{
  string arg_name=name;
  if((parg=arg_map.find(arg_name))!=arg_map.end())
        return true;
  return false;
}


bool Cmd::is_arg_value(const char *name, const char * value)
{
  string arg_name=name;
  if((parg=arg_map.find(arg_name))!=arg_map.end())
    if(!strcmp(parg->second.arg_value.sval,value))
        return true;
    else
        return false;
  return false;
}




