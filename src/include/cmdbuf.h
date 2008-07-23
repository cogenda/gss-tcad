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

#ifndef _cmdparse_h
#define _cmdparse_h
#include <stdio.h>
#include <list>
#include <stack>
#include <map>
#include <string>
#include <string.h>
#include <stdarg.h>
using namespace std;

//structures work with lex/yacc

const int BOOLVAR  = 20;
const int INTERGER = 21;
const int DOUBLE   = 22;
const int CHARS    = 23;

// value can be double ,int, bool or string
typedef union
{
    bool   bval;
    int    ival;
    double dval;
    char   sval[32];
}Arg_value;

//this structure storage: parameter_name = arg_value;
typedef struct
{
  int       arg_type;
  string    arg_label;
  Arg_value arg_value;
}Arg;

typedef multimap<string, Arg> ARG_MAP;

//this structure storage: cmd_type  Keyword line number and  arg_map
//the sequence of parameters is not critical, multimap is very efficent.
class Cmd
{
public:
  string  KeyWord;
  int     lineno;
  ARG_MAP arg_map;
  ARG_MAP::iterator    parg;
public:
  int       get_current_lineno() {return lineno;};
  // search the arg list of current cmd
  bool      allowed_args(int argc,...);
  void      arg_search_begin() {parg=arg_map.begin();}
  void      goto_next_arg()    {parg++;}
  bool      arg_search_end()   {return parg==arg_map.end();}
  const char *get_cmd_keyword()    { return KeyWord.c_str();}
  const char *get_current_arg_label()  { return parg->second.arg_label.c_str();}
  Arg_value get_current_arg_value()  { return parg->second.arg_value;}
  // get parameter value of current cmd
  double    get_number(const char *name, int alias_num, ...);
  int       get_integer(const char *name, int alias_num, ...);
  bool      get_bool(const char *name, int alias_num, ...);
  char *    get_string(const char *name, int alias_num, ...);
  bool      is_arg_exist(const char *name);
  bool      is_arg_value(const char *name, const char * value);
};


extern int yyparse(list<Cmd> & cmd_list);
extern FILE* yyin;

// because the sequence of commands must be keep,
// I choose list as data structure.
class CmdBuf
{
public:
  list<Cmd> cmd_list;
  list<Cmd>::iterator    pcmd;
  stack<list<Cmd>::iterator> cmd_iterator_stack;
public:
  int       parse_file(char *cmd_file)
  {
    yyin=fopen(cmd_file,"r");
    if(yyin==NULL)
    { printf("Can't open input file, does it exist?\n"); return 1; }
    if(yyparse(cmd_list))
    { printf("Input parse interrupt.\n"); return 1; }
    fclose(yyin);
    return 0;
  }
  // cmd search
  void      cmd_search_begin()   {cmd_iterator_stack.push(pcmd);pcmd=cmd_list.begin();}
  void      goto_next_cmd()      {pcmd++;}
  bool      cmd_search_end()
  {
    if(pcmd==cmd_list.end())
    {
      pcmd=cmd_iterator_stack.top();
      cmd_iterator_stack.pop();
      return true;
    }
    else
      return false;
  }
  list<Cmd>::iterator    get_current_cmd()    {return pcmd;}
  bool      is_current_cmd(const char *cmd_name)    {return !strcmp(pcmd->KeyWord.c_str(),cmd_name);}
  void      delete_current_cmd()              {pcmd = cmd_list.erase(pcmd);}
  void      delete_cmd(const char *cmd_name)
  {
    for(pcmd=cmd_list.begin();pcmd!=cmd_list.end();pcmd++)
      if(!strcmp(cmd_name,pcmd->KeyWord.c_str())) pcmd = cmd_list.erase(pcmd);
  }
};

#endif
