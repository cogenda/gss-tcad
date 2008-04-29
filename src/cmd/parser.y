%{

#include <stdio.h>
#include <string.h>
#include "cmdbuf.h"

extern int yylineno;
extern int yylex();
int yyerror(list<Cmd> & cmd_list, const char *s);

/* define global variable needed by yyparse here*/
Arg arg;
Cmd cmd;

void ClearCmd(Cmd &cmd);


/*#define DEBUG*/

%}

%start inputfile

%union{
    bool   bval;
    int    ival;
    double dval;
    char   sval[32];
}

%parse-param { list<Cmd> & cmd_list }

%token <sval> ASSIGN KEYWORD  BOOL_PARAMETER REAL_PARAMETER INT_PARAMETER STRING_PARAMETER
%token <ival> INT_NUMBER
%token <dval> REAL_NUMBER
%token <sval> STRING
%token <ival> BLANK COMMENT BAD_WORD
%token <bval> BOOL

%type  <dval> real_expr
%type  <ival> int_expr
%left '+' '-' 
%left '*' '/'
%%

inputfile
    :   command
        {}
    |   inputfile command
        {}
    ;

command
    :   ASSIGN STRING '=' real_expr
{
#ifdef DEBUG
	printf("line %d: --YACC assign --\n",yylineno-1);
#endif
        arg.arg_label=$2;
        arg.arg_type=DOUBLE;
	arg.arg_value.dval=$4;
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
	cmd.KeyWord=$1;
	cmd.lineno = yylineno-1;
	cmd_list.push_back(cmd);
	ClearCmd(cmd);
}    	
    |   KEYWORD
{
#ifdef DEBUG
	printf("line %d: --YACC command %s--\n",yylineno-1,$1);
#endif
	cmd.KeyWord=$1;
	cmd.lineno = yylineno-1;
	cmd_list.push_back(cmd);
	ClearCmd(cmd);
}
    |   KEYWORD parameter
{
#ifdef DEBUG
	printf("line %d: --YACC command %s--\n",yylineno-1,$1);
#endif
	cmd.KeyWord=$1;
	cmd.lineno = yylineno-1;
	cmd_list.push_back(cmd);
	ClearCmd(cmd);
}
     |  COMMENT
     |  BLANK
     ;



real_expr: REAL_NUMBER     
	 | INT_NUMBER  
	   { $$=$1;}
	 | STRING         {
	     int flag=0;	
	     double assign_value;	
	     list<Cmd>::iterator    pcmd;
	     for(pcmd=cmd_list.begin();pcmd!=cmd_list.end();pcmd++)
	      if(!strcmp(pcmd->KeyWord.c_str(),"ASSIGN"))
	      {
	        ARG_MAP::iterator    parg;
	        if((parg=pcmd->arg_map.find($1))!=pcmd->arg_map.end())
	        {
	        	assign_value=parg->second.arg_value.dval;
	        	flag=1;
                }
	      }
	     if(flag) 	
	       $$=assign_value;
	     else
	       return yyerror(cmd_list, $1); 
	   }	  
	 | '(' real_expr ')'
	   { $$=$2;}
	 | real_expr '+' real_expr
	   {$$=$1+$3;}
	 | real_expr '-' real_expr 
	   {$$=$1-$3;}
	 | real_expr '*' real_expr
	   {$$=$1*$3;}
	 | real_expr '/' real_expr 
	   {$$=$1/$3;}  
	 | '-' real_expr  %prec '*'
	   {$$=-$2;} 
	 ;  

int_expr : INT_NUMBER  
	   { $$=$1;}
	 | '(' int_expr ')'
	   { $$=$2;}
	 | int_expr '+' int_expr
	   {$$=$1+$3;}
	 | int_expr '-' int_expr 
	   {$$=$1-$3;}
	 | int_expr '*' int_expr
	   {$$=$1*$3;}
	 | int_expr '/' int_expr 
	   {$$=$1/$3;}  
	 | '-' int_expr  %prec '*'
	   {$$=-$2;} 
	 ;  
	 
parameter
     : parameter REAL_PARAMETER '=' real_expr
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN REAL--\n");
#endif
	arg.arg_label=$2;
        arg.arg_type=DOUBLE;
	arg.arg_value.dval=$4;
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | REAL_PARAMETER '=' real_expr
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN REAL--\n");
#endif
	arg.arg_label=$1;
        arg.arg_type=DOUBLE;
	arg.arg_value.dval=$3;
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | parameter INT_PARAMETER '=' int_expr
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN INT--\n");
#endif
	arg.arg_label=$2;
        arg.arg_type=INTERGER;
	arg.arg_value.ival=$4;
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | INT_PARAMETER '=' int_expr
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN INT--\n");
#endif
	arg.arg_label=$1;
        arg.arg_type=INTERGER;
	arg.arg_value.ival=$3;
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | parameter STRING_PARAMETER '=' STRING
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN STRING--\n");
#endif
	arg.arg_label=$2;
        arg.arg_type=CHARS;
	strcpy(arg.arg_value.sval,$4);
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | STRING_PARAMETER '=' STRING
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN STRING--\n");
#endif
	arg.arg_label=$1;
        arg.arg_type=CHARS;
	strcpy(arg.arg_value.sval,$3);
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | parameter STRING_PARAMETER
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN EMPTY STRING--\n");
#endif
	// in this case, assign an empty string
	arg.arg_label=$2;
        arg.arg_type=CHARS;
	strcpy(arg.arg_value.sval,"");
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | STRING_PARAMETER
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN EMPTY STRING--\n");
#endif
	// in this case, assign an empty string
	arg.arg_label=$1;
        arg.arg_type=CHARS;
	strcpy(arg.arg_value.sval,"");
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | parameter BOOL_PARAMETER '=' BOOL
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN BOOL--\n");
#endif
	arg.arg_label=$2;
        arg.arg_type=BOOLVAR;
	arg.arg_value.bval=$4;
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     | BOOL_PARAMETER '=' BOOL
{
#ifdef DEBUG
	printf("--PARAMETER ASSIGN BOOL--\n");
#endif
	arg.arg_label=$1;
        arg.arg_type=BOOLVAR;
	arg.arg_value.bval=$3;
	cmd.arg_map.insert(pair<const string, Arg>(arg.arg_label,arg));
}
     ;
     
     

%%

int yyerror(list<Cmd> & cmd_list, const char *s)
{
   printf("\nYACC report: line %d unrecognized word(s), %s.\n",yylineno,s);
   return 1;
}

void ClearCmd(Cmd &cmd)
{
   cmd.KeyWord.clear();
   cmd.arg_map.clear();
}
