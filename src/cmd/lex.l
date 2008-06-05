%{
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <ctype.h>
#include "parser.h"

/*#define DEBUG*/

int initial = 1;

%}

D       [0-9]
E       [EeDd][+-]?{D}+
SIGN    ([+-]?)

%s COMMAND PARAMETERS

%%

        if(initial==1)
        {
                initial = 0;
                BEGIN(COMMAND);
        }

^\#[^\n]*                         {
/* comment line */
#ifdef DEBUG
        printf("COMMENT LINE:");ECHO;
#endif
        return COMMENT;
}

\#[^\n]*                          {
/* inline comment */
#ifdef DEBUG
        printf("COMMENT:");ECHO;
#endif
}


\\[ \t\r]*\n                      {
/* continue line */
#ifdef DEBUG
        printf("CONTINUE LINE:");ECHO;
#endif
}


\n                                {
/* new command */
BEGIN COMMAND;
#ifdef DEBUG
        ECHO;
#endif
        return BLANK;
}

[ \t\r]*                          {
/* use to drop blankspace/tab/cr */
#ifdef DEBUG
        ECHO;
#endif
}


<COMMAND>ASSIGN        {
BEGIN  PARAMETERS;
#ifdef DEBUG
        printf("ASSIGN");
#endif
        strncpy(yylval.sval, yytext,31);
        for(unsigned int c=0;c<strlen(yylval.sval);c++)
        {
          if(islower(yylval.sval[c]))
            yylval.sval[c]=toupper(yylval.sval[c]);
        }
        return ASSIGN;
}


<COMMAND>SET           |
<COMMAND>VSOURCE       |
<COMMAND>ISOURCE       |
<COMMAND>LSOURCE       |
<COMMAND>BOUNDARY      |
<COMMAND>CONTACT       |
<COMMAND>XMESH         |
<COMMAND>YMESH         |
<COMMAND>MESH          |
<COMMAND>ELIMINATE     |
<COMMAND>SPREAD        |
<COMMAND>REGION        |
<COMMAND>SEGMENT       |
<COMMAND>ELECTRODE     |
<COMMAND>PROFILE       |
<COMMAND>METHOD        |
<COMMAND>PMIS          |
<COMMAND>SOLVE         |
<COMMAND>EXPORT        |
<COMMAND>IMPORT        |
<COMMAND>ATTACH        |
<COMMAND>REFINE        |
<COMMAND>PLOT          |
<COMMAND>PLOTMESH      |
<COMMAND>PROBE         |
<COMMAND>PLOTVTK       |
<COMMAND>PHOTOGEN      |
<COMMAND>END                               {
BEGIN  PARAMETERS;
#ifdef DEBUG
        printf("KEYWORD:");ECHO;
#endif
        strncpy(yylval.sval, yytext,31);
        for(unsigned int c=0;c<strlen(yylval.sval);c++)
        {
          if(islower(yylval.sval[c]))
            yylval.sval[c]=toupper(yylval.sval[c]);
        }
        return KEYWORD;
}


<PARAMETERS>X"."ORG         |//copyed by zhangxih ----06-12-20
<PARAMETERS>Y"."ORG         |
<PARAMETERS>ANGLE           |
<PARAMETERS>WAVELEN         |
<PARAMETERS>INTENSITY       |
<PARAMETERS>WAVE"."STA      |
<PARAMETERS>WAVE"."END      |
<PARAMETERS>WAVE"."NUM      |
<PARAMETERS>RAY"."WIDTH     |
<PARAMETERS>QUAN"."EFF      |
<PARAMETERS>WTE             |
<PARAMETERS>WTM             |
<PARAMETERS>phase"."diff    |
<PARAMETERS>power           |
<PARAMETERS>powerlo         |
<PARAMETERS>powerhi         |
<PARAMETERS>Tpeak           |
<PARAMETERS>Tradius         |
<PARAMETERS>TC              |
<PARAMETERS>Width           |
<PARAMETERS>Depth           |
<PARAMETERS>H1              |
<PARAMETERS>H2              |
<PARAMETERS>X               |
<PARAMETERS>Y               |
<PARAMETERS>X"."Min         |
<PARAMETERS>X"."Left        |
<PARAMETERS>X"."Max         |
<PARAMETERS>X"."Right       |
<PARAMETERS>Y"."Min         |
<PARAMETERS>Y"."Bottom      |
<PARAMETERS>Y"."Max         |
<PARAMETERS>Y"."Top         |
<PARAMETERS>CentreX         |
<PARAMETERS>CentreY         |
<PARAMETERS>MajorRadii      |
<PARAMETERS>MinorRadii      |
<PARAMETERS>Theta           |
<PARAMETERS>Ratio           |
<PARAMETERS>Encroach        |
<PARAMETERS>Y.Lower         |
<PARAMETERS>VOL.RAT         |
<PARAMETERS>Grading         |
<PARAMETERS>Middle          |
<PARAMETERS>Y.Middle        |
<PARAMETERS>GR1             |
<PARAMETERS>GR2             |
<PARAMETERS>N"."Peak        |
<PARAMETERS>Dose            |
<PARAMETERS>X"."CHAR        |
<PARAMETERS>Y"."CHAR        |
<PARAMETERS>XY"."RATIO      |
<PARAMETERS>Y"."JUNCTION    |
<PARAMETERS>X"."MOLE        |
<PARAMETERS>Mole"."Slope    |
<PARAMETERS>Mole"."End      |
<PARAMETERS>Z.Width         |
<PARAMETERS>LatticeTemp     |
<PARAMETERS>DopingScale     |
<PARAMETERS>Tdelay          |
<PARAMETERS>Vconst          |
<PARAMETERS>Vamp            |
<PARAMETERS>Iconst          |
<PARAMETERS>Iamp            |
<PARAMETERS>Freq            |
<PARAMETERS>Alpha           |
<PARAMETERS>TRC             |
<PARAMETERS>TFD             |
<PARAMETERS>TFC             |
<PARAMETERS>Vlo             |
<PARAMETERS>Vhi             |
<PARAMETERS>Ilo             |
<PARAMETERS>Ihi             |
<PARAMETERS>Tr              |
<PARAMETERS>Tf              |
<PARAMETERS>Pw              |
<PARAMETERS>Pr              |
<PARAMETERS>Res             |
<PARAMETERS>Cap             |
<PARAMETERS>Ind             |
<PARAMETERS>Eps             |
<PARAMETERS>Thickness       |
<PARAMETERS>QF              |
<PARAMETERS>Heat"."Transfer |
<PARAMETERS>EXT"."TEMP      |
<PARAMETERS>WorkFunction    |
<PARAMETERS>VBarrier        |
<PARAMETERS>Dispersion      |
<PARAMETERS>DivisionRatio   |
<PARAMETERS>QNFactor        |
<PARAMETERS>QPFactor        |
<PARAMETERS>relative"."tol  |
<PARAMETERS>toler"."relax   |
<PARAMETERS>possion"."tol   |
<PARAMETERS>elec"."continuty"."tol |
<PARAMETERS>hole"."continuty"."tol |
<PARAMETERS>latt"."temp"."tol      |
<PARAMETERS>elec"."energy"."tol    | 
<PARAMETERS>hole"."energy"."tol    | 
<PARAMETERS>elec"."quantum"."tol   | 
<PARAMETERS>hole"."quantum"."tol   | 
<PARAMETERS>electrode"."tol        | 
<PARAMETERS>VStart                 |
<PARAMETERS>VStep                  |
<PARAMETERS>VStepMax               |
<PARAMETERS>VStop                  |
<PARAMETERS>IStart                 |
<PARAMETERS>IStep                  |
<PARAMETERS>IStop                  |
<PARAMETERS>TStart                 |
<PARAMETERS>TStep                  |
<PARAMETERS>TStop                  |
<PARAMETERS>VAC                    |
<PARAMETERS>FStart                 |
<PARAMETERS>FMultiple              |
<PARAMETERS>FStop                  |
<PARAMETERS>AutoSave               |
<PARAMETERS>AzAngle                |
<PARAMETERS>ElAngle                               {

#ifdef DEBUG
        printf("PARAMETER:");ECHO;
#endif
        strncpy(yylval.sval, yytext,31);
        for(unsigned int c=0;c<strlen(yylval.sval);c++)
        {
          if(isupper(yylval.sval[c]))
            yylval.sval[c]=tolower(yylval.sval[c]);
        }
        return REAL_PARAMETER;
}


<PARAMETERS>IX            |
<PARAMETERS>IX            |
<PARAMETERS>IX"."Min      |
<PARAMETERS>IX"."Left     |
<PARAMETERS>IX"."Max      |
<PARAMETERS>IX"."Right    |
<PARAMETERS>IY"."Min      |
<PARAMETERS>IY"."Bottom   |
<PARAMETERS>IY"."Max      |
<PARAMETERS>IY"."Top      |
<PARAMETERS>Upper         |
<PARAMETERS>Lower         |
<PARAMETERS>N"."Spaces    |
<PARAMETERS>Division      |
<PARAMETERS>MaxIteration  |
<PARAMETERS>ServerPort                        {

#ifdef DEBUG
        printf("PARAMETER:");ECHO;
#endif
        strncpy(yylval.sval, yytext,31);
        for(unsigned int c=0;c<strlen(yylval.sval);c++)
        {
          if(isupper(yylval.sval[c]))
            yylval.sval[c]=tolower(yylval.sval[c]);
        }
        return INT_PARAMETER;
}


<PARAMETERS>Direction        |
<PARAMETERS>Label            |
<PARAMETERS>Material         |
<PARAMETERS>Location         |
<PARAMETERS>Shape            |
<PARAMETERS>Mole.Grad        |
<PARAMETERS>ModelFile        |
<PARAMETERS>Triangle         |
<PARAMETERS>Ion              |
<PARAMETERS>Carrier          |
<PARAMETERS>Type             |
<PARAMETERS>ID               |
<PARAMETERS>Variable         |
<PARAMETERS>Measure          |
<PARAMETERS>Region           |
<PARAMETERS>Mobility         |
<PARAMETERS>II.Model         |
<PARAMETERS>Optical.Model    |
<PARAMETERS>Scheme           |
<PARAMETERS>NS               |
<PARAMETERS>LS               |
<PARAMETERS>Damping          |
<PARAMETERS>II"."Type        |
<PARAMETERS>IVRecord         |
<PARAMETERS>IVFile           |
<PARAMETERS>CoreFile         |
<PARAMETERS>AscFile          |
<PARAMETERS>VTKFile          |
<PARAMETERS>Segment          |
<PARAMETERS>ProbeFile        |
<PARAMETERS>Electrode        |
<PARAMETERS>ConnectTo        |
<PARAMETERS>VScan            |
<PARAMETERS>IScan            |
<PARAMETERS>ACScan           |
<PARAMETERS>DLL              |
<PARAMETERS>Func             |
<PARAMETERS>ODE"."Formula    |
<PARAMETERS>VApp             |
<PARAMETERS>IApp             |
<PARAMETERS>Resolution       |
<PARAMETERS>PS.Out           |
<PARAMETERS>TIFF.Out         |
<PARAMETERS>Style                             {

#ifdef DEBUG
        printf("PARAMETER:");ECHO;
#endif
        strncpy(yylval.sval, yytext,31);
        for(unsigned int c=0;c<strlen(yylval.sval);c++)
        {
          if(isupper(yylval.sval[c]))
            yylval.sval[c]=tolower(yylval.sval[c]);
        }
        return STRING_PARAMETER;
}


<PARAMETERS>Fix.Lower         |
<PARAMETERS>HighFieldMobility | 
<PARAMETERS>ImpactIonization  |
<PARAMETERS>BandBandTunneling |
<PARAMETERS>Fermi             |
<PARAMETERS>PseudoTime        |
<PARAMETERS>Append            |
<PARAMETERS>EJModel           |
<PARAMETERS>Projection                        {

#ifdef DEBUG
        printf("PARAMETER:");ECHO;
#endif
        strncpy(yylval.sval, yytext,31);
        for(unsigned int c=0;c<strlen(yylval.sval);c++)
        {
          if(isupper(yylval.sval[c]))
            yylval.sval[c]=tolower(yylval.sval[c]);
        }
        return BOOL_PARAMETER;
}


<PARAMETERS>\=                                {
#ifdef DEBUG
        ECHO;
#endif
        return '=';
}

<PARAMETERS>\+                                {
#ifdef DEBUG
        ECHO;
#endif
        return '+';
}

<PARAMETERS>\-                                {
#ifdef DEBUG
        ECHO;
#endif
        return '-';
}

<PARAMETERS>\*                                {
#ifdef DEBUG
        ECHO;
#endif
        return '*';
}

<PARAMETERS>\/                                {
#ifdef DEBUG
        ECHO;
#endif
        return '/';
}

<PARAMETERS>\(                                {
#ifdef DEBUG
        ECHO;
#endif
        return '(';
}

<PARAMETERS>\)                                {
#ifdef DEBUG
        ECHO;
#endif
        return ')';
}

<PARAMETERS>{D}+                        {
        sscanf(yytext, "%d", &yylval.ival);
#ifdef DEBUG
        printf("int_value: %d",yylval.ival);
#endif
        return INT_NUMBER;
}


<PARAMETERS>{D}+"."{D}*({E})?     |
<PARAMETERS>{D}*"."{D}+({E})?     |
<PARAMETERS>{D}+({E})                       {
        sscanf(yytext, "%lf", &yylval.dval);
#ifdef DEBUG
        printf("double_value: %lf",yylval.dval);
#endif
        return REAL_NUMBER;
}



<PARAMETERS>TRUE         |
<PARAMETERS>True         |
<PARAMETERS>true         |
<PARAMETERS>ON           |
<PARAMETERS>On           |
<PARAMETERS>on                                {
        yylval.bval = true;
#ifdef DEBUG
        printf("BOOL True:%s",yytext);
#endif
        return BOOL;
}



<PARAMETERS>FALSE         |
<PARAMETERS>False         |
<PARAMETERS>false         |
<PARAMETERS>OFF           |
<PARAMETERS>Off           |
<PARAMETERS>off                                {
        yylval.bval = false;
#ifdef DEBUG
        printf("BOOL False:%s",yytext);
#endif
        return BOOL;
}


<PARAMETERS>[a-zA-Z0-9\_][a-zA-Z0-9\.\_]*     {
#ifdef DEBUG
        printf("STRING:%s",yytext);
#endif
        strncpy(yylval.sval,yytext,31);
        return STRING;
}



<PARAMETERS>\"[a-zA-Z0-9\.\_\ \-]*\"          {
        yytext[yyleng-1]='\0';
#ifdef DEBUG
        printf("STRING:%s",&yytext[1]);
#endif
        strncpy(yylval.sval, &yytext[1],31);
        return STRING;
}

<COMMAND>[a-zA-Z0-9[:punct:]]+        {
        printf("Bad String: %s",yytext);
        strncpy(yylval.sval,yytext,31);
        return BAD_WORD;
}


%%
int yywrap()
{
#ifdef DEBUG
        printf("\nline %d\n",yylineno);
#endif
        return 1;

}
