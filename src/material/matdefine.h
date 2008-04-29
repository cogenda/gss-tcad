/*****************************************************************************/
/*                                                                           */
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
/*  GSS material database Version 0.4                                        */
/*  Last update: Feb 17, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#ifndef _mat_define_h_
#define _mat_define_h_

#include <string.h>

//define material index, each should have a positive value
const int Si     = 1;
const int GaAs   = 2;
const int Ge     = 3;
const int SiGe   = 4;
const int AlGaAs = 5;
const int InGaAs = 6;
const int InAs   = 7;
const int InP    = 8;
const int InN    = 9;
const int InSb   = 10;
const int HgCdTe = 11;
const int AlInAs = 20;
const int GaAsP  = 21;
const int InGaP  = 22;
const int InAsP  = 23;
const int SiC3C  = 24;
const int SiO2S  = 25;

const int SiO2   = 51;
const int Nitride= 52;
const int Air    = 53;

const int Al     = 70; //by zhangxih 06-09-28
const int Ag     = 71;
const int Au     = 72;
const int Cu     = 73;
const int PolySi = 74; //

//the material type GSS supproted now
const int Vacuum        = 150;
const int Semiconductor = 151;
const int Conductor     = 152;
const int Insulator     = 153;
const int PML           = 154;


//----------------------------------------------------------
// the following functions are used to judge material type
//----------------------------------------------------------
inline int IsSemiconductor(char *mat_name)
{
  //return true if mat_name specifies Semiconductor.
  if(!strcmp(mat_name,"Si"))     return Si;
  if(!strcmp(mat_name,"Silicon"))return Si;
  if(!strcmp(mat_name,"GaAs"))   return GaAs;
  if(!strcmp(mat_name,"Ge"))     return Ge;
  if(!strcmp(mat_name,"SiGe"))   return SiGe;
  if(!strcmp(mat_name,"AlGaAs")) return AlGaAs;
  if(!strcmp(mat_name,"InGaAs")) return InGaAs;
  if(!strcmp(mat_name,"InAs"))   return InAs;
  if(!strcmp(mat_name,"InSb"))   return InSb;
  if(!strcmp(mat_name,"InP"))    return InP;
  if(!strcmp(mat_name,"InN"))    return InN;
  if(!strcmp(mat_name,"HgCdTe")) return HgCdTe;
  if(!strcmp(mat_name,"AlInAs")) return AlInAs;
  if(!strcmp(mat_name,"GaAsP"))  return GaAsP;
  if(!strcmp(mat_name,"InGaP"))  return InGaP;
  if(!strcmp(mat_name,"InAsP"))  return InAsP;
  if(!strcmp(mat_name,"3C-SiC")) return SiC3C;
  if(!strcmp(mat_name,"S-SiO2")) return SiO2S;
  return 0;
}

inline int IsSingleCompSemiconductor(char *mat_name)
{
  //return true if mat_name specifies Single Compound Semiconductor.
  if(!strcmp(mat_name,"SiGe"))   return SiGe;
  if(!strcmp(mat_name,"AlGaAs")) return AlGaAs;
  if(!strcmp(mat_name,"InGaAs")) return InGaAs;
  if(!strcmp(mat_name,"HgCdTe")) return HgCdTe;
  if(!strcmp(mat_name,"InGa"))   return InGaAs;
  if(!strcmp(mat_name,"AlInAs")) return AlInAs;
  if(!strcmp(mat_name,"GaAsP"))  return GaAsP;
  if(!strcmp(mat_name,"InGaP"))  return InGaP;
  if(!strcmp(mat_name,"InAsP"))  return InAsP;
  return 0;
}

inline int IsInsulator(char *mat_name)
{
  //return true if mat_name specifies Insulator.
  if(!strcmp(mat_name,"Ox"))      return SiO2;
  if(!strcmp(mat_name,"SiO2"))    return SiO2;
  if(!strcmp(mat_name,"Nitride")) return Nitride;
  if(!strcmp(mat_name,"Air"))     return Air;
  return 0;
}

inline int IsElectrode(char *mat_name)
{
  //return true if mat_name specifies (Electrode) Conductor.
  if(!strcmp(mat_name,"Elec"))     return Conductor;
  if(!strcmp(mat_name,"Al"))       return Al;
  if(!strcmp(mat_name,"Ag"))       return Ag;
  if(!strcmp(mat_name,"Silver"))   return Ag;
  if(!strcmp(mat_name,"Au"))       return Au;
  if(!strcmp(mat_name,"Gold"))     return Au;
  if(!strcmp(mat_name,"Cu"))       return Cu;
  if(!strcmp(mat_name,"Copper"))   return Cu;
  if(!strcmp(mat_name,"PolySi"))   return PolySi;
  if(!strcmp(mat_name,"TiSi2"))    return Conductor;
  return 0;
}

inline int IsVacuum(char *mat_name)
{
  //return true if mat_name specifies Vacuum. Only used for EM FEM solver.
  if(!strcmp(mat_name,"Vacuum")) return Vacuum;
  return 0;
}

inline int IsPML(char *mat_name)
{
  //return true if mat_name specifies PML. Only used for EM FEM solver.
  if(!strcmp(mat_name,"PML"))     return PML;
  return 0;
}

//----------------------------------------------------------
// since some matrial has alias ( such as "Ox" and "SiO2" )
// format them to unique name
//----------------------------------------------------------
inline const char * FormatVacuumString(const char *mat_name)
{
  if(!strcmp(mat_name,"Vacuum")) return "Vacuum";
  return mat_name;
}

inline const char * FormatPMLString(const char *mat_name)
{
  if(!strcmp(mat_name,"PML"))    return "PML";
  return mat_name;
}

inline const char * FormatElectrodeString(const char *mat_name)
{
  if(!strcmp(mat_name,"Elec"))   return "Elec";
  if(!strcmp(mat_name,"PolySi")) return "PolySi";
  if(!strcmp(mat_name,"Al"))     return "Al";
  if(!strcmp(mat_name,"Ag"))     return "Ag";
  if(!strcmp(mat_name,"Silver")) return "Ag";
  if(!strcmp(mat_name,"Au"))     return "Au";
  if(!strcmp(mat_name,"Gold"))   return "Au";
  if(!strcmp(mat_name,"Cu"))     return "Cu";
  if(!strcmp(mat_name,"Copper")) return "Cu";
  if(!strcmp(mat_name,"TiSi2"))  return "TiSi2";
  return mat_name;
}

inline const char * FormatInsulatorString(const char *mat_name)
{
  if(!strcmp(mat_name,"Ox"))      return "SiO2";
  if(!strcmp(mat_name,"SiO2"))    return "SiO2";
  if(!strcmp(mat_name,"Nitride")) return "Nitride";
  if(!strcmp(mat_name,"Si3N4"))   return "Nitride";
  if(!strcmp(mat_name,"Air"))     return "Air";
  return mat_name;
}

inline const char * FormatSemiconductorString(const char *mat_name)
{
  if(!strcmp(mat_name,"Si"))      return "Si";
  if(!strcmp(mat_name,"Silicon")) return "Si";
  if(!strcmp(mat_name,"GaAs"))    return "GaAs";
  if(!strcmp(mat_name,"Ge"))      return "Ge";
  if(!strcmp(mat_name,"SiGe"))    return "SiGe";
  if(!strcmp(mat_name,"AlGaAs"))  return "AlGaAs";
  if(!strcmp(mat_name,"InGaAs"))  return "InGaAs";
  if(!strcmp(mat_name,"InAs"))    return "InAs";
  if(!strcmp(mat_name,"InSb"))    return "InSb";
  if(!strcmp(mat_name,"InP"))     return "InP";
  if(!strcmp(mat_name,"InN"))     return "InN";
  if(!strcmp(mat_name,"HgCdTe"))  return "HgCdTe";
  if(!strcmp(mat_name,"AlInAs"))  return "AlInAs";
  if(!strcmp(mat_name,"GaAsP"))   return "GaAsP";
  if(!strcmp(mat_name,"InGaP"))   return "InGaP";
  if(!strcmp(mat_name,"InAsP"))   return "InAsP";
  if(!strcmp(mat_name,"3C-SiC"))  return "SiC3C";
  if(!strcmp(mat_name,"S-SiO2"))  return "SiO2S";
  return mat_name;
}

#endif
