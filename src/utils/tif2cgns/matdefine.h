#ifndef _matdefine_h_
#define _matdefine_h_
//define material_type
const int Si     = 1;
const int GaAs   = 2;
const int Ge     = 3;
const int SiGe   = 4;
const int AlGaAs = 5;
const int InGaAs = 6;
const int AlInAs = 7;
const int GaAsP  = 8;
const int InGaP  = 9;
const int InAsP  = 10;
const int InP    = 11;
const int InAs   = 12;
const int PolySi = 13;
const int SiC    = 14;
const int ZnSe   = 15;
const int ZnTe   = 16;
const int HgCdTe = 17;

const int SiO2   = -1;
const int Nitride= -2;
const int Al     = -10;

const int Semiconductor = 151;
const int Conductor     = 152;
const int Insulator     = 153;


int   IsInsulator(char *);
int   IsSemiconductor(char *);
int   IsElectrode(char *);

inline int IsInsulator(char *mat_name)
{
  if(!strcmp(mat_name,"Ox")) return SiO2;
  if(!strcmp(mat_name,"SiO2")) return SiO2;
  if(!strcmp(mat_name,"Nit")) return Nitride;
  return 0;
}

inline int IsSemiconductor(char *mat_name)
{
  if(!strcmp(mat_name,"Si"))     return Si;
  if(!strcmp(mat_name,"Silicon"))return Si;
  if(!strcmp(mat_name,"GaAs"))   return GaAs;
  if(!strcmp(mat_name,"Ge"))     return Ge;
  if(!strcmp(mat_name,"SiGe"))   return SiGe;
  if(!strcmp(mat_name,"AlGaAs")) return AlGaAs;
  if(!strcmp(mat_name,"InGaAs")) return InGaAs;
  if(!strcmp(mat_name,"InGa"))   return InGaAs;
  if(!strcmp(mat_name,"AlInAs")) return AlInAs;
  if(!strcmp(mat_name,"GaAsP"))  return GaAsP;
  if(!strcmp(mat_name,"InGaP"))  return InGaP;
  if(!strcmp(mat_name,"InAsP"))  return InAsP;
  if(!strcmp(mat_name,"InP"))    return InP;
  if(!strcmp(mat_name,"InAs"))   return InAs;
  if(!strcmp(mat_name,"SiC"))    return SiC;
  if(!strcmp(mat_name,"ZnSe"))   return ZnSe;
  if(!strcmp(mat_name,"ZnTe"))   return ZnTe;
  if(!strcmp(mat_name,"HgCdTe")) return HgCdTe;
 
  return 0;
}

inline int IsSingleCompSemiconductor(char *mat_name)
{
  if(!strcmp(mat_name,"SiGe"))   return SiGe;
  if(!strcmp(mat_name,"AlGaAs")) return AlGaAs;
  if(!strcmp(mat_name,"InGaAs")) return InGaAs;
  if(!strcmp(mat_name,"InGa"))   return InGaAs;
  if(!strcmp(mat_name,"AlInAs")) return AlInAs;
  if(!strcmp(mat_name,"GaAsP"))  return GaAsP;
  if(!strcmp(mat_name,"InGaP"))  return InGaP;
  if(!strcmp(mat_name,"InAsP"))  return InAsP;
  return 0;
}

inline int IsElectrode(char *mat_name)
{
  if(!strcmp(mat_name,"Elec"))     return Conductor;
  if(!strcmp(mat_name,"Poly"))     return PolySi;
  if(!strcmp(mat_name,"PolySi"))   return PolySi;
  return 0;
}

inline const char * FormatMaterialString(const char *mat_name)
{
  if(!strcmp(mat_name,"Elec"))    return "Elec";
  if(!strcmp(mat_name,"Ox"))      return "SiO2";
  if(!strcmp(mat_name,"SiO2"))    return "SiO2";
  if(!strcmp(mat_name,"Nit"))     return "Nitride";
  if(!strcmp(mat_name,"Si"))      return "Si";
  if(!strcmp(mat_name,"Silicon")) return "Si";
  if(!strcmp(mat_name,"GaAs"))    return "GaAs";
  if(!strcmp(mat_name,"Ge"))      return "Ge";
  if(!strcmp(mat_name,"SiGe"))    return "SiGe";
  if(!strcmp(mat_name,"AlGaAs"))  return "AlGaAs";
  if(!strcmp(mat_name,"InGaAs"))  return "InGaAs";
  if(!strcmp(mat_name,"InGa"))    return "InGaAs";
  if(!strcmp(mat_name,"AlInAs"))  return "AlInAs";
  if(!strcmp(mat_name,"GaAsP"))   return "GaAsP";
  if(!strcmp(mat_name,"InGaP"))   return "InGaP";
  if(!strcmp(mat_name,"InAsP"))   return "InAsP";
  if(!strcmp(mat_name,"InP"))     return "InP";
  if(!strcmp(mat_name,"InAs"))    return "InAs";
  if(!strcmp(mat_name,"Poly"))    return "PolySi";
  if(!strcmp(mat_name,"SiC"))     return "SiC";
  if(!strcmp(mat_name,"ZnSe"))    return "ZnSe";
  if(!strcmp(mat_name,"ZnTe"))    return "ZnTe";
  if(!strcmp(mat_name,"HgCdTe"))  return "HgCdTe";
  if(!strcmp(mat_name,"TiSi2"))   return "TiSi2";
  return mat_name;
}

#endif


