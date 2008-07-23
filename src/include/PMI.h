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
/*  Last update: Feb 17 2008                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

// This is the head file of physical model ineface (PMI), which contains base class of:
//   physical model interface of semiconductor   PMIS
//   physical model interface of insulator       PMII 
//   physical model interface of conductor       PMIC 
//   physical model interface of vacuum          PMIV 
//   physical model interface of PML             PMIP 
// It links the main solver and material. The solver load required parameters from
// re-implemented virtual functions.  

#ifndef _PMIS_h_
#define _PMIS_h_
#include <complex>
#include <math.h>
#include "phy_scale.h"
#include "element.h"
#include "soldata.h"
#include "adolc.h"

using namespace adtl;


/*----------------------------------------------------------------------------
 *
 *               Physical Model Interface for Semiconductor 
 *
 *---------------------------------------------------------------------------*/


/* ----------------------------------------------------------------------------
 * PMIS_Environment, this class contains data for initializing PMIS Server.
 */
class PMIS_Environment
{
public:
  PhysicalUnitScale *pscale;       // physical unit
  Node              **ppnode;      // the location of current node
  SemiAuxData       **ppsemi_data; // the physical parameter data
  PetscScalar       *pclock;       // current time
public:
  PMIS_Environment(PhysicalUnitScale *scale, Node** node, SemiAuxData **data, PetscScalar *time)
  {
    pscale      = scale;
    ppsemi_data = data;
    ppnode      = node;
    pclock      = time;
  }
};



/* ----------------------------------------------------------------------------
 * PMIS_Server, the base class of PMIS
 */
class PMIS_Server
{
public:
  PetscScalar m,cm,um,g,kg,C,V,eV,A,mA,J,W,s,K;    //phusical unit
  PetscScalar kb,e,me,eps0,mu0,h;                  //physical constant
public:
  Node            **ppnode;      // the location of current node
  SemiAuxData     **ppsemi_data; // the physical parameter data
  PetscScalar     *pclock;       // current time
  // aux function return node coordinate.
  void   ReadCoordinate (PetscScalar& x, PetscScalar& y) const  { x = (*ppnode)->x; y = (*ppnode)->y;}
  // aux function return current time.
  PetscScalar ReadTime () const                                 { return *pclock;}
  // aux function return first mole function.
  PetscScalar ReadxMoleFraction () const                        { return (*ppsemi_data)->mole_x;}
  // aux function return first mole function with upper and lower bind.
  PetscScalar ReadxMoleFraction (const PetscScalar mole_xmin, const PetscScalar mole_xmax) const
  { 
     PetscScalar mole_x=(*ppsemi_data)->mole_x;
     if(mole_x<mole_xmin) return mole_xmin;
     if(mole_x>mole_xmax) return mole_xmax;
     return mole_x;
  }
  // aux function return second mole function.
  PetscScalar ReadyMoleFraction () const                        { return (*ppsemi_data)->mole_y;}
  // aux function return second mole function with upper and lower bind.
  PetscScalar ReadyMoleFraction (const PetscScalar mole_ymin, const PetscScalar mole_ymax) const
  { 
     PetscScalar mole_y=(*ppsemi_data)->mole_y;
     if(mole_y<mole_ymin) return mole_ymin;
     if(mole_y>mole_ymax) return mole_ymax;
     return mole_y;
  }
  // aux function return doping Na
  PetscScalar ReadDopingNa () const                             { return (*ppsemi_data)->Total_Na(); }
  // aux function return doping Nd
  PetscScalar ReadDopingNd () const                             { return (*ppsemi_data)->Total_Nd(); }

public:
  PMIS_Server(const PMIS_Environment &env)
  {
    m  = env.pscale->s_meter;
    cm = env.pscale->s_centimeter;
    um = env.pscale->s_micron;
    kg = env.pscale->s_kg;
    g  = 0.001*kg;
    C  = env.pscale->s_coulomb;
    V  = env.pscale->s_volt;
    eV = env.pscale->s_eV;
    A  = env.pscale->s_A;
    mA = env.pscale->s_mA;
    J  = env.pscale->s_joule;
    s  = env.pscale->s_second;
    W  = J*s;
    K  = env.pscale->s_kelvin;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;

    ppnode      = env.ppnode;
    ppsemi_data = env.ppsemi_data;
    pclock      = env.pclock;
  }
  virtual ~PMIS_Server(){};
};


/* ----------------------------------------------------------------------------
 * PMIS_BasicParameter. The PMIS interface for basic physical parameters of 
 * semiconductor material. User should implement each pure virtual functions.
 */
class PMIS_BasicParameter : public PMIS_Server
{
public:
  PMIS_BasicParameter(const PMIS_Environment &env):PMIS_Server(env) { }
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;
  virtual PetscScalar Permittivity  ()                      const=0;
  virtual PetscScalar Permeability  ()                      const=0;
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;
  virtual ~PMIS_BasicParameter() {};
};

/* ----------------------------------------------------------------------------
 * PMIS_BandStructure. The PMIS interface for band structure of 
 * semiconductor material. User should implement each pure virtual functions.
 */
class PMIS_BandStructure : public PMIS_Server
{
public:
  PMIS_BandStructure(const PMIS_Environment &env):PMIS_Server(env) { }
  virtual PetscScalar Eg             (const PetscScalar &Tl) =0;
  virtual AutoDScalar Eg             (const AutoDScalar &Tl) =0;
  virtual PetscScalar EgNarrow       (const PetscScalar &Tl) =0;
  virtual AutoDScalar EgNarrow       (const AutoDScalar &Tl) =0;
  virtual PetscScalar EgNarrowToEc   (const PetscScalar &Tl) =0;
  virtual PetscScalar EgNarrowToEv   (const PetscScalar &Tl) =0;
  virtual AutoDScalar EgNarrowToEc   (const AutoDScalar &Tl) =0;
  virtual AutoDScalar EgNarrowToEv   (const AutoDScalar &Tl) =0;
  
  virtual PetscScalar EffecElecMass  (const PetscScalar &Tl) =0;
  virtual PetscScalar EffecHoleMass  (const PetscScalar &Tl) =0;
  virtual PetscScalar Nc             (const PetscScalar &Tl) =0;
  virtual AutoDScalar Nc             (const AutoDScalar &Tl) =0;
  virtual PetscScalar Nv             (const PetscScalar &Tl) =0;
  virtual AutoDScalar Nv             (const AutoDScalar &Tl) =0;
  virtual PetscScalar nie            (const PetscScalar &Tl) =0;
  virtual AutoDScalar nie            (const AutoDScalar &Tl) =0;

  virtual PetscScalar TAUN           (const PetscScalar &Tl) =0;
  virtual PetscScalar TAUP           (const PetscScalar &Tl) =0;
  virtual PetscScalar Gamman         () {return 1.0;}
  virtual PetscScalar Gammap         () {return 1.0;}
  
  virtual PetscScalar R_Direct     (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;
  virtual AutoDScalar R_Direct     (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;
  
  virtual PetscScalar R_Auger      (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;
  virtual AutoDScalar R_Auger      (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;
  
  virtual PetscScalar R_SHR        (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;
  virtual AutoDScalar R_SHR        (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;
   
  virtual PetscScalar R_Surf       (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl, const PetscScalar &len) =0;
  virtual AutoDScalar R_Surf       (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl, const PetscScalar &len) =0;
  
  virtual PetscScalar Recomb       (const PetscScalar &p, const PetscScalar &n, const PetscScalar &Tl) =0;
  virtual AutoDScalar Recomb       (const AutoDScalar &p, const AutoDScalar &n, const AutoDScalar &Tl) =0;

  virtual PetscScalar ElecEnergyRelaxTime(const PetscScalar &Tn, const PetscScalar &Tl) =0;
  virtual AutoDScalar ElecEnergyRelaxTime(const AutoDScalar &Tn, const AutoDScalar &Tl) =0;
  virtual PetscScalar HoleEnergyRelaxTime(const PetscScalar &Tp, const PetscScalar &Tl) =0;
  virtual AutoDScalar HoleEnergyRelaxTime(const AutoDScalar &Tp, const AutoDScalar &Tl) =0;

  
  virtual PetscScalar SchottyJsn (PetscScalar n,PetscScalar Tl,PetscScalar Vb)=0;
  virtual PetscScalar SchottyJsp (PetscScalar p,PetscScalar Tl,PetscScalar Vb)=0;
  virtual PetscScalar SchottyBarrierLowerring (PetscScalar eps, PetscScalar E)=0;
  virtual PetscScalar pdSchottyJsn_pdn(PetscScalar n,PetscScalar Tl,PetscScalar Vb)=0;
  virtual PetscScalar pdSchottyJsp_pdp(PetscScalar p,PetscScalar Tl,PetscScalar Vb)=0;
  virtual PetscScalar pdSchottyJsn_pdTl(PetscScalar n,PetscScalar Tl,PetscScalar Vb)=0;
  virtual PetscScalar pdSchottyJsp_pdTl(PetscScalar p,PetscScalar Tl,PetscScalar Vb)=0;

  virtual PetscScalar ThermalVn (PetscScalar Tl)=0;
  virtual PetscScalar ThermalVp (PetscScalar Tl)=0;
  virtual PetscScalar pdThermalVn_pdTl (PetscScalar Tl)=0;
  virtual PetscScalar pdThermalVp_pdTl (PetscScalar Tl)=0;

  virtual PetscScalar BB_Tunneling(const PetscScalar &Tl, const PetscScalar &E) =0;
  virtual AutoDScalar BB_Tunneling(const AutoDScalar &Tl, const AutoDScalar &E) =0;


  virtual ~PMIS_BandStructure() {};
};


/* ----------------------------------------------------------------------------
 * PMIS_Mobility. The PMIS interface for semiconductor mobility.
 * User should implement each pure virtual functions.
 */
class PMIS_Mobility : public PMIS_Server
{
public:
  PMIS_Mobility(const PMIS_Environment &env):PMIS_Server(env) { }
  virtual PetscScalar ElecMob (const PetscScalar &p,  const PetscScalar &n,  const PetscScalar &Tl, 
                               const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tn) const=0;
  virtual PetscScalar HoleMob (const PetscScalar &p,  const PetscScalar &n,  const PetscScalar &Tl, 
                               const PetscScalar &Ep, const PetscScalar &Et, const PetscScalar &Tp) const=0;
  virtual AutoDScalar ElecMob (const AutoDScalar &p,  const AutoDScalar &n,  const AutoDScalar &Tl, 
                               const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tn) const=0;
  virtual AutoDScalar HoleMob (const AutoDScalar &p,  const AutoDScalar &n,  const AutoDScalar &Tl, 
                               const AutoDScalar &Ep, const AutoDScalar &Et, const AutoDScalar &Tp) const=0;
  virtual ~PMIS_Mobility() {}; 
};


/* ----------------------------------------------------------------------------
 * PMIS_Avalanche. The PMIS interface for semiconductor impact ionization parameter.
 * User should implement each pure virtual functions.
 */
class PMIS_Avalanche : public PMIS_Server
{
public:
  PMIS_Avalanche(const PMIS_Environment &env):PMIS_Server(env) { }
  //generation procedure for DDM
  virtual PetscScalar ElecGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const=0;
  virtual PetscScalar HoleGenRate (const PetscScalar &Tl,const PetscScalar &Ep,const PetscScalar &Eg) const=0;
  //Automatic Differentiation version for DDM
  virtual AutoDScalar ElecGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const=0;
  virtual AutoDScalar HoleGenRate (const AutoDScalar &Tl,const AutoDScalar &Ep,const AutoDScalar &Eg) const=0;
  
  //generation procedure for EBM
  virtual PetscScalar ElecGenRateEBM (const PetscScalar &Tn,const PetscScalar &Tl,const PetscScalar &Eg) const=0;
  virtual PetscScalar HoleGenRateEBM (const PetscScalar &Tp,const PetscScalar &Tl,const PetscScalar &Eg) const=0;
  //Automatic Differentiation version for EBM
  virtual AutoDScalar ElecGenRateEBM (const AutoDScalar &Tn,const AutoDScalar &Tl,const AutoDScalar &Eg) const=0;
  virtual AutoDScalar HoleGenRateEBM (const AutoDScalar &Tp,const AutoDScalar &Tl,const AutoDScalar &Eg) const=0;
  virtual ~PMIS_Avalanche() {};
};


/* ----------------------------------------------------------------------------
 * PMIS_Thermal. The PMIS interface for semiconductor thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIS_Thermal : public PMIS_Server
{
public:
  PMIS_Thermal(const PMIS_Environment &env):PMIS_Server(env) { }
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;
  virtual AutoDScalar HeatCapacity  (const AutoDScalar &Tl) const=0;
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;
  virtual AutoDScalar HeatConduction(const AutoDScalar &Tl) const=0;
  virtual ~PMIS_Thermal() {};
};


/* ----------------------------------------------------------------------------
 * PMIS_Optical. The PMIS interface for semiconductor optical refraction index
 * User should implement each pure virtual functions.
 */
class PMIS_Optical : public PMIS_Server
{
  public:
  PMIS_Optical(const PMIS_Environment &env):PMIS_Server(env) { }
  virtual complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0, PetscScalar Tl=1) const=0;
};  



/*----------------------------------------------------------------------------
 *
 *               Physical Model Interface for Insulator 
 *
 *---------------------------------------------------------------------------*/
 
/* ----------------------------------------------------------------------------
 * PMII_Environment, this class contains data for initializing PMII Server.
 */
class PMII_Environment
{
public:
  PhysicalUnitScale *pscale;      // physical unit
  Node              **ppnode;     // the location of current node
  ISAuxData         **ppis_data;  // the physical parameter data
  PetscScalar       *pclock;      // current time
public:
  PMII_Environment(PhysicalUnitScale *scale, Node** node, ISAuxData **data, PetscScalar *time)
  {
    pscale      = scale;
    ppis_data   = data;
    ppnode      = node;
    pclock      = time;
  }
};


/* ----------------------------------------------------------------------------
 * PMII_Server, the base class of PMII
 */
class PMII_Server
{
public:
  PetscScalar m,cm,um,g,kg,C,V,eV,A,mA,J,W,s,K;    //phusical unit
  PetscScalar kb,e,me,eps0,mu0,h;                  //physical constant
public:
  Node              **ppnode;    // the location of current node
  ISAuxData         **ppis_data; // the physical parameter data
  PetscScalar       *pclock;     // current time
  // aux function return node coordinate.
  void   ReadCoordinate (PetscScalar& x, PetscScalar& y) const { x = (*ppnode)->x; y = (*ppnode)->y;}
  // aux function return current time.
  PetscScalar ReadTime () const                                { return *pclock;}
public:
  PMII_Server(const PMII_Environment &env)
  {
    m  = env.pscale->s_meter;
    cm = env.pscale->s_centimeter;
    um = env.pscale->s_micron;
    kg = env.pscale->s_kg;
    g  = 0.001*kg;
    C  = env.pscale->s_coulomb;
    V  = env.pscale->s_volt;
    eV = env.pscale->s_eV;
    A  = env.pscale->s_A;
    mA = env.pscale->s_mA;
    J  = env.pscale->s_joule;
    s  = env.pscale->s_second;
    W  = J*s;
    K  = env.pscale->s_kelvin;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;

    ppnode      = env.ppnode;
    ppis_data   = env.ppis_data;
    pclock      = env.pclock;
  }
  virtual ~PMII_Server(){};
};


/* ----------------------------------------------------------------------------
 * PMII_BasicParameter. The PMII interface for basic physical parameters of 
 * insulator material. User should implement each pure virtual functions.
 */
class PMII_BasicParameter : public PMII_Server
{
public:
  PMII_BasicParameter(const PMII_Environment &env):PMII_Server(env) { }
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;
  virtual PetscScalar Permittivity  ()                      const=0;
  virtual PetscScalar Permeability  ()                      const=0;
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;
  virtual ~PMII_BasicParameter() {};
};


/* ----------------------------------------------------------------------------
 * PMII_Thermal. The PMII interface for insulator thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMII_Thermal : public PMII_Server
{
public:
  PMII_Thermal(const PMII_Environment &env):PMII_Server(env) { }
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;
  virtual ~PMII_Thermal() {};
};


/* ----------------------------------------------------------------------------
 * PMII_Optical. The PMII interface for insulator optical refraction index
 * User should implement each pure virtual functions.
 */
class PMII_Optical : public PMII_Server
{
  public:
  PMII_Optical(const PMII_Environment &env):PMII_Server(env) { }
  virtual complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0, PetscScalar Tl=1) const=0;
}; 



/*----------------------------------------------------------------------------
 *
 *               Physical Model Interface for Conductor
 *
 *---------------------------------------------------------------------------*/
 
/* ----------------------------------------------------------------------------
 * PMIC_Environment, this class contains data for initializing PMIC Server.
 */
class PMIC_Environment
{
public:
  PhysicalUnitScale *pscale;      // physical unit
  Node              **ppnode;     // the location of current node
  ELAuxData         **ppel_data;  // the physical parameter data
  PetscScalar       *pclock;      // current time
public:
  PMIC_Environment(PhysicalUnitScale *scale, Node** node, ELAuxData **data, PetscScalar *time)
  {
    pscale      = scale;
    ppel_data   = data;
    ppnode      = node;
    pclock      = time;
  }
};


/* ----------------------------------------------------------------------------
 * PMIC_Server, the base class of PMIC
 */
class PMIC_Server
{
public:
  PetscScalar m,cm,um,g,kg,C,V,eV,A,mA,J,W,s,K;    //phusical unit
  PetscScalar kb,e,me,eps0,mu0,h;                  //physical constant
public:
  Node              **ppnode;    // the location of current node
  ELAuxData         **ppel_data; // the physical parameter data
  PetscScalar            *pclock;     // current time
  // aux function return node coordinate.
  void   ReadCoordinate (PetscScalar& x, PetscScalar& y) const { x = (*ppnode)->x; y = (*ppnode)->y;}
  // aux function return current time.
  PetscScalar ReadTime () const                                { return *pclock;}
public:
  PMIC_Server(const PMIC_Environment &env)
  {
    m  = env.pscale->s_meter;
    cm = env.pscale->s_centimeter;
    um = env.pscale->s_micron;
    kg = env.pscale->s_kg;
    g  = 0.001*kg;
    C  = env.pscale->s_coulomb;
    V  = env.pscale->s_volt;
    eV = env.pscale->s_eV;
    A  = env.pscale->s_A;
    mA = env.pscale->s_mA;
    J  = env.pscale->s_joule;
    s  = env.pscale->s_second;
    W  = J*s;
    K  = env.pscale->s_kelvin;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;

    ppnode      = env.ppnode;
    ppel_data   = env.ppel_data;
    pclock      = env.pclock;
  }
  virtual ~PMIC_Server(){};
};


/* ----------------------------------------------------------------------------
 * PMIC_BasicParameter. The PMIC interface for basic physical parameters of 
 * conductor material. User should implement each pure virtual functions.
 */
class PMIC_BasicParameter : public PMIC_Server
{
public:
  PMIC_BasicParameter(const PMIC_Environment &env):PMIC_Server(env) { }
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;
  virtual PetscScalar Permittivity  ()                      const=0;
  virtual PetscScalar Permeability  ()                      const=0;
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;
  virtual ~PMIC_BasicParameter() {};
};


/* ----------------------------------------------------------------------------
 * PMIC_Thermal. The PMIC interface for conductor thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIC_Thermal : public PMIC_Server
{
public:
  PMIC_Thermal(const PMIC_Environment &env):PMIC_Server(env) { }
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;
  virtual ~PMIC_Thermal() {};
};


/* ----------------------------------------------------------------------------
 * PMIC_Optical. The PMIC interface for conductor optical refraction index
 * User should implement each pure virtual functions.
 */
class PMIC_Optical : public PMIC_Server
{
  public:
  PMIC_Optical(const PMIC_Environment &env):PMIC_Server(env) { }
  virtual complex<PetscScalar> RefractionIndex(PetscScalar lamda, PetscScalar Eg=0, PetscScalar Tl=1) const=0;
}; 


/*----------------------------------------------------------------------------
 *
 * Physical Model Interface for Vacuum, only used in EM FEM solver.
 *
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PMIV_Environment, this class contains data for initializing PMIV Server.
 */
class PMIV_Environment
{
public:
  PhysicalUnitScale *pscale;      // physical unit
  Node              **ppnode;     // the location of current node
  VacuumAuxData     **ppvac_data; // the physical parameter data
  PetscScalar       *pclock;      // current time
  
public:
  PMIV_Environment(PhysicalUnitScale *scale, Node** node, VacuumAuxData **data, PetscScalar *time)
  {
    pscale      = scale;
    ppvac_data  = data;
    ppnode      = node;
    pclock      = time;
  }
};


/*----------------------------------------------------------------------------
 * PMIV_Server, the base class of PMIV
 */
class PMIV_Server
{
public:
  PetscScalar m,cm,um,g,kg,C,V,eV,A,mA,J,W,s,K;    //phusical unit
  PetscScalar kb,e,me,eps0,mu0,h;                  //physical constant
public:
  Node              **ppnode;     // the location of current node
  VacuumAuxData     **ppvac_data; // the physical parameter data
  PetscScalar       *pclock;      // current time
  // aux function return node coordinate.
  void   ReadCoordinate (PetscScalar& x, PetscScalar& y) const  { x = (*ppnode)->x; y = (*ppnode)->y;}
  // aux function return current time.
  PetscScalar ReadTime () const                                 { return *pclock;}
public:
  PMIV_Server(const PMIV_Environment &env)
  {
    m  = env.pscale->s_meter;
    cm = env.pscale->s_centimeter;
    um = env.pscale->s_micron;
    kg = env.pscale->s_kg;
    g  = 0.001*kg;
    C  = env.pscale->s_coulomb;
    V  = env.pscale->s_volt;
    eV = env.pscale->s_eV;
    A  = env.pscale->s_A;
    mA = env.pscale->s_mA;
    J  = env.pscale->s_joule;
    s  = env.pscale->s_second;
    W  = J*s;
    K  = env.pscale->s_kelvin;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;

    ppnode      = env.ppnode;
    ppvac_data  = env.ppvac_data;
    pclock      = env.pclock;
  }
  virtual ~PMIV_Server(){};
};


/*----------------------------------------------------------------------------
 * PMIV_BasicParameter. The PMIV interface for basic physical parameters of 
 * vacuum. User should implement each pure virtual functions.
 */
class PMIV_BasicParameter : public PMIV_Server
{
public:
  PMIV_BasicParameter(const PMIV_Environment &env):PMIV_Server(env) { }
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;
  virtual PetscScalar Permittivity  ()                      const=0;
  virtual PetscScalar Permeability  ()                      const=0;
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;
  virtual ~PMIV_BasicParameter() {};
};


/* ----------------------------------------------------------------------------
 * PMIV_Thermal. The PMIV interface for vacuum thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIV_Thermal : public PMIV_Server
{
public:
  PMIV_Thermal(const PMIV_Environment &env):PMIV_Server(env) { }
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;
  virtual ~PMIV_Thermal() {};
};



/*----------------------------------------------------------------------------
 *
 * Physical Model Interface for PML Boundary, only used in EM FEM solver.
 *
 *---------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 * PMIP_Environment, this class contains data for initializing PMIP Server.
 */ 
class PMIP_Environment
{
public:
  PhysicalUnitScale *pscale;      // physical unit
  Node              **ppnode;     // the location of current node
  PMLAuxData        **ppml_data;  // the physical parameter data
  PetscScalar       *pclock;      // current time
public:
  PMIP_Environment(PhysicalUnitScale *scale, Node** node, PMLAuxData **data, PetscScalar *time)
  {
    pscale      = scale;
    ppml_data   = data;
    ppnode      = node;
    pclock      = time;
  }
};


/* ----------------------------------------------------------------------------
 * PMIP_Server, the base class of PMIP
 */
class PMIP_Server
{
public:
  PetscScalar m,cm,um,g,kg,C,V,eV,A,mA,J,W,s,K;    //phusical unit
  PetscScalar kb,e,me,eps0,mu0,h;                  //physical constant
public:
  Node              **ppnode;    // the location of current node
  PMLAuxData        **ppml_data; // the physical parameter data
  PetscScalar       *pclock;     // current time
  // aux function return node coordinate.
  void   ReadCoordinate (PetscScalar& x, PetscScalar& y) const  { x = (*ppnode)->x; y = (*ppnode)->y;}
  // aux function return node coordinate.
  PetscScalar ReadTime () const                                 { return *pclock;}
public:
  PMIP_Server(const PMIP_Environment &env)
  {
    m  = env.pscale->s_meter;
    cm = env.pscale->s_centimeter;
    um = env.pscale->s_micron;
    kg = env.pscale->s_kg;
    g  = 0.001*kg;
    C  = env.pscale->s_coulomb;
    V  = env.pscale->s_volt;
    eV = env.pscale->s_eV;
    A  = env.pscale->s_A;
    mA = env.pscale->s_mA;
    J  = env.pscale->s_joule;
    s  = env.pscale->s_second;
    W  = J*s;
    K  = env.pscale->s_kelvin;

    kb   = 1.3806503e-23*J/K;
    e    = 1.602176462e-19*C;
    me   = 9.10938188e-31*kg;
    eps0 = 8.854187818e-12*C/V/m;
    mu0  = 12.56637061e-7*pow(s,2)/C*V/m;
    h    = 6.62606876e-34*J*s;

    ppnode      = env.ppnode;
    ppml_data   = env.ppml_data;
    pclock      = env.pclock;
  }
  virtual ~PMIP_Server(){};
};


/*----------------------------------------------------------------------------
 * PMIP_BasicParameter. The PMIP interface for basic physical parameters of 
 * PML. User should implement each pure virtual functions.
 */
class PMIP_BasicParameter : public PMIP_Server
{
public:
  PMIP_BasicParameter(const PMIP_Environment &env):PMIP_Server(env) { }
  virtual PetscScalar Density       (const PetscScalar &Tl) const=0;
  virtual PetscScalar Permittivity  ()                      const=0;
  virtual PetscScalar Permeability  ()                      const=0;
  virtual PetscScalar Affinity      (const PetscScalar &Tl) const=0;
  virtual ~PMIP_BasicParameter() {};
};


/* ----------------------------------------------------------------------------
 * PMIP_Thermal. The PMIP interface for PML thermal parameter.
 * User should implement each pure virtual functions.
 */
class PMIP_Thermal : public PMIP_Server
{
public:
  PMIP_Thermal(const PMIP_Environment &env):PMIP_Server(env) { }
  virtual PetscScalar HeatCapacity  (const PetscScalar &Tl) const=0;
  virtual PetscScalar HeatConduction(const PetscScalar &Tl) const=0;
  virtual ~PMIP_Thermal() {};
};

#endif
