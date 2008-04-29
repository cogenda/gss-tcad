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
/*  Last update: Nov 27, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _bc_h_
#define _bc_h_
#include "bzonedata.h"
#include "vsource.h"
#include "isource.h"
#include "lsource.h"
#include "phy_scale.h"
#include "typedef.h"
#include "petsc.h"
#include <complex>
using namespace std;

//boundary type
const int NeumannBoundary             = 0x0001;
const int IF_Semiconductor_Vacuum     = 0x0002;
const int IF_Insulator_Vacuum         = 0x0004;
const int IF_Electrode_Vacuum         = 0x0008;

const int IF_Electrode_Insulator      = 0x0011;//not the same as gate contact. 
const int InsulatorInterface          = 0x0012;
const int IF_Insulator_Semiconductor  = 0x0012;//the same as InsulatorInterface
const int IF_Insulator_Insulator      = 0x0013;
const int IF_Electrode_Electrode      = 0x0014;
const int HomoInterface               = 0x0015;
const int HeteroInterface             = 0x0016;
const int IF_Electrode_Semiconductor  = 0x0017;//electrode contact


const int OhmicContact        = 0x0101;
const int SchottkyContact     = 0x0102;
const int InsulatorContact    = 0x0104;
const int GateContact         = 0x0108;
const int ChargedContact      = 0x0111;
const int AbsorbingBoundary   = 0x1001;//by zhangxih 06-09-27
const int SourceBoundary      = 0x1002;//by zhangxih 06-09-27
//boundary mask
const int NeumannMask         = 0x000F;
const int InterfaceMask       = 0x00F0;
const int ContactMask         = 0x0F00;
//stimulant type of electrode
const int VoltageBC = 2005;
const int CurrentBC = 2006;


/* ----------------------------------------------------------------------------
 * BaseBC:  The base class of boundary condition
 */
class BaseBC
{
public:
  int           BCType;                  //support 5 types of Contact and insulator interface now
  Segment*      psegment;                //the pointer to the correlative segment
  ZONE   *      pzone;                   //the pointer to the correlative zone
  virtual       BaseBC* Get_this_pointer()=0; //for dynamic_cast operator
  virtual       void           Set_R(PetscScalar r) {}
  virtual       PetscScalar    Get_R() {return 0;}
  virtual       PetscScalar    Get_C() {return 0;}
  virtual       PetscScalar    Get_L() {return 0;}
  virtual       PetscScalar    Get_Vapp() {return 0;}
  virtual       PetscScalar    Get_Iapp() {return 0;}
  virtual       PetscScalar    Get_Potential() {return 0;}
  virtual       PetscScalar    Get_Potential_new() {return 0;}
  virtual       PetscScalar    Get_Potential_old() {return 0;}
  virtual       void    Set_Vapp(PetscScalar V) {}
  virtual       void    Set_Iapp(PetscScalar I) {}
  virtual       void    Set_Potential(PetscScalar V) {}
  virtual       void    Set_Potential_new(PetscScalar V) {}
  virtual       void    Set_Potential_old(PetscScalar V) {}
  virtual       void    Set_Current(PetscScalar I) {}
  virtual       void    Set_Current_new(PetscScalar I) {}
  virtual       void    Set_Current_old(PetscScalar I) {}
  virtual       PetscScalar    Get_Current() {return 0;}
  virtual       PetscScalar    Get_Current_new() {return 0;}
  virtual       PetscScalar    Get_Current_old() {return 0;}
  virtual       void    Set_Cap_Current(PetscScalar Ic) {}
  virtual       void    Set_electrode_type(int type) {}

  virtual       PetscScalar    Get_Vac()  {return 0;}
  virtual       void    Set_Vac(PetscScalar V)  {}
  virtual       complex<PetscScalar>    Get_Pac()  {return 0;}
  virtual       void    Set_Pac(complex<PetscScalar> P)  {}
  virtual       complex<PetscScalar>    Get_Iac()  {return 0;}
  virtual       void    Set_Iac(complex<PetscScalar> I)  {}
  virtual       ~BaseBC(){}
};

/* ----------------------------------------------------------------------------
 * NeumannBC:  The data structure of Neumann boundary
 */
class NeumannBC : public BaseBC
{
public:
  PetscScalar       T_external;          //external temperature
  PetscScalar       heat_transfer;       //heat transfer
  NeumannBC* Get_this_pointer()           //for dynamic_cast operator
  {return this;}
  NeumannBC();
};

/* ----------------------------------------------------------------------------
 * OhmicBC:  The data structure of Ohmic boundary
 */
class OhmicBC : public BaseBC
{
public:
  int              electrode_type;
  PetscScalar      R,C,L;                   //extern electrical parameters for electrode
  PetscScalar      Vapp;                    //the application voltage
  PetscScalar      Iapp;                    //the application current
  PetscScalar      Vac;                     //the application voltage for AC sweep.
  complex<PetscScalar>      Pac;            //the electrode potential for AC sweep.
  complex<PetscScalar>      Iac;            //the electrode current for AC sweep.
  PetscScalar      potential;               //the potential of this electrode
  PetscScalar      potential_new;           //the potential for current iterative cycle
  PetscScalar      potential_old;           //the potential for last step
  PetscScalar      current;                 //the current flow out of this electrode
  PetscScalar      current_new;             //the current for current iterative cycle
  PetscScalar      current_old;             //the current for last step
  PetscScalar      cap_current;             //the current which flow through lumped capacitance to ground
  PetscScalar      T_external;              //external temperature
  PetscScalar      heat_transfer;           //heat transfer
  Interface        *pinterface;
  //connector to another electrode. i.e. cmos structure
  string           connect_elec;            //connect bc label
  int              inner_connect;           //connect bc index  
  int              connect_zone;            //connect bc belongs to which zone
  BZoneData        *pzonedata;              //pointer to zonedata
  
  vector<VSource*> vsrc;                   //vsource point array belongs to this bc
  vector<ISource*> isrc;                   //isource point array belongs to this bc
  OhmicBC * Get_this_pointer()       //for dynamic_cast operator
  { return this;}
  void             Set_R(PetscScalar r) {R=r;}
  PetscScalar      Get_R() {return R;}
  PetscScalar      Get_C() {return C;}
  PetscScalar      Get_L() {return L;}
  PetscScalar      Get_Vapp() {return Vapp;}
  PetscScalar      Get_Iapp() {return Iapp;}
  void             Set_Vapp(PetscScalar V) {Vapp=V;}
  void             Set_Iapp(PetscScalar I) {Iapp=I;}
  PetscScalar      Get_Potential() {return potential;}
  PetscScalar      Get_Potential_new() {return potential_new;}
  PetscScalar      Get_Potential_old() {return potential_old;}
  void      Set_Potential(PetscScalar V) {potential=V;}
  void      Set_Potential_new(PetscScalar V) {potential_new=V;}
  void      Set_Potential_old(PetscScalar V) {potential_old=V;}
  void      Set_Current(PetscScalar I) {current=I;}
  void      Set_Current_new(PetscScalar I) {current_new=I;}
  void      Set_Current_old(PetscScalar I) {current_old=I;}
  PetscScalar      Get_Current()     {return current;}
  PetscScalar      Get_Current_new() {return current_new;}
  PetscScalar      Get_Current_old() {return current_old;}
  void      Set_Cap_Current(PetscScalar Ic) {cap_current=Ic;}
  void      Set_electrode_type(int type) {electrode_type=type;}
  PetscScalar      Get_Vac()  {return Vac;}
  void      Set_Vac(PetscScalar V)  {Vac=V;}
  complex<PetscScalar>    Get_Pac()  {return Pac;}
  void      Set_Pac(complex<PetscScalar> P)  {Pac=P;}
  complex<PetscScalar>    Get_Iac()  {return Iac;}
  void      Set_Iac(complex<PetscScalar> I)  {Iac=I;}
  OhmicBC();
};

/* ----------------------------------------------------------------------------
 * SchottkyBC:  The data structure of Schottky boundary
 */
class SchottkyBC : public BaseBC
{
public:
  int              electrode_type;
  PetscScalar      R,C,L;                   //extern electrical parameters for electrode
  PetscScalar      WorkFunction;            //workfunction of electrode material
  PetscScalar      Vapp;                    //the application voltage
  PetscScalar      Iapp;                    //the application current
  PetscScalar      Vac;                     //the application voltage for AC sweep.
  complex<PetscScalar>      Pac;            //the electrode potential for AC sweep.
  complex<PetscScalar>      Iac;            //the electrode current for AC sweep.
  PetscScalar      potential;               //the potential of this electrode
  PetscScalar      potential_new;           //the potential for current iterative cycle
  PetscScalar      potential_old;           //the potential for last step
  PetscScalar      current;                 //the current flow out of this electrode
  PetscScalar      current_new;             //the current for current iterative cycle
  PetscScalar      current_old;             //the current for last step
  PetscScalar      cap_current;             //the current which flow through lumped capacitance to ground
  PetscScalar      T_external;              //external temperature
  PetscScalar      heat_transfer;           //heat transfer
  Interface        *pinterface;
  //connector to another electrode. i.e. cmos structure
  string           connect_elec;            //connect bc label
  int              inner_connect;           //connect bc index  
  int              connect_zone;            //connect bc belongs to which zone
  BZoneData        *pzonedata;              //pointer to zonedata
  
  vector<VSource*> vsrc;                   //vsource point array belongs to this bc
  vector<ISource*> isrc;                   //isource point array belongs to this bc
  SchottkyBC *Get_this_pointer()           //for dynamic_cast operator
  {return this;}
  void             Set_R(PetscScalar r) {R=r;}
  PetscScalar      Get_R() {return R;}
  PetscScalar      Get_C() {return C;}
  PetscScalar      Get_L() {return L;}
  PetscScalar      Get_Vapp() {return Vapp;}
  PetscScalar      Get_Iapp() {return Iapp;}
  void      Set_Vapp(PetscScalar V) {Vapp=V;}
  void      Set_Iapp(PetscScalar I) {Iapp=I;}
  PetscScalar      Get_Potential() {return potential;}
  PetscScalar      Get_Potential_new() {return potential_new;}
  PetscScalar      Get_Potential_old() {return potential_old;}
  void      Set_Potential(PetscScalar V) {potential=V;}
  void      Set_Potential_new(PetscScalar V) {potential_new=V;}
  void      Set_Potential_old(PetscScalar V) {potential_old=V;}
  void      Set_Current(PetscScalar I) {current=I;}
  void      Set_Current_new(PetscScalar I) {current_new=I;}
  void      Set_Current_old(PetscScalar I) {current_old=I;}
  PetscScalar      Get_Current()     {return current;}
  PetscScalar      Get_Current_new() {return current_new;}
  PetscScalar      Get_Current_old() {return current_old;}
  void      Set_Cap_Current(PetscScalar Ic) {cap_current=Ic;}
  void      Set_electrode_type(int type) {electrode_type=type;}
  PetscScalar      Get_Vac()  {return Vac;}
  void      Set_Vac(PetscScalar V)  {Vac=V;}
  complex<PetscScalar>    Get_Pac()  {return Pac;}
  void      Set_Pac(complex<PetscScalar> P)  {Pac=P;}
  complex<PetscScalar>    Get_Iac()  {return Iac;}
  void      Set_Iac(complex<PetscScalar> I)  {Iac=I;}
  SchottkyBC();
};


/* ----------------------------------------------------------------------------
 * GateBC:  The data structure of MOS Gate.
 */
class GateBC : public BaseBC
{
public:
  int              electrode_type;
  PetscScalar      WorkFunction;          //workfunction of electrode material
  PetscScalar      R,C,L;                 //extern electrical parameters for electrode
  PetscScalar      Vapp;                  //the application voltage
  PetscScalar      Vac;                     //the application voltage for AC sweep.
  complex<PetscScalar>      Pac;            //the electrode potential for AC sweep.
  complex<PetscScalar>      Iac;            //the electrode current for AC sweep.
  PetscScalar      potential;             //the potential of this electrode
  PetscScalar      potential_new;           //the potential for current iterative cycle
  PetscScalar      potential_old;           //the potential for last step
  PetscScalar      current;               //the current flow out of this electrode
  PetscScalar      current_new;           //the current for current iterative cycle
  PetscScalar      current_old;             //the current for last step
  PetscScalar      cap_current;           //the current which flow through lumped capacitance to ground
  PetscScalar      T_external;            //external temperature
  PetscScalar      heat_transfer;         //heat transfer
  Interface        *pinterface;
  vector<VSource*> vsrc;           //vsource point array belongs to this bc
  GateBC*   Get_this_pointer()     //for dynamic_cast operator
  { return this;}
  void             Set_R(PetscScalar r) {R=r;}
  PetscScalar      Get_R() {return R;}
  PetscScalar      Get_C() {return C;}
  PetscScalar      Get_L() {return L;}
  PetscScalar      Get_Vapp() {return Vapp;}
  void      Set_Vapp(PetscScalar V) {Vapp=V;}
  PetscScalar      Get_Potential() {return potential;}
  PetscScalar      Get_Potential_new() {return potential_new;}
  PetscScalar      Get_Potential_old() {return potential_old;}
  void      Set_Potential(PetscScalar V) {potential=V;}
  void      Set_Potential_new(PetscScalar V) {potential_new=V;}
  void      Set_Potential_old(PetscScalar V) {potential_old=V;}
  void      Set_Current(PetscScalar I) {current=I;}
  void      Set_Current_new(PetscScalar I) {current_new=I;}
  void      Set_Current_old(PetscScalar I) {current_old=I;}
  PetscScalar      Get_Current()     {return current;}
  PetscScalar      Get_Current_new() {return current_new;}
  PetscScalar      Get_Current_old() {return current_old;}
  void      Set_Cap_Current(PetscScalar Ic) {cap_current=Ic;}
  void      Set_electrode_type(int type) {electrode_type=type;}
  PetscScalar      Get_Vac()  {return Vac;}
  void      Set_Vac(PetscScalar V)  {Vac=V;}
  complex<PetscScalar>    Get_Pac()  {return Pac;}
  void      Set_Pac(complex<PetscScalar> P)  {Pac=P;}
  complex<PetscScalar>    Get_Iac()  {return Iac;}
  void      Set_Iac(complex<PetscScalar> I)  {Iac=I;}
  GateBC();
};


/* ----------------------------------------------------------------------------
 * InsulatorContactBC:  The data structure of simplified MOS Gate
 */
class InsulatorContactBC : public BaseBC
{
public:
  int              electrode_type;
  PetscScalar      WorkFunction;            //workfunction of electrode material
  PetscScalar      Thick;                   //for simple SiSiO2 interface
  PetscScalar      eps;                     //relative dielectric permittivity of insulator
  PetscScalar      R,C,L;                   //extern electrical parameters for electrode
  PetscScalar      Vapp;                    //the application voltage
  PetscScalar      Vac;                     //the application voltage for AC sweep.
  complex<PetscScalar>      Pac;            //the electrode potential for AC sweep.
  complex<PetscScalar>      Iac;            //the electrode current for AC sweep.
  PetscScalar      potential;               //the potential of this electrode
  PetscScalar      potential_new;           //the potential for current iterative cycle
  PetscScalar      potential_old;           //the potential for last step
  PetscScalar      current;                 //the current flow out of this electrode
  PetscScalar      current_new;             //the current for current iterative cycle
  PetscScalar      current_old;             //the current for last step
  PetscScalar      cap_current;             //the current which flow through lumped capacitance to ground
  PetscScalar      QF;
  PetscScalar      T_external;              //external temperature
  PetscScalar      heat_transfer;           //heat transfer
  vector<VSource*> vsrc;                   //vsource point array belongs to this bc
  InsulatorContactBC *  Get_this_pointer()    //for dynamic_cast operator
  {return this;}
  void             Set_R(PetscScalar r) {R=r;}
  PetscScalar      Get_R() {return R;}
  PetscScalar      Get_C() {return C;}
  PetscScalar      Get_L() {return L;}
  PetscScalar      Get_Vapp() {return Vapp;}
  void             Set_Vapp(PetscScalar V) {Vapp=V;}
  PetscScalar      Get_Potential() {return potential;}
  PetscScalar      Get_Potential_new() {return potential_new;}
  PetscScalar      Get_Potential_old() {return potential_old;}
  void             Set_Potential(PetscScalar V) {potential=V;}
  void             Set_Potential_new(PetscScalar V) {potential_new=V;}
  void             Set_Potential_old(PetscScalar V) {potential_old=V;}
  void             Set_Current(PetscScalar I) {current=I;}
  void             Set_Current_new(PetscScalar I) {current_new=I;}
  void             Set_Current_old(PetscScalar I) {current_old=I;}
  PetscScalar      Get_Current()     {return current;}
  PetscScalar      Get_Current_new() {return current_new;}
  PetscScalar      Get_Current_old() {return current_old;}
  void             Set_Cap_Current(PetscScalar Ic) {cap_current=Ic;}
  void             Set_electrode_type(int type) {electrode_type=type;}
  PetscScalar      Get_Vac()  {return Vac;}
  void             Set_Vac(PetscScalar V)  {Vac=V;}
  complex<PetscScalar>    Get_Pac()  {return Pac;}
  void             Set_Pac(complex<PetscScalar> P)  {Pac=P;}
  complex<PetscScalar>    Get_Iac()  {return Iac;}
  void             Set_Iac(complex<PetscScalar> I)  {Iac=I;}
  InsulatorContactBC();
};


/* ----------------------------------------------------------------------------
 * ChargedContactBC:  The data structure of float matel boundary 
 */
class ChargedContactBC : public BaseBC
{
public:
  PetscScalar      QF;                      //the free charge of boundary
  PetscScalar      potential;               //the potential of float gate
  PetscScalar      potential_new;           //the potential for current iterative cycle
  PetscScalar      potential_old;           //the potential for last step
  Interface *pinterface;
  PetscScalar      Get_Potential() {return potential;}
  PetscScalar      Get_Potential_new() {return potential_new;}
  PetscScalar      Get_Potential_old() {return potential_old;}
  void             Set_Potential(PetscScalar V) {potential=V;}
  void             Set_Potential_new(PetscScalar V) {potential_new=V;}
  void             Set_Potential_old(PetscScalar V) {potential_old=V;}
  ChargedContactBC* Get_this_pointer()   {return this;} //for dynamic_cast operator
  ChargedContactBC();
};



/* ----------------------------------------------------------------------------
 * ElectrodeElectrodeBC:  The data structure of Interface between Electrode and Electrode
 */
class ElectrodeElectrodeBC : public BaseBC
{
public:
  Interface *pinterface;
  ElectrodeElectrodeBC* Get_this_pointer()   //for dynamic_cast operator
  {return this;}
  ElectrodeElectrodeBC(){};
};


/* ----------------------------------------------------------------------------
 * ElectrodeInsulatorBC:  The data structure of interface between region of 
 * (Ohmic or Schottky) Electrode and Insulator. 
 */
class ElectrodeInsulatorBC : public BaseBC
{
public:
  Interface *pinterface;
  ElectrodeInsulatorBC* Get_this_pointer()   {return this;} //for dynamic_cast operator
  ElectrodeInsulatorBC(){};
};


/* ----------------------------------------------------------------------------
 * InsulatorInterfaceBC:  The data structure of Interface between SiO2 and Si bulk
 */
class InsulatorInterfaceBC : public BaseBC
{
public:
  PetscScalar      QF;                      //InsulatorInterface fixed charge density
  Interface *pinterface;
  InsulatorInterfaceBC* Get_this_pointer()   //for dynamic_cast operator
  {return this;}
  InsulatorInterfaceBC();
};


/* ----------------------------------------------------------------------------
 * InsulatorInsulatorBC:  The data structure of Interface between Insulator and Insulator
 */
class InsulatorInsulatorBC : public BaseBC
{
public:
  Interface *pinterface;
  InsulatorInsulatorBC* Get_this_pointer()   //for dynamic_cast operator
  {return this;}
  InsulatorInsulatorBC(){};
};


/* ----------------------------------------------------------------------------
 * HomoInterfaceBC:  The data structure of interface between region of same
 * semiconductor material, mainly for further parallel code.
 */
class HomoInterfaceBC : public BaseBC
{
public:
  Interface *pinterface;
  HomoInterfaceBC* Get_this_pointer()   {return this;} //for dynamic_cast operator
  HomoInterfaceBC(){};
};


/* ----------------------------------------------------------------------------
 * HeteroInterfaceBC:  The data structure of Heterojuction
 */
class HeteroInterfaceBC : public BaseBC
{
public:
  PetscScalar      QF;                      //InsulatorInterface fixed charge density
  Interface *pinterface;
  HeteroInterfaceBC* Get_this_pointer()   //for dynamic_cast operator
  {return this;}
  HeteroInterfaceBC();
};


/* ----------------------------------------------------------------------------
 * AbsorbingBC:  The data structure of absorbing boundary
 */
class AbsorbingBC : public BaseBC //by zhangxih ---06-09-27
{
public:
  PetscScalar       T_external;          //external temperature
  PetscScalar       heat_transfer;       //heat transfer
  AbsorbingBC* Get_this_pointer()         //for dynamic_cast operator
  {return this;}
  AbsorbingBC();
};


/* ----------------------------------------------------------------------------
 * SourceBC:  The data structure of source boundary
 */
class SourceBC : public BaseBC //by zhangxih ---06-09-27
{
public:
  PetscScalar       T_external;          //external temperature
  PetscScalar       heat_transfer;       //heat transfer
  SourceBC* Get_this_pointer()    //for dynamic_cast operator
  {return this;}
  SourceBC();
};


/* ----------------------------------------------------------------------------
 * DABC:  The dynamic array of boundary condition
 */
class DABC
{
protected:
  int  zone_num;
  int   bc_num;
  vector<BaseBC*> bc_point_array;                   // dynamic array for bc point
  Segment** psegment_array;
  PhysicalUnitScale *scale;
  ZoneInterface *pzintface;
  double  lattice_temperature;
protected:
  //these functions set BC by infos from CmdParseBuf
  int SetBCNeumannBoundary(list<Cmd>::iterator);
  int SetBCOhmicContact(list<Cmd>::iterator);
  int SetBCInsulatorInterface(list<Cmd>::iterator);
  int SetBCHeterojunction(list<Cmd>::iterator);
  int SetBCInsulatorContact(list<Cmd>::iterator);
  int SetBCSchottkyContact(list<Cmd>::iterator);
  int SetBCGateContact(list<Cmd>::iterator);
  int SetBCChargedContact(list<Cmd>::iterator);
  int SetBCAbsorbingBoundary(list<Cmd>::iterator);
  int SetBCSourceBoundary(list<Cmd>::iterator);
  
  int SetElectrodeOhmicContact(list<Cmd>::iterator,int,ZONE*);
  int SetElectrodeSchottkyContact(list<Cmd>::iterator,int,ZONE*);
  int SetElectrodeGateContact(list<Cmd>::iterator,int,ZONE*);
  int SetFloatMetal(list<Cmd>::iterator,int,ZONE*);
public:  
  int      size();
  void     clear();
  BaseBC * Get_pointer(int );
  int      Get_bc_index  (int, const char * );
  int      Get_bc_index(const char *, const char *);
  int      Get_bc_index  (const char * );
  int      Get_bc_index_nocase  (const char * );
  
  int      is_electrode(int );                       
  int      is_electrode(const char * );   
  BaseBC * Get_electrode_pointer(const char * );
  BaseBC * Get_electrode_pointer_nocase(const char * );
  int      is_electrode_label(int, const char * );
  int      is_electrode_label_nocase(int, const char * );
  const char * format_electrode_name(const char *);   
  const char * format_electrode_name_nocase(const char *);                    
  void     Update_Vapp(PetscScalar );                 //set electrode application voltage by vsrc
  void     Clear_Vapp();
  void     Update_Iapp(PetscScalar );                 //set electrode application current by vsrc
  void     Clear_Iapp();
  void     Set_Vapp(const char *, PetscScalar);       //set electrode application voltage by user specified value
  void     Set_Vapp_nocase(const char *, PetscScalar);
  void     Set_Vac(const char *, PetscScalar);
  void     Set_Iapp(const char *, PetscScalar);
  void     Attach_Vapp(const char *, vector<VSource *> &);
  void     Attach_Iapp(const char *, vector<ISource *> &);
  void     Set_electrode_type(const char *,int);
  const BaseBC & operator[](int );
  int      EmitTo(int zone_index, PetscScalar source_x, PetscScalar source_y, PetscScalar ex, PetscScalar ey); 
public:
  int      InitBC(int,Segment**,int,ZONE*,ZoneInterface &,CmdBuf *,double, PhysicalUnitScale *);
public:
  DABC();
  ~DABC();
};

#endif
