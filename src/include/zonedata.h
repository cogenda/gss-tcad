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
/*  Last update: Feb 26, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#ifndef _zonedata_h_
#define _zonedata_h_
#include "bzonedata.h"
#include "petscsnes.h"
#include "soldata.h"
#include "solverdef.h"
#include "bc.h"

//----------------------------------------------------------------
//The Vacuum zone,only used for EM solver
//----------------------------------------------------------------
class VacuumZone: public BZoneData
{
public:
  VacuumData    *fs;
  VacuumAuxData *aux;
  MatVacuum     *mt;
  int           Init(ZONE*,double,PhysicalUnitScale *);
  void          report();
public:
   VacuumZone()  {mt=0; fs=0; aux=0;}
  ~VacuumZone()  {delete mt; delete [] fs; delete [] aux;}
};


//----------------------------------------------------------------
//The PML zone,only used for EM solver
//----------------------------------------------------------------
class PMLZone: public BZoneData
{
public:
  PMLData       *fs;
  PMLAuxData    *aux;
  MatPML        *mt;
  int           Init(ZONE*,double,PhysicalUnitScale *);
  void          report();
public:
   PMLZone()    {mt=0; fs=0; aux=0;}
  ~PMLZone()    {delete mt; delete [] fs; delete [] aux;}
};

//----------------------------------------------------------------
//an insulator zone
//----------------------------------------------------------------
class ISZone: public BZoneData
{
public:
  MatInsulator    *mt;
  ISData          *fs;
  ISAuxData       *aux;
  vector<int>   electrode;
public:
  void          F1_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F1_ddm_gatebc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1_ddm_semiconductor_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, SMCZone *,int);
  void          F1_ddm_electrode_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, ElZone *,int);
  void          F1_ddm_insulator_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, ISZone *,int);
  void          F1_ddm_chargebc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1_gate_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1_charge_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          J1_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J1_ddm_gatebc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1_ddm_semiconductor_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1_ddm_electrode_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          J1_ddm_insulator_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J1_ddm_chargebc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1_gate_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1_charge_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          F1_efield_update(PetscScalar *,vector<int> &, DABC &,vector<BZoneData *>);

  void          AC1_ddm_inner (int, PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &);
  void          AC1_ddm_gatebc(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &);
  void          AC1_ddm_chargebc(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &);
  void          AC1_ddm_semiconductor_insulator_interface(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &, SMCZone *,int);
  void          AC1_ddm_electrode_insulator_interface(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &, ElZone *,int);
  void          AC1_ddm_insulator_insulator_interface(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &, ISZone *,int);
  void          AC1_gate_electrode(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &, DABC &, PetscScalar);
  void          AC1_charge_electrode(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &, DABC &);
  void          AC1_gate_electrode_current(int ,PetscScalar, PetscScalar *, Vec *, vector<int> &, DABC &, PetscScalar,PetscScalar &, PetscScalar &);

  void          F1_mix_ddm_gatebc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          J1_mix_ddm_gatebc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          F1_mix_gate_electrode_current(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1_mix_gate_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);

  //--------------------------------------------------------------------------------------------------------
  void          F2_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F2_ddm_neumannbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2_ddm_gatebc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2_ddm_gatebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F2_ddm_chargebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          F2_ddm_semiconductor_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F2_ddm_electrode_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          F2_ddm_insulator_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F2_gate_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2_charge_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          J2_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J2_ddm_neumannbc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2_ddm_gatebc_segment(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2_ddm_gatebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J2_ddm_chargebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J2_ddm_semiconductor_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J2_ddm_electrode_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          J2_ddm_insulator_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J2_gate_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J2_charge_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          F2_efield_update(PetscScalar *,vector<int> &, DABC &,vector<BZoneData *>);

  void          F2_mix_ddm_gatebc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2_mix_ddm_gatebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J2_mix_ddm_gatebc_segment(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2_mix_ddm_gatebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F2_mix_gate_electrode_current(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2_mix_gate_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);

  //--------------------------------------------------------------------------------------------------------
  void          F3_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F3_ddm_neumannbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F3_ddm_gatebc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F3_ddm_gatebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F3_ddm_chargebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          F3_ddm_semiconductor_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F3_ddm_electrode_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          F3_ddm_insulator_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F3_gate_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F3_charge_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          J3_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J3_ddm_neumannbc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J3_ddm_gatebc_segment(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J3_ddm_gatebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J3_ddm_chargebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J3_ddm_semiconductor_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J3_ddm_electrode_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          J3_ddm_insulator_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J3_gate_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J3_charge_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          F3_efield_update(PetscScalar *,vector<int> &, DABC &,vector<BZoneData *>);

  //--------------------------------------------------------------------------------------------------------
  void          F1Q_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F1Q_ddm_gatebc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1Q_ddm_semiconductor_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, SMCZone *,int);
  void          F1Q_ddm_electrode_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, ElZone *,int);
  void          F1Q_ddm_insulator_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, ISZone *,int);
  void          F1Q_ddm_chargebc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1Q_gate_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1Q_charge_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          J1Q_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J1Q_ddm_gatebc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1Q_ddm_semiconductor_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1Q_ddm_electrode_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ElZone *,int);
  void          J1Q_ddm_insulator_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J1Q_ddm_chargebc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1Q_gate_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1Q_charge_electrode(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          F1Q_efield_update(PetscScalar *,vector<int> &, DABC &,vector<BZoneData *>);
  
  //--------------------------------------------------------------------------------------------------------  
public:
  int           import_solution(char *,char *,DABC &,PhysicalUnitScale *);
  int           export_solution(char *,char *,DABC &,PhysicalUnitScale *);
public:
  int           Init(ZONE*,double,PhysicalUnitScale *);
  void          report();
public:
  ISZone()      {mt=0; fs=0; aux=0;}	
  ~ISZone()     {delete mt; delete [] fs; delete [] aux;}
};


//----------------------------------------------------------------
//an electrode zone
//----------------------------------------------------------------
class ElZone: public BZoneData
{
public:
  ELData        *fs;
  ELAuxData     *aux;
  MatConductor  *mt;
public:
  vector<int>   electrode;
public:
  //--------------------------------------------------------------------------------------------------------  
  void          F1_ddm_inner(int , PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & );
  void          F1_ddm_om_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,SMCZone *, int );
  void          F1_ddm_stk_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,SMCZone *, int );
  void          F1_ddm_gate_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,ISZone *, int );
  void          F1_ddm_charge_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,ISZone *, int );
  void          J1_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J1_ddm_om_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1_ddm_stk_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1_ddm_gate_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J1_ddm_charge_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  //--------------------------------------------------------------------------------------------------------  
  void          AC1_ddm_inner (int, PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &);
  void          AC1_ddm_om_contact(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &, SMCZone *,int);
  void          AC1_ddm_stk_contact(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &, SMCZone *,int);
  void          AC1_ddm_gate_contact(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &, ISZone *,int);
  void          AC1_ddm_charge_contact(int ,PetscScalar, PetscScalar *, Vec *, Mat *, vector<int> &,DABC &, ISZone *,int);
  //--------------------------------------------------------------------------------------------------------  
  void          F2_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F2_ddm_neumannbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2_ddm_ombc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F2_ddm_stkbc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F2_ddm_gatebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F2_ddm_chargebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F2_ddm_elec_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J2_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J2_ddm_neumannbc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2_ddm_ombc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J2_ddm_stkbc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J2_ddm_gatebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J2_ddm_chargebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &,ISZone *,int);
  void          J2_ddm_elec_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  //--------------------------------------------------------------------------------------------------------  
  void          F3_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F3_ddm_neumannbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F3_ddm_ombc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F3_ddm_stkbc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F3_ddm_gatebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F3_ddm_chargebc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F3_ddm_elec_insulator_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J3_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J3_ddm_neumannbc(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J3_ddm_ombc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J3_ddm_stkbc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J3_ddm_gatebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J3_ddm_chargebc_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &,DABC &,ISZone *,int);
  void          J3_ddm_elec_insulator_interface(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  //--------------------------------------------------------------------------------------------------------  
  void          F1Q_ddm_inner(int , PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & );
  void          F1Q_ddm_om_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,SMCZone *, int );
  void          F1Q_ddm_stk_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,SMCZone *, int );
  void          F1Q_ddm_gate_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,ISZone *, int );
  void          F1Q_ddm_charge_contact(int ,PetscScalar *,PetscScalar *, ODE_Formula &, vector<int> & , DABC &,ISZone *, int );
  void          J1Q_ddm_inner(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &);
  void          J1Q_ddm_om_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1Q_ddm_stk_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1Q_ddm_gate_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J1Q_ddm_charge_contact(int ,PetscScalar *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
public:
  int           import_solution(char *,char *,DABC &,PhysicalUnitScale *);
  int           export_solution(char *,char *,DABC &,PhysicalUnitScale *);
public:
  int           Init(ZONE*,double,PhysicalUnitScale *);
  void          report();
public:
  ElZone()      {mt=0; fs=0; aux=0;}
  ~ElZone()     {delete mt; delete [] fs; delete [] aux;}
};


//----------------------------------------------------------------
//the semiconductor zone
//----------------------------------------------------------------
class SMCZone: public BZoneData
{
public:
  MatSemiconductor  *mt;
  SemiData          *fs;
  SemiAuxData       *aux;
  bool          HighFieldMobility;
  bool          ImpactIonization;
  bool          BandBandTunneling;
  bool          Fermi;
  bool          IncompleteIonization;
  bool          QuantumMechanical;
  bool          EJModel;
  bool          pad1;
  bool          pad2;
  int           IIType;
  double        QNFactor;
  double        QPFactor;


public:
  vector<int>   electrode;

  void          F1E_Tri_ddm(Tri *,PetscScalar *,PetscScalar *, vector<int> &);
  void          F1E_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F1E_ddm_ombc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1E_ddm_stkbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1E_ddm_insulator_gate(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &);
  void          F1E_ddm_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F1E_ddm_homojunction(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F1E_ddm_heterojunction(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F1E_om_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1E_stk_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1E_ins_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1E_Tri_ddm(Tri *,PetscScalar *,Mat *, vector<int> &);
  void          J1E_ddm_inner(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &);
  void          J1E_ddm_ombc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1E_ddm_stkbc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1E_ddm_insulator_gate(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          J1E_ddm_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J1E_ddm_homojunction(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1E_ddm_heterojunction(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J1E_om_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1E_stk_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1E_ins_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1E_efield_update(PetscScalar *, vector<int> &, DABC &, vector<BZoneData *>);
  void          F1E_om_electrode_Trace(int ,  PetscScalar &, Vec &, Vec &, Vec &,Mat *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);
  void          F1E_stk_electrode_Trace(int , PetscScalar &, Vec &, Vec &, Vec &,Mat *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);


  void          AC1_ddm_inner(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &);
  void          AC1_ddm_ombc(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &, DABC &);
  void          AC1_ddm_stkbc(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &, DABC &);
  void          AC1_ddm_insulator_gate(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &, DABC &);
  void          AC1_ddm_interface(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &, DABC &,ISZone *,int);
  void          AC1_ddm_homojunction(int ,PetscScalar,PetscScalar *, Mat *, Vec *, Mat *, vector<int> &, DABC &,SMCZone *,int);
  void          AC1_ddm_heterojunction(int ,PetscScalar,PetscScalar *, Mat *, Vec *, Mat *, vector<int> &, DABC &,SMCZone *,int);
  void          AC1_om_electrode(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &, DABC &, PetscScalar);
  void          AC1_stk_electrode(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &, DABC &, PetscScalar);
  void          AC1_ins_electrode(int, PetscScalar, PetscScalar *, Mat *,  Vec *, Mat *, vector<int> &, DABC &, PetscScalar);
  void          AC1_om_electrode_current(int, PetscScalar, PetscScalar *,Vec *, Mat *, vector<int> &, DABC &, PetscScalar, PetscScalar &, PetscScalar &);
  void          AC1_stk_electrode_current(int, PetscScalar, PetscScalar *,Vec *, Mat *, vector<int> &, DABC &, PetscScalar, PetscScalar &, PetscScalar &);
  void          AC1_ins_electrode_current(int, PetscScalar, PetscScalar *,Vec *, Mat *, vector<int> &, DABC &, PetscScalar, PetscScalar &, PetscScalar &);


  void          F1E_mix_ddm_ombc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1E_mix_ddm_stkbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1E_mix_ddm_insulator_gate(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &);
  void          J1E_mix_ddm_ombc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1E_mix_ddm_stkbc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1E_mix_ddm_insulator_gate(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          F1E_mix_om_electrode_current (int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1E_mix_stk_electrode_current(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1E_mix_ins_electrode_current(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1E_mix_om_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);
  void          F1E_mix_stk_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);
  void          F1E_mix_ins_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);

  void          F2E_Tri_ddm(Tri *,PetscScalar *,PetscScalar *, vector<int> &);
  void          F2E_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F2E_ddm_neumannbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2E_ddm_ombc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2E_ddm_ombc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F2E_ddm_stkbc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2E_ddm_stkbc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F2E_ddm_insulator_gate(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &);
  void          F2E_ddm_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F2E_ddm_homojunction(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F2E_ddm_heterojunction(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F2E_om_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2E_stk_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2E_ins_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J2E_Tri_ddm(Tri *,PetscScalar *,Mat *, vector<int> &);
  void          J2E_ddm_inner(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &);
  void          J2E_ddm_neumannbc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2E_ddm_ombc_segment(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2E_ddm_ombc_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J2E_ddm_stkbc_segment(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2E_ddm_stkbc_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J2E_ddm_insulator_gate(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          J2E_ddm_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J2E_ddm_homojunction(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J2E_ddm_heterojunction(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J2E_om_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J2E_stk_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J2E_ins_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2E_efield_update(PetscScalar *, vector<int> &, DABC &, vector<BZoneData *>);
  void          F2E_om_electrode_Trace(int ,  PetscScalar &, Vec &, Vec &, Vec &,Mat *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);
  void          F2E_stk_electrode_Trace(int , PetscScalar &, Vec &, Vec &, Vec &,Mat *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);

  void          F2E_mix_ddm_ombc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2E_mix_ddm_ombc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F2E_mix_ddm_stkbc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F2E_mix_ddm_stkbc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F2E_mix_ddm_insulator_gate(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &);
  void          J2E_mix_ddm_ombc_segment(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2E_mix_ddm_ombc_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J2E_mix_ddm_stkbc_segment(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J2E_mix_ddm_stkbc_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J2E_mix_ddm_insulator_gate(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          F2E_mix_om_electrode_current (int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2E_mix_stk_electrode_current(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2E_mix_ins_electrode_current(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F2E_mix_om_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);
  void          F2E_mix_stk_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);
  void          F2E_mix_ins_electrode_Load(int , double &, PetscScalar &, Vec &, Vec &, PetscScalar *,Mat *, ODE_Formula &, vector<int> &, DABC &,PetscScalar);


  void          F3E_Tri_ddm(Tri *,PetscScalar *,PetscScalar *, vector<int> &);
  void          J3E_Tri_ddm(Tri *,PetscScalar *,Mat *, vector<int> &);
  void          F3E_ddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F3E_ddm_neumannbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F3E_ddm_ombc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F3E_ddm_ombc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F3E_ddm_stkbc_segment(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F3E_ddm_stkbc_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          F3E_ddm_insulator_gate(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &);
  void          F3E_ddm_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F3E_ddm_homojunction(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  //void          F3E_ddm_heterojunction(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          F3E_om_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F3E_stk_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F3E_ins_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J3E_ddm_inner(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &);
  void          J3E_ddm_neumannbc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J3E_ddm_ombc_segment(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J3E_ddm_ombc_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J3E_ddm_stkbc_segment(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J3E_ddm_stkbc_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &,ElZone *,int);
  void          J3E_ddm_insulator_gate(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          J3E_ddm_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J3E_ddm_homojunction(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,SMCZone *,int);
  void          J3E_om_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J3E_stk_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J3E_ins_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F3E_efield_update(PetscScalar *, vector<int> &, DABC &, vector<BZoneData *>);

  void          F1Q_Tri_qddm(Tri *,PetscScalar *,PetscScalar *, vector<int> &);
  void          F1Q_qddm_inner(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &);
  void          F1Q_qddm_ombc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1Q_qddm_stkbc(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &,DABC &);
  void          F1Q_qddm_insulator_gate(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &);
  void          F1Q_qddm_interface(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          F1Q_om_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1Q_stk_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1Q_ins_electrode(int ,PetscScalar *,PetscScalar *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1Q_Tri_qddm(Tri *,PetscScalar *,Mat *, vector<int> &);
  void          J1Q_qddm_inner(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &);
  void          J1Q_qddm_ombc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1Q_qddm_stkbc(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &,DABC &);
  void          J1Q_qddm_insulator_gate(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &);
  void          J1Q_qddm_interface(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &,ISZone *,int);
  void          J1Q_om_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1Q_stk_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          J1Q_ins_electrode(int ,PetscScalar *,Mat *,Mat *,ODE_Formula &, vector<int> &, DABC &, PetscScalar);
  void          F1Q_efield_update(PetscScalar *, vector<int> &, DABC &, vector<BZoneData *>);
public:
  int           import_doping(char *, PhysicalUnitScale *);
  int           export_doping(char *, PhysicalUnitScale *);
  int           import_mole(char *,int);
  int           export_mole(char *,int);
  int           import_solution(char *, char *, DABC &, PhysicalUnitScale *);
  int           export_solution(char *, char *, DABC &, PhysicalUnitScale *);
public:
  int           Init(ZONE*,double, PhysicalUnitScale *);
  void          report();
public:
  SMCZone()
  {
    mt=0;
    fs=0;
    aux=0;
  }
  ~SMCZone()
  {
    delete mt;
    delete [] fs;
    delete [] aux;
  }
};



#endif
