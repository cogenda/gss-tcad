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
/*  Last update: Oct 17 2005                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _material_h_
#define _material_h_

#include "cmdbuf.h"
#include "phy_scale.h"
#include "element.h"
#include "petsc.h"
#include "solverdef.h"
#include "soldata.h"
#include "../material/matdefine.h"
#include "PMI.h"

int init_PMI_Semiconductor(CmdBuf *pcmdbuf, PMISDefine & PMIS_define);

class PhysicalConst
{
public:
  PetscScalar kb;      // Boltzmann constant
  PetscScalar e;       // elementary charge
  PetscScalar me;      // electron mass
  PetscScalar eps0;    // electric constant
  PetscScalar mu0;     // magnetic constant
  PetscScalar h;       // Planck constant
  PetscScalar hbar;    // Reduced Planck constant
public:
  PhysicalConst(PhysicalUnitScale *unit_scale);
};


//----------------------------------------------------------------
class MatSemiconductor: public PhysicalConst
{
public:
  PhysicalUnitScale   *pscale;
  // variables take to tcl script
public:
  Node                *pnode;
  SemiAuxData         *plocal_semi_data;
  PetscScalar          time;
  void                *dp;
  void *              (*set_ad_num)(const unsigned int);
  PMIS_BasicParameter *basic;
  PMIS_BandStructure  *band;
  PMIS_Mobility       *mob;
  PMIS_Avalanche      *gen;
  PMIS_Thermal        *thermal;
  PMIS_Optical        *optical;
  void mapping(Node*, SemiAuxData*, PetscScalar);
public:
  MatSemiconductor(const char *, const char *,PMISDefine &PMIS_define,PhysicalUnitScale *unit_scale);
  ~MatSemiconductor();
};


//----------------------------------------------------------------
class MatInsulator: public PhysicalConst
{
public:
  PhysicalUnitScale   *pscale;
  Node                *pnode;
  ISAuxData           *plocal_is_data;
  PetscScalar          time;
  void                *dp;
  PMII_BasicParameter *basic;
  PMII_Thermal        *thermal;
  PMII_Optical        *optical;
  void mapping(Node*, ISAuxData*, PetscScalar);
public:
  MatInsulator(const char *, const char *,PhysicalUnitScale *unit_scale);
  ~MatInsulator();
};


//----------------------------------------------------------------//```
class MatConductor: public PhysicalConst
{
public:
  PhysicalUnitScale   *pscale;
  Node                *pnode;
  ELAuxData           *plocal_el_data;
  PetscScalar          time;
  void                *dp;
  PMIC_BasicParameter *basic;
  PMIC_Thermal        *thermal;
  PMIC_Optical        *optical;
  void mapping(Node*,ELAuxData*,PetscScalar);
public:
  MatConductor(const char *, const char *,PhysicalUnitScale *unit_scale);
  ~MatConductor();
};


//----------------------------------------------------------------//```
class MatVacuum: public PhysicalConst
{
public:
  PhysicalUnitScale   *pscale;
  Node                *pnode;
  VacuumAuxData       *plocal_vac_data;
  PetscScalar          time;
  void                *dp;
  PMIV_BasicParameter *basic;
  PMIV_Thermal        *thermal;
  void mapping(Node*,VacuumAuxData*,PetscScalar);
public:
  MatVacuum(const char *, const char *,PhysicalUnitScale *unit_scale);
  ~MatVacuum();
};


//----------------------------------------------------------------
class MatPML: public PhysicalConst
{
public:
  PhysicalUnitScale   *pscale;
  Node                *pnode;
  PMLAuxData           *plocal_pml_data;
  PetscScalar          time;
  void                *dp;
  PMIP_BasicParameter *basic;
  PMIP_Thermal        *thermal;
  void mapping(Node*,PMLAuxData*,PetscScalar);
public:
  MatPML(const char *, const char *,PhysicalUnitScale *unit_scale);
  ~MatPML();
};

#endif
