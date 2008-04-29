#ifndef _soldata_h_
#define _soldata_h_
#include "petsc.h"

//The solution and auxiliary data structure for each zone node

//----------------------------------------------------------------
class   SemiData
{
public:
  PetscScalar        n,p;           //density
  PetscScalar        P;             //electrostatic potential
  PetscScalar        T;             //lattice temperature
  PetscScalar        Tn,Tp;         //carrier temperature
  PetscScalar        Eqc,Eqv;       //quantum conduction/valence band 
  PetscScalar        n_last,p_last; //solution at previous step
  PetscScalar        P_last;
  PetscScalar        T_last;
  PetscScalar        Tn_last,Tp_last;
  PetscScalar        Eqc_last,Eqv_last;       
public:
  SemiData()
  {
    n=n_last=0.0;
    p=p_last=0.0;
    P=P_last=0.0;
    T=T_last=300.0;
    Tn=Tn_last=300.0;
    Tn=Tn_last=300.0;
    Eqc=Eqc_last=0.0;
    Eqv=Eqv_last=0.0;
  }
};


class   SemiAuxData
{
public:
  PetscScalar         density;
  PetscScalar         affinity;
  PetscScalar         Ec,Ev,Eg;
  PetscScalar         phin,phip;
  PetscScalar         phi_intrinsic;
  PetscScalar         Nc,Nv;
  PetscScalar         eps,mu;
  PetscScalar         Na,Nd;         //doping profile
  PetscScalar         mole_x,mole_y;
  PetscScalar         Ex,Ey;
  PetscScalar         mun,mup;
  PetscScalar         OpEx,OpEy,OpEz;
  PetscScalar         OpHx,OpHy,OpHz;
  PetscScalar         OptG;//dot not plot this variable
  PetscScalar         RealOptG;//for transient case
  SemiAuxData()
  {
    Na=Nd=0;
    Ex=Ey=0;
    phi_intrinsic=phin=phip=0;
    mole_x=mole_y=0;
    OpEx=OpEy=OpEz=0.0;
    OpHx=OpHy=OpHz=0.0;
    OptG=0.0;
    RealOptG=0.0;
  }
};



//----------------------------------------------------------------
class    ISData
{
public:
  PetscScalar          P;      //potential
  PetscScalar          T;      //lattice temperature
  PetscScalar          P_last; //solution at n-1 step
  PetscScalar          T_last;

public:
  ISData()
  {
    P=P_last=0.0;
    T=T_last=1.0;
  }
};

class   ISAuxData
{
public:
  PetscScalar          density;
  PetscScalar          eps,mu;
  PetscScalar          affinity;
  PetscScalar          Ex,Ey;
  PetscScalar          OpEx,OpEy,OpEz;
  PetscScalar          OpHx,OpHy,OpHz;
  ISAuxData()
  {
    Ex=Ey=0;
    OpEx=OpEy=OpEz=0.0;
    OpHx=OpHy=OpHz=0.0;
  }
};


//----------------------------------------------------------------
class  ELData
{
public:
  PetscScalar          P;      //potential
  PetscScalar          T;      //lattice temperature
  PetscScalar          P_last; //solution at n-1 step
  PetscScalar          T_last; //solution at n-1 step
  ELData()
  {
    P=P_last=0;
    T=T_last=1.0;
  }
} ;

class  ELAuxData
{
public:
  PetscScalar          density;
  PetscScalar          eps,mu;
  PetscScalar          affinity;
  PetscScalar          Ex,Ey;
  PetscScalar          OpEx,OpEy,OpEz;
  PetscScalar          OpHx,OpHy,OpHz;
  ELAuxData()
  {
    OpEx=OpEy=OpEz=0.0;
    OpHx=OpHy=OpHz=0.0;
    Ex=Ey=0.0;
  }
};


//----------------------------------------------------------------
class  VacuumData
{
public:
  PetscScalar          P;      //potential
  PetscScalar          P_last; //solution at n-1 step
  VacuumData()
  {
    P=P_last=0;
  }
} ;

class  VacuumAuxData
{
public:
  PetscScalar         eps,mu;
  PetscScalar         Ex,Ey;
  PetscScalar         OpEx,OpEy,OpEz;
  PetscScalar         OpHx,OpHy,OpHz;
  VacuumAuxData()
  {
    Ex=Ey=0;
    OpEx=OpEy=OpEz=0.0;
    OpHx=OpHy=OpHz=0.0;
  }
};


//----------------------------------------------------------------
class  PMLData
{
public:
  PetscScalar          P;      //potential
  PetscScalar          P_last; //solution at n-1 step
  PMLData()
  {
    P=P_last=0;
  }
} ;

class  PMLAuxData
{
public:
  PetscScalar         eps,mu;
  PetscScalar         Ex,Ey;
  PetscScalar         OpEx,OpEy,OpEz;
  PetscScalar         OpHx,OpHy,OpHz;
  PMLAuxData()
  {
    Ex=Ey=0;
    OpEx=OpEy=OpEz=0.0;
    OpHx=OpHy=OpHz=0.0;
  }
};
#endif
