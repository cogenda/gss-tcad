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
/*  Last update: Jan 19, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _solverdef_h_
#define _solverdef_h_
#include "typedef.h"
#include <stdio.h>
#include <vector>
#include <string>
using namespace std;

//define cgns file type
const int MODEL = 101;
const int CORE  = 102;

//define carrier type
const int N_Type  = 201;
const int P_Type  = 202;
const int PN_Both = 203;
const int PN_None = 204;

//define solver type
const int POISSON = 300;

const int DDML1   = 310;
const int DDML1E  = 311;
const int DDML1AC = 312;
const int DDML1MIX= 313;
const int QDDML1E = 315;

const int DDML2   = 320;
const int DDML2E  = 321;
const int DDML2MIX= 322;

const int EBML3E  = 331;

const int EMFEM   = 391;

//define scheme
const int DDM_Gummel   = 501;
const int DDM_Newton   = 502;

//define solve type
const int EQUILIBRIUM    =  601;
const int STEADYSTATE    =  602;
const int DCSWEEP        =  603;
const int DCSWEEP_VSCAN  =  604;
const int DCSWEEP_ISCAN  =  605;
const int TRANSIENT      =  606;
const int ACSWEEP        =  607;
const int TRACE          =  608;

//define parameter index
const int MeshObject     = 700;
const int ElecDensity    = 701;
const int HoleDensity    = 702;
const int EFieldX        = 703;
const int EFieldY        = 704;
const int EField         = 705;
const int Potential      = 706;
const int Phi_Intrinsic  = 707;
const int PhiN           = 708;
const int PhiP           = 709;
const int ValenceBand    = 710;
const int ConductionBand = 711;
const int QuantumEc      = 712;
const int QuantumEv      = 713;

const int Doping      = 720;
const int DopingNa    = 721;
const int DopingNd    = 722;
const int NetDoping   = 723;
const int Phosphorus  = 724;
const int Arsenic     = 725;
const int Antimony    = 726;
const int Boron       = 727;

const int Temperature = 730;
const int ElecTemp    = 731;
const int HoleTemp    = 732;

const int OpticalEx   = 750;//```
const int OpticalEy   = 751;
const int OpticalEz   = 752;
const int OpticalHx   = 753;
const int OpticalHy   = 754;
const int OpticalHz   = 755;
const int OpticalG    = 756;

//define order for ODE solver
const int BDF1 = 1;
const int BDF2 = 2;

const int Linear    = 851;
const int SignedLog = 852;

//define Newton's method
const int Basic        = 915;
const int LineSearch   = 917;
const int TrustRegion  = 910;

//define Newton damping
const int DampingNo         = 930;
const int DampingBankRose   = 931;
const int DampingPotential  = 932;

//define Impact Ionization Type
const int EdotJ   = 0;
const int GradQf  = 1;
const int EVector = 2;
const int ESide   = 3;
const int SoftII  = 4;
const int HardII  = 5;
const int TempII  = 6; 

class SolveDefine
{
public:
  int     Solver;       //DDM model level 1,2...
  int     Scheme;       //solver scheme: DDM_Newton, DDM_Gummel...
  int     Type;         //solver type:   transient , steadystate...
  int     NS;           //nonlinear solver scheme: basic, line search, trust region...
  int     Damping;      //Newton damping
  string  LS;           //linear solver scheme: LU, BCGS, GMRES ...
  int     MaxIteration; //Max Iteration number
  int     port;         //the net port used for mix simulation
  int     BDF_Type;     //transient solver: ODE formula
  bool    Projection;   //use solution projection in DC sweep
  // semiconductor region advanced model specification
  bool    HighFieldMobility;
  bool    BandBandTunneling;
  bool    ImpactIonization;
  bool    Fermi;
  bool    IncompleteIonization;
  bool    QuantumMechanical;
  bool    EJModel;
  double  QNFactor;
  double  QPFactor;
  int     IIType;
  // nonlinear solver convergence criteria
  int         maxit;
  double      relative_toler;
  double      toler_relax;
  double      possion_abs_toler;
  double      elec_continuty_abs_toler;
  double      hole_continuty_abs_toler;
  double      heat_equation_abs_toler;
  double      elec_energy_abs_toler;
  double      hole_energy_abs_toler;
  double      elec_quantum_abs_toler;
  double      hole_quantum_abs_toler;
  double      electrode_abs_toler;
  
  // parameters for DC, TRAN, TRACE etc.
  int               Electrode_VScan;
  vector<string>    Electrode_VScan_Name;
  int               Electrode_IScan;
  string            Electrode_IScan_Name;
  int               Electrode_ACScan;
  string            Electrode_ACScan_Name;
  vector<int>       Electrode_Record;
  vector<int>       Electrode_Record_Index;
  vector<string>    Electrode_Record_Name;
  string    IVFile;
  bool      IVFileAppend;
  double    VStart,VStep,VStepMax,VStop;
  double    IStart,IStep,IStop;
  double    TStart,TStep,TStop;
  bool      AutoStep;
  bool      Predict;
  double    VAC;
  double    FStart,FMultiple,FStop;
  SolveDefine()
  {
    Solver = DDML1;
    Scheme = DDM_Newton;
    NS = LineSearch;
    Damping = DampingNo;
    Projection = false;
    LS = "gmres";
    Fermi=false;
    IncompleteIonization=false;
    ImpactIonization = false;
    EJModel = false;
    BandBandTunneling = false;
    BDF_Type = BDF1;
    AutoStep=true;
    Predict=true;
  }
};

class PMISDefine
{
public:
  vector<string>  Region;
  vector<string>  Material;
  vector<string>  Mobility;
  vector<string>  IIModel;
  vector<string>  OpticalModel;
};

typedef struct
{
  int     TimeDependent;
  int     BDF_Type;
  int     BDF2_restart;
  PetscScalar    clock;
  PetscScalar    dt;
  PetscScalar    dt_last;
  PetscScalar    dt_last_last;
}
ODE_Formula;

typedef struct
{
  int     Variable;
  int     Measure;
  double  Dispersion;
  double  DivisionRatio;
  char    tri_cmd[32];
}
RefineDefine;

typedef struct
{
  bool    core;
  bool    ascii;
  bool    vtk;
  bool    acceptor;
  bool    donor;
  char    CoreFile[32];
  char    AscFile[32];
  char    VTKFile[32];
  char    AcceptorFile[32];
  char    DonorFile[32];
}
ExportDefine;

typedef struct
{
  int     file_type;
  char    CoreFile[32];
  char    ModelFile[32];
}
ImportDefine;

typedef struct
{
  char    Electrode[32];
  int     bc_index;
  int     electrode_type;
  vector<string> vsrc_name;
  vector<int>    vsrc_index;
  vector<string> isrc_name;
  vector<int>    isrc_index;
}
AttachDefine;

typedef struct
{
  int    MeshOnly;            // draw mesh only (2D)
  int    Variable;            // which Variable to plot
  string VariableName;        // the name of the variable to plot.
  int    Resolution;          // x window resolution
  int    Measure;             // log or linear plot
  int    GeneratePS;          // generate postscript file
  int    GenerateTIFF;        // generate TIFF file
  string PSFileName;
  string TIFFFileName;
  int    NoMesh;              // do not plot the mesh
  int    ConSurf;
  int    NoSurf;
  int    Style;               // STYLE_BW | STYLE_WIREFRAME | STYLE_COLOR | (default) | STYLE_GRAYLEVEL
  int    ztick;               // # of z axis ticks
  int    Persp;               // Perspective distance

  int    Con;                 // contour plot
  int    NumCon;              // # of contours

  double az;                // azimuthal angle
  double el;                // elevation angle
}
PlotDefine;


class ProbeDefine
{
public:
  string  Region;
  string  Segment;
  int     zone_index;
  int     bc_index;
  string  ProbeFile;
  FILE*   pFile;
  bool    Append;
  bool    flag;
  int     Variable;            // which Variable to probe
  string  VariableName;        // the name of the variable to probe.

  ProbeDefine()
  {
    pFile=NULL;
    flag=false;
  }
};

#endif
