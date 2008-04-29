// ----------------------------------------------------------------------------
//
// FILE:     grafix3d.h
//
// PROJECT:
//
// AUTHOR:   Kevin M. Kramer
//           Brian Chu, sing yun
//           Neil K. H
//           Eric R. Keiter
//           J. Joshua Feng
//
// REVISION HISTORY
// DATE      INITIALS   DESCRIPTION
// --------  --------   -------------------------------------------------------
// ??/??/??    KMK      initial implementation as triplot
// ??/??/??    BSC      change to sgplot (grid and gray level)
// ??/??/??    NKH      frame and grid
// ??/??/??    ERK      2D contour plot
// 04/20/97    JJF      correct the criterion for ordering / modules
//
// DESCRIPTION: TriPlot class member functions
//
// ----------------------------------------------------------------------------

#ifndef __grafix3d_h
#define __grafix3d_h

#include "termnlgy.h"
///////////////////////////////////////////////////////////////////////////////
///////////////  MISCELLANEOUS CONSTANTS, MACROS AND VARIABLES  ///////////////
///////////////////////////////////////////////////////////////////////////////

// definition of switches/flags
#define F_NO_MESH           0x0001
#define F_CON_SURF          0x0002
#define F_NO_SURF           0x0004
#define F_SCALE             0x0008
#define F_GRAYLEVEL         0x0010
#define F_COLOR             0x0020
#define F_WIREFRAME         0x0040

#define MIN_PIXEL_PER_LINE  50  // number of pixels needed for label width
#define LABELLINEFACTOR     4   // 1/4 tick spacing line offset
#define TICKFACTOR          2   // 1/2 line offset tick length
#define CHARFACTOR          2   // 1/2 of line offset(space label char fr tick)
#define RIGHTJUST           1   // right justify text3D output
#define LEFTJUST            0   // left justify text3D output
#define UPORIGIN            1   // upper left is char box origin
#define DOWNORIGIN          0   // lower left is char box origin

// font used on screen and in printing
#define CHARPIXELWIDE       7   // width of graphics char
#define CHARPIXELHIGH       14  // hight of graphics char(var if mode dep)
#define AXISRJCHARS         10  // # of chars to move left if right justify
#define AXISLJCHARS         5   // # of chars to move right if left justify
#define ZLABELUP            2   // # of chars to move up if labeling z axis
#define ZLABEL              1   // flag true if exponentLabel() should label z

///////////////////////////////////////////////////////////////////////////////

typedef UINT AINODE[4];                 // ainode

typedef struct tagPNT3D                 // pnt3D
{ double x;                               // x coordinate
  double y;                               // y coordinate
  double z;                               // z coordinate
} PNT3D;

typedef struct tagPNT2D                 // pnt2D
{ int x;                                // x coordinate
  int y;                                // y coordinate
} PNT2D;


class TriPlot3D
{ private:
    UINT cNodes;                        // number of nodes
    UINT cTriangles;                    // number of triangles
    AINODE *aainode;                    // node index of triangles
    PNT3D  *apnt3dInput;                // input data
    PNT3D  *apnt3dXform;                // transformed data
    PNT3D  *apnt3dMesh;                 // transformed mesh data
    PNT2D  *apnt2d;                     // screen data - surface plot
    PNT2D  *apnt2dMesh;                 // screen data - mesh
    UINT   *aiPlotOrder;                // plotting order
    BOOL   flogz;                       // log scaled data
  private:
    BOOL fScale;                        // Scale x, y, z data to fit the screen
    UINT Persp;                         // Perspective distance; 0: nonperspective
    void Xform3D_3D (PNT3D &, PNT3D &); // 3d to 3d transformation
    void Xform3D_2D (PNT3D &, PNT2D &); // 2d to 3d transformation
    void XformMesh  (PNT3D &, PNT3D &); // 3d to 3d transformation for mesh

  private:
    double AZ;                            // azimuthal rotation angle
    double EL;                            // elevation rotation angle
    double Cxx, Cxy, Cxz;                 // coefficients of xform matrix
    double Cyx, Cyy, Cyz;                 // coefficients of xform matrix
    double Czx, Czy, Czz;                 // coefficients of xform matrix
    PNT3D pnt3dMin;                     // minimum values
    PNT3D pnt3dMax;                     // maximum values
    PNT3D pnt3dScale;                   // scaling values

    PNT3D apnt3dCorner[8];              // Neil: without transformation

    PNT3D pnt3dMinXform;                // minimum values after transformation
    PNT3D pnt3dMaxXform;                // maximum values after transformation
    PNT3D pnt3dScaleXform;              // scaling values after transformation

  public:
    UINT cx;
    UINT cy;
    UINT cxLBorder;
    UINT cyLBorder;
    UINT cxRBorder;
    UINT cyRBorder;
    UINT nResFactor;                    // 1:screening  20:printing

  private:
    BOOL fGrayLevel;
    BOOL fColor;
    BOOL fWireframe;
    int  cMaxColor;
    void OpenWindow (int);                  // open graphics window
    void CloseWindow ();                    // close graphics window

  // ------------------ Contour stuff by Chu and Keiter -----------------
  private:
    int  NumContour;
    void DrawContour (PNT3D, PNT3D, PNT3D, double, BOOL);
    void line3D (PNT3D, PNT3D);          // plot 3D line

  // ------------------ Frame/Grid stuff by Neil KH ---------------------
  private:
    int  iztick;                         // number of ticks on z-axis
    void DrawFrame (UINT);               // draw background scale grid frame
    UINT DrawGridLines (UINT, int = -1); // draw gridlines (-1:autodetermine)
    void DrawFrameLabels (UINT, UINT, int = -1, int = -1);
                                         // draw label lines, tick marks, labels
    void text3D (PNT3D, char *, int = LEFTJUST, int = UPORIGIN);
                                         // put a textstring at a 3D point
    void exponentLabel(PNT3D, int, int = LEFTJUST, int = UPORIGIN, BOOL = FALSE);
                                         // write exponent label for x, y or z

  public:
    void GetViewAspect (double &, double &);        // prompt user for view aspect
    void DeterminePlotOrder ();             // determine plot order
    void Xform3DData ();                    // transform data: 3D rotation
    void Xform2DData ();                 // transform data: yz plane -> screen
    void PlotMesh ();                    // plot mesh
    void PlotData (BOOL, BOOL);          // plot data/contour lines (SYC & ERK)
    void PlotGrid ();                    // plot grid (NKH)

  public:
    TriPlot3D (UINT, UINT, AINODE *, PNT3D *, double, double, int, int,
               int = 0, int = -1, UINT = 0, BOOL=false);  // Constructor
   ~TriPlot3D ();                                         // Destructor
};

#endif
