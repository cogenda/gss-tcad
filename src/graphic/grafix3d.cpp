// ----------------------------------------------------------------------------
//
// FILE:     grafix3d.cpp
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
// 06/17/97    JJF      modification/separate into modules/correct the sorting
//
// DESCRIPTION: TriPlot class member functions
//
// ----------------------------------------------------------------------------
#include "config.h"
#if (defined(HAVE_WIN32) || defined(HAVE_X11))

#if defined(HAVE_WIN32)
  #include <io.h>
  #include <windows.h>
  #include "wgraph.h"
#elif defined(HAVE_X11)
  #include <unistd.h>
  #include "xgraph.h"
  #define  getch getchar
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>

#include "grafix3d.h"               // TriPlot3D class header file

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::TriPlot3D
// Scope:        Public, Constructor
TriPlot3D::TriPlot3D (UINT CTriangles, UINT CNodes,
                      AINODE *Aainode, PNT3D *Apnt3dInput,
                      double Az, double El,
                      int NResolution, int iSwitches,
                      int NumCon, int Iztick, UINT Perspective,
                      BOOL logdata)
{
  NumContour = NumCon;
  iztick     = Iztick;
  Persp      = Perspective;
  fScale     = iSwitches & F_SCALE;
  fGrayLevel = iSwitches & F_GRAYLEVEL;
  fColor     = iSwitches & F_COLOR;
  fWireframe = iSwitches & F_WIREFRAME;
  flogz      = logdata;
  
  // Get the View Aspect from the user
  GetViewAspect(Az, El);

  cMaxColor = GRSetLUT(fGrayLevel);    // the style of showing

  // Open a window
  OpenWindow(NResolution);

  cTriangles = CTriangles;             // number of triangles
  cNodes = CNodes;                     // number of nodes

  aainode     = Aainode;
  apnt3dInput = Apnt3dInput;
  apnt3dXform = new PNT3D[cNodes];
  apnt3dMesh  = new PNT3D[cNodes];
  apnt2d      = new PNT2D[cNodes];
  apnt2dMesh  = new PNT2D[cNodes];
  aiPlotOrder = new UINT[cTriangles];

  Xform3DData();                      // Transform the data 3D->3D
  Xform2DData();                      // Transform the data 3D->2D
  DeterminePlotOrder();               // Sort the order of appearance
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::~TriPlot3D
// Scope:        Public, Destructor
TriPlot3D::~TriPlot3D ()
{
  CloseWindow();

  delete [] apnt3dXform;
  delete [] apnt3dMesh;
  delete [] apnt2d;
  delete [] apnt2dMesh;
  delete [] aiPlotOrder;
}

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::GetViewAspect
// Scope:        Private
// Description:  This function uses the azimuthal and the elevation angle to
//               compute the elements of the transformation/rotation matrix.
//               After rotation, the data will be mapped to new yz plane,
//               and the x-component will determine the plot order.
//               Ref: Chap 3 of "Modern Quantum Mechanics" by Sakurai.
//                    which is using Euler angles.
void TriPlot3D::GetViewAspect (double & az, double & el)
{
  // try to convert AZ range (0 to 360), EL range (-90 to 90) in the future

  while(az<0)   az+=360;
  while(az>360) az-=360;
  while(el<0)   el=0;
  while(el>80)  el=80;
  AZ = az;
  EL = el;

  double AZrad = AZ * M_PI / 180.0;
  double ELrad = EL * M_PI / 180.0;

  double sinAZ = sin(AZrad);
  double cosAZ = cos(AZrad);
  double sinEL = sin(ELrad);
  double cosEL = cos(ELrad);

  Cxx =   cosAZ * cosEL;      Cxy =   sinAZ * cosEL;      Cxz =   sinEL;
  Cyx = - sinAZ;              Cyy =   cosAZ;              Cyz =   0.0;
  Czx = - cosAZ * sinEL;      Czy = - sinAZ * sinEL;      Czz =   cosEL;
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::Xform3DData
// Scope:        Private
// Description:  This function transforms the data.
void TriPlot3D::Xform3DData ( )
{
  // pnt3dMin.x = 1e300;   pnt3dMax.x = -1e300;
  // pnt3dMin.y = 1e300;   pnt3dMax.y = -1e300;
  // pnt3dMin.z = 1e300;   pnt3dMax.z = -1e300;

  pnt3dMin.x = pnt3dMax.x = apnt3dInput[0].x;
  pnt3dMin.y = pnt3dMax.y = apnt3dInput[0].y;
  pnt3dMin.z = pnt3dMax.z = apnt3dInput[0].z;

  UINT i;
  for (i = 1; i < cNodes; i++)
  { if (apnt3dInput[i].x < pnt3dMin.x) pnt3dMin.x = apnt3dInput[i].x;
    if (apnt3dInput[i].y < pnt3dMin.y) pnt3dMin.y = apnt3dInput[i].y;
    if (apnt3dInput[i].z < pnt3dMin.z) pnt3dMin.z = apnt3dInput[i].z;
    if (apnt3dInput[i].x > pnt3dMax.x) pnt3dMax.x = apnt3dInput[i].x;
    if (apnt3dInput[i].y > pnt3dMax.y) pnt3dMax.y = apnt3dInput[i].y;
    if (apnt3dInput[i].z > pnt3dMax.z) pnt3dMax.z = apnt3dInput[i].z;
  }

  if (pnt3dMax.z == pnt3dMin.z) pnt3dMax.z = 2.0 * (1.0 + fabs(pnt3dMin.z));
  pnt3dScale.x = pnt3dMax.x - pnt3dMin.x;
  pnt3dScale.y = pnt3dMax.y - pnt3dMin.y;
  pnt3dScale.z = pnt3dMax.z - pnt3dMin.z;

  for (i = 0; i < cNodes; ++i)
  { // Rotate the point and scale the point
    // Result will be in apnt3dXform, apnt3dInput unchanged for contour plot
    Xform3D_3D(apnt3dInput[i], apnt3dXform[i]);
    // Same as above. mesh -> x = 0, i.e plotted first
    XformMesh(apnt3dInput[i], apnt3dMesh[i]);
  }

  PNT3D pnt3dA, pnt3dB;
  pnt3dA.x = pnt3dMin.x;
  pnt3dA.y = pnt3dMin.y;
  pnt3dA.z = pnt3dMin.z;
  apnt3dCorner[0] = pnt3dA;
  Xform3D_3D(pnt3dA, pnt3dB);
  pnt3dMinXform.x = pnt3dMaxXform.x = pnt3dB.x;
  pnt3dMinXform.y = pnt3dMaxXform.y = pnt3dB.y;
  pnt3dMinXform.z = pnt3dMaxXform.z = pnt3dB.z;
  // 8 corners of the cube
  for (i = 1; i < 8; ++i)
  { pnt3dA.x = (i & 0x01) ? pnt3dMax.x : pnt3dMin.x;
    pnt3dA.y = (i & 0x02) ? pnt3dMax.y : pnt3dMin.y;
    pnt3dA.z = (i & 0x04) ? pnt3dMax.z : pnt3dMin.z;

    apnt3dCorner[i] = pnt3dA;

    Xform3D_3D(pnt3dA, pnt3dB);
    if (pnt3dB.x < pnt3dMinXform.x) pnt3dMinXform.x = pnt3dB.x;
    if (pnt3dB.y < pnt3dMinXform.y) pnt3dMinXform.y = pnt3dB.y;
    if (pnt3dB.z < pnt3dMinXform.z) pnt3dMinXform.z = pnt3dB.z;
    if (pnt3dB.x > pnt3dMaxXform.x) pnt3dMaxXform.x = pnt3dB.x;
    if (pnt3dB.y > pnt3dMaxXform.y) pnt3dMaxXform.y = pnt3dB.y;
    if (pnt3dB.z > pnt3dMaxXform.z) pnt3dMaxXform.z = pnt3dB.z;
  }
  pnt3dScaleXform.x = pnt3dMaxXform.x - pnt3dMinXform.x;
  pnt3dScaleXform.y = pnt3dMaxXform.y - pnt3dMinXform.y;
  pnt3dScaleXform.z = pnt3dMaxXform.z - pnt3dMinXform.z;

  if (!fScale)
  { if (pnt3dScaleXform.y > pnt3dScaleXform.z)
      pnt3dScaleXform.z = pnt3dScaleXform.y;
    else
      pnt3dScaleXform.y = pnt3dScaleXform.z;
  }

}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::Xform2DData
// Scope:        Public
// Description:  This function transforms the data.
void TriPlot3D::Xform2DData ( )
{
  for(UINT i = 0; i < cNodes; ++i)
  { Xform3D_2D(apnt3dXform[i], apnt2d[i]);
    Xform3D_2D(apnt3dMesh[i], apnt2dMesh[i]);
  }
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::Xform3D_3D
// Scope:        Private
// Description:  This function peforms the 3D rotation ... note the x
//               component does not take into account the orignal z value.
void TriPlot3D::Xform3D_3D ( PNT3D &pnt3dIn, PNT3D &pnt3dOut )
{
  double x = (pnt3dIn.x - pnt3dMin.x) / pnt3dScale.x;
  double y = (pnt3dIn.y - pnt3dMin.y) / pnt3dScale.y;
  double z = (pnt3dIn.z - pnt3dMin.z) / pnt3dScale.z;

  pnt3dOut.x = Cxx * x + Cxy * y + Cxz * z;
  pnt3dOut.y = Cyx * x + Cyy * y + Cyz * z;
  pnt3dOut.z = Czx * x + Czy * y + Czz * z;
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::Xform3D_2D
// Scope:        Private
// Description:  This function peforms the 3D to 2D transformation
//               with perspective ratio.
void TriPlot3D::Xform3D_2D (PNT3D &pnt3dIn, PNT2D &pnt2dOut)
{
  double x = (pnt3dIn.y - pnt3dMinXform.y) / pnt3dScaleXform.y;
  double y = (pnt3dIn.z - pnt3dMinXform.z) / pnt3dScaleXform.z;

  if (Persp)
  { double ratio = (double) Persp + 3.0;
    x -= .5;  x *= (1 + pnt3dIn.x / (pnt3dScaleXform.x*ratio - pnt3dIn.x));
    y -= .5;  y *= (1 + pnt3dIn.x / (pnt3dScaleXform.x*ratio - pnt3dIn.x));
    x += .5;
    y += .5;
  }

  pnt2dOut.x = (int) ((double) cx * x) + cxLBorder;
  pnt2dOut.y = (int) ((double) cy * (1.0 - y)) + cyLBorder;
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::XformMesh
// Scope:        Private
// Description:  This function peforms the 3D to 3D transformation for the
//               mesh triangles.
void TriPlot3D::XformMesh (PNT3D &pnt3dIn, PNT3D &pnt3dOut)
{
  double x = (pnt3dIn.x - pnt3dMin.x) / pnt3dScale.x;
  double y = (pnt3dIn.y - pnt3dMin.y) / pnt3dScale.y;

  // They are ploted first
  pnt3dOut.x =   Cxx * x + Cxy * y + (EL < 0 ? Cxz : 0.);
  pnt3dOut.y =   Cyx * x + Cyy * y + (EL < 0 ? Cyz : 0.);
  pnt3dOut.z =   Czx * x + Czy * y + (EL < 0 ? Czz : 0.);
}

///////////////////////////////////////////////////////////////////////////////
// Function:     QuickSortCompare
// Scope:        Private
// Description:  This function compares two surface elements in order to
//               determine which should be plotted first.  This function is
//               called by qsort.
//               The triangle will be plotted according to their original
//               coordinates in x-y plane projections.
static AINODE *aainode;
static PNT3D  *apnt3d;
static double Nx, Ny;
static int QuickSortCompare (const void *p1, const void *p2)
{
  UINT i1 = *((UINT *) p1);
  UINT i2 = *((UINT *) p2);
  UINT i1A = aainode[i1][0];
  UINT i1B = aainode[i1][1];
  UINT i1C = aainode[i1][2];
  UINT i1D = aainode[i1][3];
  UINT i2A = aainode[i2][0];
  UINT i2B = aainode[i2][1];
  UINT i2C = aainode[i2][2];
  UINT i2D = aainode[i2][3];
  double x1A = Nx * apnt3d[i1A].x + Ny * apnt3d[i1A].y;
  double x1B = Nx * apnt3d[i1B].x + Ny * apnt3d[i1B].y;
  double x1C = Nx * apnt3d[i1C].x + Ny * apnt3d[i1C].y;
  double x2A = Nx * apnt3d[i2A].x + Ny * apnt3d[i2A].y;
  double x2B = Nx * apnt3d[i2B].x + Ny * apnt3d[i2B].y;
  double x2C = Nx * apnt3d[i2C].x + Ny * apnt3d[i2C].y;
  double m1  = (x1A > x1B) ? x1A : x1B;
  double m2  = (x2A > x2B) ? x2A : x2B;
  m1  = (m1 > x1C) ? m1 : x1C;
  m2  = (m2 > x2C) ? m2 : x2C;
  if (i1D != -1u)
  { double x1D = Nx * apnt3d[i1D].x + Ny * apnt3d[i1D].y;
    m1 = (m1 > x1D) ? m1 : x1D;
  }
  if (i2D != -1u)
  { double x2D = Nx * apnt3d[i2D].x + Ny * apnt3d[i2D].y;
    m2 = (m2 > x2D) ? m2 : x2D;
  }
  // conform to the usage of qsort: ascending order
  return (m1 > m2) ? 1 : -1;
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::DeterminePlotOrder
// Scope:        Private
// Description:  This function determines the plot order.
void TriPlot3D::DeterminePlotOrder ()
{
  for(UINT i = 0; i < cTriangles; ++i) aiPlotOrder[i] = i;
  ::aainode = aainode;
  ::apnt3d  = apnt3dInput;
  // ::apnt3d  = apnt3dXform;
  ::Nx = cos(AZ * M_PI / 180.0);
  ::Ny = sin(AZ * M_PI / 180.0);
  qsort(aiPlotOrder, cTriangles, sizeof(UINT), QuickSortCompare);
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::OpenWindow
// Scope:        Private
// Description:  This function opens the graphics window
void TriPlot3D::OpenWindow (int NResolution)
{

  GRInitGraphics (CHARPIXELWIDE, CHARPIXELHIGH);

  nResFactor = 1;
  cxLBorder = 70;
  cyLBorder = 30;
  cxRBorder = 55;
  cyRBorder = 45;


  if      (NResolution == RES_640x480)
  { cx = 640 - cxRBorder - cxLBorder;
    cy = 480 - cyRBorder - cyLBorder;
  }
  else if (NResolution == RES_1024x768)
  { cx = 1024 - cxRBorder - cxLBorder;
    cy = 768  - cyRBorder - cyLBorder;
  }
  else
  { cx = 800 - cxRBorder - cxLBorder;
    cy = 600 - cyRBorder - cyLBorder;
  }


}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::CloseWindow
// Scope:        Public
// Description:  This function closes the graphics window
void TriPlot3D::CloseWindow ()
{
  // closegraph();
  //GRFreeGraphics();
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::PlotMesh
// Scope:        Public
// Description:  This function plots the mesh
void TriPlot3D::PlotMesh ()
{
  GRSetColor(GR_WHITE);
  for(UINT i = 0; i < cTriangles; ++i)
  { UINT iA = aainode[i][0];
    UINT iB = aainode[i][1];
    UINT iC = aainode[i][2];
    UINT iD = aainode[i][3];

    int an[10];
    an[0] = apnt2dMesh[iA].x;
    an[1] = apnt2dMesh[iA].y;
    an[2] = apnt2dMesh[iB].x;
    an[3] = apnt2dMesh[iB].y;
    an[4] = apnt2dMesh[iC].x;
    an[5] = apnt2dMesh[iC].y;

    if (iD == -1u)
    { an[6] = an[0];
      an[7] = an[1];
      GRDrawPoly(4, an);
    }
    else
    { an[6] = apnt2dMesh[iD].x;
      an[7] = apnt2dMesh[iD].y;
      an[8] = an[0];
      an[9] = an[1];
      GRDrawPoly(5, an);
    }

  }
}


///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::PlotData
// Scope:        Public
// Description:  This function plots the data
void TriPlot3D::PlotData ( BOOL fShowSurf, BOOL fShowCont )
{
  int    an[8];        // Array to hold the points co-ordinates (data)
  PNT3D  apn[4];       // Array to hold the points co-ordinates (contour)
  double   z;
  UINT   c;
  int    color = 0;

  GRSetColor(GR_WHITE);
  GRSetFillColor(GR_RED, PS_SOLID);

  //Find out the length between each contour level.
  double contLen = (pnt3dScale.z) / (NumContour + 1);

  //Loop through all the triangles/rectangles
  for(UINT i = 0; i < cTriangles; ++i)
  { UINT j = aiPlotOrder[i];

    //Use the following will make the graph not shown in order
    //UINT j = i;

    //Obtain the index of which node to use.
    UINT iA = aainode[j][0];
    UINT iB = aainode[j][1];
    UINT iC = aainode[j][2];
    UINT iD = aainode[j][3];

    if (fShowSurf)
    { an[0] = apnt2d[iA].x;
      an[1] = apnt2d[iA].y;
      an[2] = apnt2d[iB].x;
      an[3] = apnt2d[iB].y;
      an[4] = apnt2d[iC].x;
      an[5] = apnt2d[iC].y;

      if (iD == -1u)
        c = 3;
      else
      { an[6] = apnt2d[iD].x;
        an[7] = apnt2d[iD].y;
        c = 4;
      }

      if (!fWireframe)
      { if (iD == -1u)
        { z = apnt3dInput[iA].z + apnt3dInput[iB].z + apnt3dInput[iC].z;
          z /= 3.0;
        }
        else
        { z = apnt3dInput[iA].z + apnt3dInput[iB].z +
              apnt3dInput[iC].z + apnt3dInput[iD].z;
          z /= 4.0;
        }

        // -------------- set the fillcolor by the hight -----------------
        z -= pnt3dMin.z;
        z /= pnt3dScale.z;
        if (fGrayLevel)
          color = (int) (z * cMaxColor);

        else if (fColor)
        {
                if(z>=1.0)
                        color = cMaxColor+Rainbow_Color_NUM;
                else
                        color = cMaxColor+1+(int)(z*Rainbow_Color_NUM);
        }
        else
        { if      (z >= 2.0/3.0) color = GR_BLACK;
          else if (z >= 1.0/3.0) color = GR_LIGHTGRAY;
          else                   color = GR_DARKGRAY;
        }
        GRSetFillColor(color, PS_SOLID);
    }
      GRFillPoly(c, an);
    }

    if (fShowCont)
    { // ------------------ CONTOUR PLOT -------------------------------
      // ---------- Obtain the triangle/rectangle points ---------------
      apn[0] = apnt3dInput[iA];
      apn[1] = apnt3dInput[iB];
      apn[2] = apnt3dInput[iC];

      GRSetColor(GR_WHITE);

      //Test if it is a triangle OR a rectangle
      if (iD != -1u)
      { // A rectangle, add one more point.
        apn[3] = apnt3dInput[iD];
        // Separate the rectangle to 2 triangles
        // *** Assume the following:
        // (0)--(3)
        //  | \  |
        //  |  \ |
        // (1)--(2)
        // The first triangle we will draw it later. So
        // just pass the second triangle
        DrawContour(apn[0], apn[2], apn[3], contLen, fShowSurf && fShowCont);
      }
      //Else, a triangle, draw it.
      DrawContour(apn[0], apn[1], apn[2], contLen, fShowSurf && fShowCont);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////  Frame/Grid probably by Neil KH  ////////////////////////
///////////////////////////////////////////////////////////////////////////////
//                         z
//                         |
//
//                         +------------+
//                        /|4          /|6
//                       / |          / |
//                      /  |         /  |
//                     +------------+   |
//                     |5  |        |7  |
//                     |   +------- |---+  -- y
//                     |  / 0       |  / 2
//                     | /          | /
//                     |/           |/
//                     +------------+
//                      1            3
//                   /
//                  x

// table for selection of corner to be removed depending on rotation angles
// and for drawing frame

UINT aauFrame[8][3] = {
// MAX          xy       xz       yz           // zyx
// ---        ------   ------   ------         // ---
/*  0  */  {  0x7645,  0x7623,  0x7513  },     // 000
/*  1  */  {  0x6764,  0x6732,  0x6457  },     // 001
/*  2  */  {  0x5467,  0x5401,  0x5731  },     // 010
/*  3  */  {  0x4576,  0x4510,  0x4620  },     // 011
/*  4  */  {  0x3201,  0x3267,  0x3157  },     // 100
/*  5  */  {  0x2310,  0x2376,  0x2046  },     // 101
/*  6  */  {  0x1023,  0x1045,  0x1375  },     // 110
/*  7  */  {  0x0132,  0x0154,  0x0264  }      // 111
};

// **** Possible CRASH if elevation not between 0 and 90 deg. ****
// ****   due to incomplete initialization of array           ****
//
// table for label line min and max corners
UINT labelCorners[8][6] = {
// iCorner  xLmax  xLmin  yLmax  yLmin  zLmax  zLmin
// -------  -----  -----  -----  -----  -----  -----
/*  0  */  { 5,     4,     6,     4,     6,     2 },
/*  1  */  { 5,     4,     7,     5,     4,     0 },
/*  2  */  { 7,     6,     6,     4,     7,     3 },
/*  3  */  { 7,     6,     7,     5,     5,     1 },
/*  4  */  { 1,     0,     2,     0,     6,     2 },
/*  5  */  { 1,     0,     3,     1,     4,     0 },
/*  6  */  { 3,     2,     2,     0,     7,     3 },
/*  7  */  { 3,     2,     3,     1,     5,     1 }
};

// table for label stuff which is azimuth quadrant dependant
// for x,y,z grow, -1 means grow negative, 1 means grow positive
// for z source of offsets, 0 means use x, 1 means use y
// z source is also used for label justification
//   0 means x is left justified and y is right justified
//   1 means y is left justified and x is right justified
int labelDir[8][5] = {
// iCorner  xGrow  yGrow          zyGrow  zSource
// -------  -----  -----  ------  ------  -------
/*  0  */  {-1,    -1,     1,     1,      0 },
/*  1  */  {-1,     1,     0,    -1,      1 },
/*  2  */  { 1,    -1,     0,     1,      1 },
/*  3  */  { 1,     1,     1,    -1,      0 },
/*  4  */  {-1,    -1,     1,     1,      0 },
/*  5  */  {-1,     1,     0,    -1,      1 },
/*  6  */  { 1,    -1,     0,     1,      1 },
/*  7  */  { 1,     1,     1,    -1,      0 }
};

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::PlotGrid             (N. K. Habermehl, J.Josh Feng)
// Scope:        Public
// Description:  This function plots the grid and labels.
void TriPlot3D::PlotGrid ()
{
  UINT ActualLines;
  UINT iCorner;                       // the missing corner in the background

  // if structure eliminates corner, break up azimuthal rotation into 4 zones
  if      (AZ <=  90.0 && AZ >=  0.0) iCorner = (EL>=0) ? 7 : 3;
  else if (AZ <= 180.0 && AZ >  90.0) iCorner = (EL>=0) ? 6 : 2;
  else if (AZ <= 270.0 && AZ > 180.0) iCorner = (EL>=0) ? 4 : 0; // 2D included
  else if (AZ <= 360.0 && AZ > 270.0) iCorner = (EL>=0) ? 5 : 1;

  DrawFrame(iCorner);
  if ((AZ == 270) && (EL == 90))
    ActualLines = 0;                  // 2D case
  else
    ActualLines = DrawGridLines(iCorner, iztick);
  DrawFrameLabels(iCorner, ActualLines);
}

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::DrawFrame                   (N. K. Habermehl)
// Scope:        Public
// Description:  This function ....
void TriPlot3D::DrawFrame (UINT iCorner)
{
#if DEBUG
        printf("in TriPlot3D::DrawFrame\n");
#endif

  for (UINT i = 0; i < 3; i++)
  { UINT u = aauFrame[iCorner][i];
    UINT i1 =  u & 0x000F;
    UINT i2 = (u & 0x00F0) >> 4;
    UINT i3 = (u & 0x0F00) >> 8;
    UINT i4 = (u & 0xF000) >> 12;

    line3D(apnt3dCorner[i1], apnt3dCorner[i2]);
    line3D(apnt3dCorner[i2], apnt3dCorner[i3]);
    line3D(apnt3dCorner[i3], apnt3dCorner[i4]);
    line3D(apnt3dCorner[i4], apnt3dCorner[i1]);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::DrawGridLines               (N. K. Habermehl)
// Scope:        Public
// Description:  This function draws the number of gridlines specified or
//                      automatically calculates the number of grid lines if
//                      a negative number was passed.
// Return:       The number of grid lines actually drawn, not a negative number.
UINT TriPlot3D::DrawGridLines (UINT iCorner, int cNumOfGridLine)
{
  UINT xz = aauFrame[iCorner][1];
  UINT yz = aauFrame[iCorner][2];

  // ------------- draw the gridlines on the xz plane ------------------
  PNT3D pnt3dB1 = apnt3dCorner[(xz & 0xF000) >> 12];
  PNT3D pnt3dT1 = apnt3dCorner[(xz & 0x000F)      ];
  PNT3D pnt3dB2 = apnt3dCorner[(xz & 0x0F00) >> 8 ];
  /* PNT3D pnt3dT2 = apnt3dCorner[(xz & 0x00F0) >> 4 ]; */

  if (cNumOfGridLine < 0)       // neg number, automatic calculation of
  {                             // number of gridlines
    PNT3D pnt3dB1Rot, pnt3dT1Rot;
    PNT2D pnt2dB1Rot, pnt2dT1Rot;
    Xform3D_3D(pnt3dB1, pnt3dB1Rot);
    Xform3D_3D(pnt3dT1, pnt3dT1Rot);
    Xform3D_2D(pnt3dB1Rot, pnt2dB1Rot);
    Xform3D_2D(pnt3dT1Rot, pnt2dT1Rot);
    UINT ZPixelAvail = abs(pnt2dT1Rot.y - pnt2dB1Rot.y);
    cNumOfGridLine = (ZPixelAvail / MIN_PIXEL_PER_LINE / nResFactor) - 1;
  }

  UINT i;
  PNT3D pnt3d1, pnt3d2;
  pnt3d1.x = pnt3dB1.x;
  pnt3d1.y = pnt3dB1.y;

  pnt3d2.x = pnt3dB2.x;
  pnt3d2.y = pnt3dB2.y;
  double dz = (pnt3dT1.z - pnt3dB1.z) / (double) (cNumOfGridLine + 1);
  for (i = 1; i <= cNumOfGridLine; ++i)
  { double z = pnt3dB1.z + (double) i * dz;
    pnt3d1.z = z;
    pnt3d2.z = z;
    line3D(pnt3d1, pnt3d2);
  }

  // ------------------ draw the gridlines on the yz plane --------------
  pnt3dB2 = apnt3dCorner[(yz & 0x0F00) >> 8 ];
  pnt3d2.x = pnt3dB2.x;
  pnt3d2.y = pnt3dB2.y;
  for (i = 1; i <= cNumOfGridLine; ++i)
  { double z = pnt3dB1.z + (double) i * dz;
    pnt3d1.z = z;
    pnt3d2.z = z;
    line3D(pnt3d1, pnt3d2);
  }

  return cNumOfGridLine;
}

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::DrawFrameLabels      (N. K. Habermehl, J.Josh Feng)
// Scope:        Public
// Description:  This function draws tickmarks, and scale labels.
// Parameters:   iCorner - the corner which was eliminated,
//                      in effect specifies quadrant perspective
//               ActualLines - the actual number of z lines that were drawn
//               xTicks, yTicks - the number of tick marks that should be
//                      drawn, default -1 is sentinal value for automatic
//                      generation, no check for too many lines, could lead
//                      to bunched up lables.
//
void TriPlot3D::DrawFrameLabels(UINT iCorner, UINT ActualLines,
                                int xTicks, int yTicks)
{
  char labelBuff[30]; // buffer to hold label strings for outtextxy()
  double mulFac;        // scaling factor for all axis labels = 10 to the tenExp
  int  tenExp;        // exponent of 10 to be displayed
  UINT mid;           // half the number of tick marks per axis
  double fabsMax;       // maximum absolute value of 2 doubles

  // get label line max and min corners from table
  UINT xLmax = labelCorners[iCorner][0];
  UINT xLmin = labelCorners[iCorner][1];
  UINT yLmax = labelCorners[iCorner][2];
  UINT yLmin = labelCorners[iCorner][3];
  UINT zLmax = labelCorners[iCorner][4];
  UINT zLmin = labelCorners[iCorner][5];

  // convert to 3D points
  PNT3D pnt3dxLmax = apnt3dCorner[xLmax];
  PNT3D pnt3dxLmin = apnt3dCorner[xLmin];
  PNT3D pnt3dyLmax = apnt3dCorner[yLmax];
  PNT3D pnt3dyLmin = apnt3dCorner[yLmin];
  PNT3D pnt3dzLmax = apnt3dCorner[zLmax];
  PNT3D pnt3dzLmin = apnt3dCorner[zLmin];

  // automatic tick line setup for x axis
  PNT3D pnt3dxLmaxRot, pnt3dxLminRot;
  PNT2D pnt2dxLmaxRot, pnt2dxLminRot;
  Xform3D_3D(pnt3dxLmax, pnt3dxLmaxRot);
  Xform3D_3D(pnt3dxLmin, pnt3dxLminRot);
  Xform3D_2D(pnt3dxLmaxRot, pnt2dxLmaxRot);
  Xform3D_2D(pnt3dxLminRot, pnt2dxLminRot);
  UINT XPixelAvail = (UINT) sqrt(               // Pathag theorem
   (pnt2dxLmaxRot.y - pnt2dxLminRot.y) * (pnt2dxLmaxRot.y - pnt2dxLminRot.y)
                                       +
   (pnt2dxLmaxRot.x - pnt2dxLminRot.x) * (pnt2dxLmaxRot.x - pnt2dxLminRot.x));
  UINT xNumOfTickMarks = (XPixelAvail / MIN_PIXEL_PER_LINE / nResFactor) + 1;
  double yLineOffset = (pnt3dxLmax.x - pnt3dxLmin.x) /
                     (xNumOfTickMarks * LABELLINEFACTOR);

  // automatic tick line setup for y axis
  PNT3D pnt3dyLmaxRot, pnt3dyLminRot;
  PNT2D pnt2dyLmaxRot, pnt2dyLminRot;
  Xform3D_3D(pnt3dyLmax, pnt3dyLmaxRot);
  Xform3D_3D(pnt3dyLmin, pnt3dyLminRot);
  Xform3D_2D(pnt3dyLmaxRot, pnt2dyLmaxRot);
  Xform3D_2D(pnt3dyLminRot, pnt2dyLminRot);
  UINT YPixelAvail = (UINT) sqrt(               // Pathag theorem
   (pnt2dyLmaxRot.y - pnt2dyLminRot.y) * (pnt2dyLmaxRot.y - pnt2dyLminRot.y)
                                       +
   (pnt2dyLmaxRot.x - pnt2dyLminRot.x) * (pnt2dyLmaxRot.x - pnt2dyLminRot.x));
  UINT yNumOfTickMarks = (YPixelAvail / MIN_PIXEL_PER_LINE / nResFactor) + 1;
  double xLineOffset = (pnt3dyLmax.y - pnt3dyLmin.y) /
                     (yNumOfTickMarks * LABELLINEFACTOR);

  double zLineOffset;

  // automatic or manual z lines were done in DrawGridLines() and passed
  // here in ActualLines
  UINT zNumOfTickMarks = ActualLines + 2;

  double xTickLength = xLineOffset / TICKFACTOR;
  double yTickLength = yLineOffset / TICKFACTOR;
  double zTickLength;                             // set in each case
  PNT3D tickOnLine, tickOff, labelOff;          // temp tick line points
  if (xTicks >= 0) xNumOfTickMarks = xTicks;    // use parameters here so
  if (yTicks >= 0) yNumOfTickMarks = yTicks;    // auto generation sets
                                                // offsets

  // draw the label lines.  I used a case structure because the label lines,
  // tick marks, and labels do not stay in one place relative to the data
  // as the data is rotated.  Rather, the labels must "ratchet around" the
  // data in order to maintain a "constant" 2D position.  I intend to
  // place each label as I place each tick mark, so I am starting out by
  // placing the tick marks in acending numerical order.  This is another
  // reason to adopt a case structure, as up/down, left/right, and
  // greater/lesser occure in various combinations throughout the
  // rotation.

  if (labelDir[iCorner][4] == 1)
  { zTickLength = yTickLength;
    zLineOffset = yLineOffset;
  }
  else
  { zTickLength = xTickLength;
    zLineOffset = xLineOffset;
  }

  // --------------- x tick and label stuff -----------------------
  fabsMax = fabs(pnt3dMax.x);
  if (fabs(pnt3dMin.x) > fabsMax) fabsMax = fabs(pnt3dMin.x);
  tenExp = int(log10(fabsMax)); // test Max and scale exponent
  mulFac = 1 / pow(10, tenExp);
  tickOnLine = pnt3dxLmin;
  tickOff    = pnt3dxLmin;
  tickOff.y += xTickLength * labelDir[iCorner][0];
  labelOff   = tickOff;
  labelOff.y += (xLineOffset / CHARFACTOR) * labelDir[iCorner][0];
  double xTickInc = (pnt3dxLmax.x - pnt3dxLmin.x) / (xNumOfTickMarks - 1);
  mid = xNumOfTickMarks / 2;
  for (UINT i = 0; i < xNumOfTickMarks; i++)
  { line3D(tickOnLine, tickOff);        // draw x tick lines
    sprintf(labelBuff, "%3.2f", mulFac * (pnt3dxLmin.x + i*xTickInc));
    // write label text to screen
    text3D(labelOff, labelBuff, labelDir[iCorner][4], (EL > 0));
    if (i == mid)
      exponentLabel(labelOff, tenExp, labelDir[iCorner][4], (EL > 0));
    tickOnLine.x += xTickInc;
    tickOff.x    += xTickInc;
    labelOff.x   += xTickInc;
  }

  // ---------------------- y tick and label stuff --------------------
  fabsMax = fabs(pnt3dMax.y);
  if (fabs(pnt3dMin.y) > fabsMax) fabsMax = fabs(pnt3dMin.y);
  tenExp = int(log10(fabsMax)); // test Max and scale exponent
  mulFac = 1 / pow(10, tenExp);
  tickOnLine = pnt3dyLmin;
  tickOff    = pnt3dyLmin;
  tickOff.x += yTickLength * labelDir[iCorner][1];
  labelOff   = tickOff;
  labelOff.x += (yLineOffset / CHARFACTOR) * labelDir[iCorner][1];
  double yTickInc = (pnt3dyLmax.y - pnt3dyLmin.y) / (yNumOfTickMarks - 1);
  mid = yNumOfTickMarks / 2;
  for (UINT j = 0; j < yNumOfTickMarks; j++)
  { line3D(tickOnLine, tickOff);        // draw y tick lines
    sprintf(labelBuff, "%3.2f", mulFac * (pnt3dyLmin.y + j*yTickInc));
    // write label text to screen
    text3D(labelOff, labelBuff, labelDir[iCorner][2], (EL > 0));
    if (j == mid)
      exponentLabel(labelOff, tenExp, labelDir[iCorner][2], (EL > 0));
    tickOnLine.y += yTickInc;
    tickOff.y    += yTickInc;
    labelOff.y   += yTickInc;
  }
  
    
  if ((AZ != 270) || (EL != 90))   // NOT 2D case
  {
    // -------------------- z tick and label stuff ----------------------
    fabsMax = fabs(pnt3dMax.z);
    if (fabs(pnt3dMin.z) > fabsMax) fabsMax = fabs(pnt3dMin.z);
    if(!flogz)
    {
      tenExp = int(log10(fabsMax)); // test Max and scale exponent
      mulFac = 1 / pow(10, tenExp);
    }
    else
    {
      tenExp = 0;
      mulFac = 1;
    }  
    tickOnLine = pnt3dzLmin;
    tickOff    = pnt3dzLmin;
    if (labelDir[iCorner][4] == 1)
    { tickOff.x += zTickLength * labelDir[iCorner][3];
      labelOff   = tickOff;
      labelOff.x += (zLineOffset / CHARFACTOR) * labelDir[iCorner][3];
    }
    else
    { tickOff.y += zTickLength * labelDir[iCorner][3];
      labelOff   = tickOff;
      labelOff.y += (zLineOffset / CHARFACTOR) * labelDir[iCorner][3];
    }
    double zTickInc = (pnt3dzLmax.z - pnt3dzLmin.z) / (zNumOfTickMarks - 1);
    for (UINT k = 0; k < zNumOfTickMarks; k++)
    { line3D(tickOnLine, tickOff);        // draw z tick lines
      if(!flogz)
        sprintf(labelBuff, "%5.4f", mulFac * (pnt3dzLmin.z + k*zTickInc));
      else
      {
        double value = pnt3dzLmin.z + k*zTickInc;
        if(value>0)
         sprintf(labelBuff, "%3.2e", pow(10,value)-1.0);
        else
         sprintf(labelBuff, "%3.2e", 1.0-pow(10,-value)); 
      }
      text3D(labelOff, labelBuff, RIGHTJUST, (EL < 0));  // write label text to screen
      tickOnLine.z += zTickInc;
      tickOff.z    += zTickInc;
      labelOff.z   += zTickInc;
    }
    if (EL < 0) labelOff.z = pnt3dzLmin.z;
    exponentLabel(labelOff, tenExp, LEFTJUST, (EL > 0), ZLABEL);
  }
}

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::text3D               (N. K. Habermehl, J.Josh Feng)
// Scope:        Public
// Description:  This function plots a text string at a 3D point.
// Parameters:   pnt3dT - reference point to draw char string
//               *p_textstring - pointer to string buffer
//               justify - 1 if right justify
//                       - 0 if left justify
//                       - default to left justify
//               origin  - 1 if char box origin is upper left
//                       - 0 if char box origin is lower left
//                       - default to upper left origin
void TriPlot3D::text3D (PNT3D pnt3dT, char *p_textstring,
                         int justify, int origin)
{
  // 3D to 3D transformation
  PNT3D pnt3dTout;
  Xform3D_3D(pnt3dT, pnt3dTout);

  // 3D to 2D transformation
  PNT2D pnt2dT;
  Xform3D_2D(pnt3dTout, pnt2dT);

  // output the text
  if (justify == RIGHTJUST) pnt2dT.x -= nResFactor * CHARPIXELWIDE * strlen(p_textstring);
  if (origin == DOWNORIGIN) pnt2dT.y -= nResFactor * CHARPIXELHIGH;

  GRText(pnt2dT.x, pnt2dT.y, p_textstring);
}

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::exponentLabel        (N. K. Habermehl, J.Josh Feng)
// Scope:        Public
// Description:  This function places the exponent label near for x and y
//                 axes according to RIGHTJUST or LEFTJUST.  When ZLABEL
//                 is present (true) the label is placed directly above.
//
// Parameters:   start   - 3D starting point for the label
//               tenExp  - power of ten to be displayed
//               justify - LEFTJUST, start at point and work right
//                         RIGHTJUST, move left of point and work right again
//               dooZ    - if present (true) ignore justify and put label above
//                         start, defaults to false if not present
//
void TriPlot3D::exponentLabel(PNT3D start, int tenExp,
                              int justify, int Vjustify, BOOL dooZ)
{
  char p_exp[10];

  sprintf(p_exp, "%d", tenExp);

  // 3D to 3D transformation
  PNT3D startout3D;
  Xform3D_3D(start, startout3D);

  // 3D to 2D transformation
  PNT2D sout2D;
  Xform3D_2D(startout3D, sout2D);

  Vjustify = Vjustify ? -1 : 1;
  if (dooZ == ZLABEL)
  { sout2D.y += (ZLABELUP * nResFactor * CHARPIXELHIGH * Vjustify);
    if(!flogz)
    {
      GRText(3 * nResFactor * CHARPIXELWIDE, sout2D.y, "Z x10");
      GRText(8 * nResFactor * CHARPIXELWIDE, sout2D.y - int (.5 * nResFactor * CHARPIXELHIGH),p_exp);
    }
    else
      GRText(3 * nResFactor * CHARPIXELWIDE, sout2D.y, "Z log");
  }
  else if (justify == RIGHTJUST)
  { sout2D.y -= 2 * nResFactor * CHARPIXELHIGH * Vjustify;
    GRText(sout2D.x - (AXISRJCHARS * nResFactor * CHARPIXELWIDE),
              sout2D.y, "x10 um");
    GRText(sout2D.x - (AXISRJCHARS * nResFactor * CHARPIXELWIDE) + (3 * nResFactor * CHARPIXELWIDE),
              sout2D.y - int (0.5 * nResFactor * CHARPIXELHIGH), p_exp);
  }
  else
  { sout2D.y -= 2 * nResFactor * CHARPIXELHIGH * Vjustify;
    GRText(sout2D.x + (AXISLJCHARS * nResFactor * CHARPIXELWIDE),
              sout2D.y, "x10 um");
    GRText(sout2D.x + (AXISLJCHARS * nResFactor * CHARPIXELWIDE) + (3 * nResFactor * CHARPIXELWIDE),
              sout2D.y - int (0.5 * nResFactor * CHARPIXELHIGH), p_exp);
  }
}


///////////////////////////////////////////////////////////////////////////////
//////////////////////  Coutour Plot stuff  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Function: swap
// Scope: Public
// Description:  This function swap two variables with same data structure
void swap(PNT3D *val1,  PNT3D *val2)
{
  PNT3D temp = *val1;
  *val1 = *val2;
  *val2 = temp;
}

///////////////////////////////////////////////////////////////////////////////
// Function:    DrawContour
// Scope:       private
// Description: This function calculate the contour lines
void TriPlot3D::DrawContour ( PNT3D a, PNT3D b, PNT3D c,
                              double contLen, BOOL fConSurf )
{
#if DEBUG
        printf("in DrawContour\n");
#endif

  // a
  // |\
  // | \
  // b--c

  // Variable to hold which level is using now.
  double level;
  double factor;

  // Array to hold the points.
  PNT3D temp1, temp2;
  PNT3D tempa = a;
  PNT3D tempb = b;
  PNT3D tempc = c;

  // Sort the three points tempa, tempb, tempc in increasing z components;
  if (tempa.z > tempb.z) swap(&tempa, &tempb);
  if (tempb.z > tempc.z)
  { if(tempc.z < tempa.z) swap(&tempa, &tempc);
    swap(&tempb, &tempc);
  }
  //return if no grad.
  if(fabs(contLen)<1e-14) return;
  // Loop though all the levels
  for (level = pnt3dMin.z; level <= pnt3dMax.z; level += contLen)
  {
    // Check if the step size is in range: 1 < level < MaxContour
    //fprintf(stderr,"%e %e %e %e\n",level,tempa.z,tempb.z,tempc.z);
    if ((level == tempa.z) && (level == tempc.z))
    { line3D(tempa, tempb);
      line3D(tempb, tempc);
      line3D(tempc, tempa);
      break;
    }
    else if ((level >= tempa.z) && (level <= tempc.z))
    { // OK, test each edge, please.
      if      ((level == tempa.z) && (tempa.z == tempb.z) && (tempb.z != tempc.z))
      { temp1.x = tempa.x;    temp2.x = tempb.x;
        temp1.y = tempa.y;    temp2.y = tempb.y;
      }
      else if ((level == tempc.z) && (tempc.z == tempb.z) && (tempa.z != tempb.z))
      { temp1.x = tempb.x;    temp2.x = tempc.x;
        temp1.y = tempb.y;    temp2.y = tempc.y;
      }
      else if (level == tempb.z)
      { // Now tempb (level) MUST be between tempa and tempc
        // so a<level=b<c

        factor = (level - tempa.z) / (tempc.z - tempa.z);
        temp1.x = tempa.x + (tempc.x - tempa.x) * factor;
        temp1.y = tempa.y + (tempc.y - tempa.y) * factor;
        temp2.x = tempb.x;
        temp2.y = tempb.y;
      }
      else if (level < tempb.z)
      { // tempb and tempc above the tempa
        // a<level<b<c

        factor = (level - tempa.z) / (tempb.z - tempa.z);
        temp1.x = tempa.x + (tempb.x - tempa.x) * factor;
        temp1.y = tempa.y + (tempb.y - tempa.y) * factor;

        factor = (level - tempa.z) / (tempc.z - tempa.z);
        temp2.x = tempa.x + (tempc.x - tempa.x) * factor;
        temp2.y = tempa.y + (tempc.y - tempa.y) * factor;

      }
      else
      { // tempa and tempb are below tempc
        // a<b<level<c

        factor = (level - tempa.z) / (tempc.z - tempa.z);
        temp1.x = tempa.x + (tempc.x - tempa.x) * factor;
        temp1.y = tempa.y + (tempc.y - tempa.y) * factor;

        factor = (level - tempb.z) / (tempc.z - tempb.z);
        temp2.x = tempb.x + (tempc.x - tempb.x) * factor;
        temp2.y = tempb.y + (tempc.y - tempb.y) * factor;

      }
      temp1.z = (EL > 0) ? pnt3dMin.z : pnt3dMax.z;
      if (fConSurf) temp1.z = level;
      temp2.z = temp1.z;
      line3D(temp1, temp2);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
// Function:     TriPlot3D::line3D
// Scope:        Public
// Description:  This function plots a line between two 3D points.
void TriPlot3D::line3D (PNT3D pnt3dA, PNT3D pnt3dB)
{
  // 3D to 3D transformation
  PNT3D pnt3dAout, pnt3dBout;
  Xform3D_3D(pnt3dA, pnt3dAout);
  Xform3D_3D(pnt3dB, pnt3dBout);

  // 3D to 2D transformation
  PNT2D pnt2dA, pnt2dB;
  Xform3D_2D(pnt3dAout, pnt2dA);
  Xform3D_2D(pnt3dBout, pnt2dB);

  // draw the line
  GRLine(pnt2dA.x, pnt2dA.y, pnt2dB.x, pnt2dB.y);
}

#endif
