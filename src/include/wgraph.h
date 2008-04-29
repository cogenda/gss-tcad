/*-----------------------------------------------------------------------------

 FILE:     WGRAPH.H

 PROJECT:

 REVISION HISTORY
 DATE      INITIALS   DESCRIPTION
 --------  --------   --------------------------------------------------------
 04/15/97    KMK      initial implementation

 DESCRIPTION:

 This module contains the graphics functions.

-----------------------------------------------------------------------------*/
#ifdef HAVE_WIN32

#ifndef __wgraph_h
#define __wgraph_h

// definition of command line argument structure
#define RES_640x480         1       /* 640 x 480 resolution */
#define RES_800x600         2       /* 800 x 600 resolution */
#define RES_1024x768        3       /* 1024 x 768 resolution */

#define GR_BLACK        0
#define GR_BLUE         1
#define GR_GREEN        2
#define GR_CYAN         3
#define GR_RED          4
#define GR_MAGENTA      5
#define GR_BROWN        6
#define GR_LIGHTGRAY    7
#define GR_DARKGRAY     8
#define GR_LIGHTBLUE    9
#define GR_LIGHTGREEN   10
#define GR_LIGHTCYAN    11
#define GR_LIGHTRED     12
#define GR_LIGHTMAGENTA 13
#define GR_YELLOW       14
#define GR_WHITE        15

#define Rainbow_Color_NUM    20
#define Rainbow_0            16
#define Rainbow_1            17
#define Rainbow_2            18
#define Rainbow_3            19
#define Rainbow_4            20
#define Rainbow_5            21
#define Rainbow_6            22
#define Rainbow_7            23
#define Rainbow_8            24
#define Rainbow_9            25
#define Rainbow_10           26
#define Rainbow_11           27
#define Rainbow_12           28
#define Rainbow_13           29
#define Rainbow_14           30
#define Rainbow_15           31
#define Rainbow_16           32
#define Rainbow_17           33
#define Rainbow_18           34
#define Rainbow_19           35


typedef void (*PFREDRAW) ();              /* pfredraw */

int  GRPrintOpen      (const char *, int);
void GRPrintClose     (void);
void GRInitGraphics   (int = 7, int = 14);
void GRFreeGraphics   (void);
void GROpenGraphWin   (char *, char *, int, int, int, int, PFREDRAW);
void GRSetColor       (int);
void GRSetFillColor   (int, int=0);
void GRPutPixel       (int, int, int);
void GRLine           (int, int, int, int);
void GRCircle         (int, int, int);
void GRDrawPoly       (int, int *);
void GRFillPoly       (int, int *);
void GRText           (int, int, char *);
int  GRSetLUT         (int);

#endif

#endif
