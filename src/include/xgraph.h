/*-----------------------------------------------------------------------------

 FILE:     xgraph.h

 PROJECT:

 AUTHOR:   Kevin M. Kramer
           J. Joshua Feng

 REVISION HISTORY
 DATE      INITIALS   DESCRIPTION
 --------  --------   --------------------------------------------------------
 12/15/93    KMK      initial implementation
 02/04/97    JJF      Modification (Adding LUT functions)
 12/12/05    gdiso    Support rainbow color

 DESCRIPTION:

 This module contains functions implement a small subset of the Borland BGI
 functions under X windows. Hence this module should not be included when
 compiling under Borland C++ enviroments.

-----------------------------------------------------------------------------*/
#ifdef HAVE_X11

#ifndef __xgraph_h
#define __xgraph_h

// definition of command line argument structure
#define RES_640x480         1       /* 640 x 480 resolution */
#define RES_800x600         2       /* 800 x 600 resolution */
#define RES_1024x768        3       /* 1024 x 768 resolution */

#define GR_BLACK             0
#define GR_BLUE              1
#define GR_GREEN             2
#define GR_CYAN              3
#define GR_RED               4
#define GR_MAGENTA           5
#define GR_BROWN             6
#define GR_LIGHTGRAY         7
#define GR_DARKGRAY          8
#define GR_LIGHTBLUE         9
#define GR_LIGHTGREEN        10
#define GR_LIGHTCYAN         11
#define GR_LIGHTRED          12
#define GR_LIGHTMAGENTA      13
#define GR_YELLOW            14
#define GR_WHITE             15

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

#define DETECT            0

#define PS_SOLID          1

#define EXPOSE            0
#define CONFIGURENOTIFY   1
#define KEYPRESS          2
#define BUTTONPRESS       3
#define OTHEREVENT        4
#define BACKSPACE         5
#define ESCAPE            6
#define LEFTARROW         7
#define RIGHTARROW        8
#define UPARROW           9
#define DOWNARROW         10
#define BUTTONMOTION      11
#define BUTTONRELEASE     12

#define grOk              0
#define grNoInitGraph    -1
#define grNotDetected    -2
#define grFileNotFound   -3
#define grInvalidDriver  -4
#define grNoLoadMem      -5

typedef void  (*PFREDRAW) ();

int   GRPrintOpen     (const char *, int = 20);
void  GRPrintClose    (void);
void  GRInitGraphics  (int = 7, int = 14);
void  GRFreeGraphics  (void);
void  GROpenGraphWin  (char *, char *, int, int, int, int, PFREDRAW);
void  GRClearGraphWin (void);
void  GRSetColor      (int);
void  GRSetFillColor  (int, int = PS_SOLID);
void  GRPutPixel      (int, int, int);
void  GRLine          (int, int, int, int);
void  GRCircle        (int, int, int);
void  GRDrawPoly      (int, int *);
void  GRFillPoly      (int, int *);
void  GRText          (int, int, char *);
int   GRGetMaxX       (void);
int   GRGetMaxY       (void);
int   GRSetLUT        (int);
void  setcmdargs      (int, char **);
void  flushdisplay    (void);
int   getevent        (short int *, int *, int *, int *, int *);
int   GRSaveScreen    (const char *, int, int);
#endif /* __xgraph_h */

#endif

