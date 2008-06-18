/*-----------------------------------------------------------------------------

 FILE:     xgraph.cpp

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
#include "config.h"
#ifdef HAVE_X11

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/keysym.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xgraph.h"

typedef unsigned char     BYTE;         /* b   */
typedef unsigned char     UCHAR;        /* uch */
typedef int               BOOL;         /* f   */
typedef int               BOOLEAN;      /* f   */
typedef unsigned int      WORD;         /* w   */
typedef unsigned int      UINT;         /* u   */
typedef unsigned long     DWORD;        /* dw  */
typedef char *            PSTR;         /* psz */
typedef int  *            PINT;         /* pn  */
typedef int  *            PBOOL;        /* pf  */
typedef unsigned int *    PWORD;        /* pw  */
typedef long *            PLONG;        /* pl  */
typedef unsigned long *   PDWORD;       /* pdw */
typedef void *            PVOID;        /* pv  */
typedef double            REAL;         /*     */
typedef double *          PREAL;        /*     */


static int nScreen;
static int cDisplayWidth;
static int cDisplayHeight;
static int cXresolution;
static int cYresolution;
static int nXorigin;
static int nYorigin;
static int cColorDepth;
static int cMaxColor;
/* static int cStyle; */
static int nFillColor;
static int nCurrentColor;
static int cArg;
static int cTextWidth;
static int cTextHeight;
static UINT cBorderWidth;
static DWORD dwForeground;
static DWORD dwBackground;
static char *pszDisplayName;
static char **apszArg;
static GC **gc;
static Display *display;
static Window window;
static Visual *visual;
static XFontStruct *xfont;
static Colormap colormap;
static Pixmap pixmapIcon;
static XSizeHints sizehints;
static XPoint axpoint[128];  /* used by poly filled */
static XEvent xevent;

static FILE *filePS = NULL;

#define ICON_BITMAP_WIDTH  40
#define ICON_BITMAP_HEIGHT 40
static char szIconBitmap[] = /* SGFramework Icon */
  { 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0xff, 0xff, 0xff, 0xff, 0xff,
    0xff, 0xff, 0xff, 0xff, 0xff,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0xff, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00, 0x00
  };

static unsigned int rc[][3] =
  {
    0, 0, 65535,
    0, 29298, 65535,
    0, 41377, 65535,
    0, 50629, 65535,
    0, 58596, 65535,
    0, 65535, 65535,
    0, 65535, 56797,
    0, 65535, 46003,
    0, 65535, 32639,
    0, 65535, 0,
    32125, 65535, 0,
    46260, 65535, 0,
    56540, 65535, 0,
    65535, 65535, 0,
    65535, 59881, 0,
    65535, 53199, 0,
    65535, 45746, 0,
    65535, 38036, 0,
    65535, 26471, 0,
    65535, 0, 0
  };

void GetGC (GC *, char *);
void GetGCRainbow (GC *, unsigned short, unsigned short, unsigned short);
/*****************************************************************************/
/************************  POSTSCRIPT FUNCTIONS  *****************************/
/*****************************************************************************/


/******************************************************************************
 * Description:  This function opens the postscipt file.
 */
int GRPrintOpen ( const char *pszName, int  nResFactor )
{
  char szFileName[128];
  int i;

  sprintf(szFileName, "%s", pszName);
  filePS = fopen(szFileName, "wt");
  if(!filePS) return 1;

  fprintf(filePS,
          "%%!PS-Adobe-1.0\n"
          "%%%%Creator: SGFramework\n"
          "%%%%Title:  %s\n"
          "%%%%Pages: (atend)\n"
          "%%%%BoundingBox: 72 162 540 640 \n"
          "%%%%EndComments\n\n", pszName);
  fprintf(filePS,
          "%%! SGFramework prolog Version 1.1\n"
          "%%%%\n"
          "\n"
          "/setlc {/linecolor exch def} def\n"
          "/setfc {/fillcolor exch def} def\n"
          "\n"
          "/colors [\n");

  fprintf(filePS,
          "  [1.0 1.0 1.0] %% white\n"
          "  [0.0 0.0 1.0] %% blue\n"
          "  [0.0 1.0 0.0] %% green\n"
          "  [0.0 1.0 1.0] %% cyan\n"
          "  [1.0 0.0 0.0] %% red\n"
          "  [1.0 1.0 0.2] %% magenta\n"
          "  [0.5 0.4 0.3] %% brown\n"
          "  [0.7 0.7 0.7] %% light gray\n"
          "  [0.3 0.3 0.3] %% dark gray\n"
          "  [0.2 0.2 1.0] %% light blue\n"
          "  [0.2 1.0 0.2] %% light green\n"
          "  [0.5 1.0 1.0] %% light cyan\n"
          "  [1.0 0.2 0.2] %% light red\n"
          "  [1.0 1.0 0.5] %% light magenta\n"
          "  [1.0 1.0 0.0] %% yellow\n"
          "  [0.0 0.0 0.0] %% black\n");
  fprintf(filePS,
          "  [0.000000 0.000000 1.000000]\n"
          "  [0.000000 0.447059 1.000000]\n"
          "  [0.000000 0.631373 1.000000]\n"
          "  [0.000000 0.772549 1.000000]\n"
          "  [0.000000 0.894118 1.000000]\n"
          "  [0.000000 1.000000 1.000000]\n"
          "  [0.000000 1.000000 0.866667]\n"
          "  [0.000000 1.000000 0.701961]\n"
          "  [0.000000 1.000000 0.498039]\n"
          "  [0.000000 1.000000 0.000000]\n"
          "  [0.490196 1.000000 0.000000]\n"
          "  [0.705882 1.000000 0.000000]\n"
          "  [0.862745 1.000000 0.000000]\n"
          "  [1.000000 1.000000 0.000000]\n"
          "  [1.000000 0.913725 0.000000]\n"
          "  [1.000000 0.811765 0.000000]\n"
          "  [1.000000 0.698039 0.000000]\n"
          "  [1.000000 0.580392 0.000000]\n"
          "  [1.000000 0.403922 0.000000]\n"
          "  [1.000000 0.000000 0.000000]\n" );

  fprintf(filePS,
          "  [0.0 0.0 0.0] %% black\n"
          "] def\n"
          "\n"
          "/usecolor {colors exch get aload pop setrgbcolor} def\n"
          "\n"
          "/drawpoly\n"
          "{ 2 sub /points exch def \n"
          "  newpath moveto lineto \n"
          "  { points 0 ne {lineto /points points 1 sub def} {exit} ifelse\n"
          "  } loop\n"
          "  linecolor usecolor stroke\n"
          "} def\n"
          "\n"
          "/fillpoly\n"
          "{ 2 sub /points exch def\n"
          "  newpath moveto lineto\n"
          "  { points 0 ne {lineto /points points 1 sub def} {exit} ifelse\n"
          "  } loop\n"
          "  closepath gsave fillcolor usecolor fill grestore\n"
          "  linecolor usecolor stroke\n"
          "} def\n"
          "\n"
          "/line {newpath moveto lineto linecolor usecolor stroke} def\n"
          "\n"
          "/circle\n"
          "{ /r exch def \n"
          "  /y exch def\n"
          "  /x exch def\n"
          "  newpath\n"
          "  x r add y moveto\n"
          "  x y r 0 360 arc\n"
          "  linecolor usecolor stroke\n"
          "} def\n"
          "\n"
          "/putpixel {0.1 circle} def\n"
          "\n"
          "/outtextxy\n"
          "{ 3 1 roll moveto\n"
          "  gsave 1 -1 scale\n"
          "  linecolor usecolor\n"
          "  show grestore\n"
          "} def\n"
          "\n"
          "/initialize\n");

  fprintf(filePS, "{ /Helvetica findfont %d scalefont setfont\n", cTextHeight * nResFactor);
  fprintf(filePS, "  %f -%f scale\n", 1./nResFactor, 1./nResFactor);
  fprintf(filePS, "  1.0 setlinewidth\n");
  fprintf(filePS, "  0 -%d translate\n", 792 * nResFactor);
  fprintf(filePS,
          "  15 setlc\n"
          "  0  setfc\n"
          "} def\n"
          "\n"
          "%%%%EndProlog\n"
          "\n"
          "initialize\n");

  cTextWidth  *= nResFactor;
  cTextHeight *= nResFactor;
  return 0;
}


/******************************************************************************
 * Description:  This function closes the postscipt file.
 */
void GRPrintClose ( )
{
  fprintf(filePS, "showpage\n");
  fclose(filePS);
  filePS = NULL;
}

/*****************************************************************************/
/**************************  PUBLIC FUNCTIONS  *******************************/
/*****************************************************************************/


/******************************************************************************
 * Description:  This function initializes the cArg and apszArg variables.
 */
void setcmdargs ( int c, char **apsz )
{
  cArg = c;
  apszArg = apsz;
}


/******************************************************************************
 * Description:  This function initializes the graphics system
 */
void GRInitGraphics ( int TextWidth, int TextHeight )
{
  /* initialize variables */
  /* Font can be "7x14", "9x15", ....  */
  char szFontName[8];           /* font name */
  sprintf(szFontName, "%dx%d", TextWidth, TextHeight);
  pszDisplayName = NULL;
  cBorderWidth = 4;

  /* ------------------  connect to X server ------------------------ */
  if ((display = XOpenDisplay(pszDisplayName)) == NULL)
  {
    fprintf(stdout, "ERROR:  cannot connect to X server %s\n",
            XDisplayName(pszDisplayName));
    exit(2);
  }

#if defined(DEBUG)
  fprintf(stdout, "connected to X server\n");
#endif

  /* fetch the default screen handle as well as the screen dimensions */
  nScreen = DefaultScreen(display);
  cDisplayWidth  = DisplayWidth(display, nScreen);
  cDisplayHeight = DisplayHeight(display, nScreen);

#if defined(DEBUG)
  fprintf(stdout, "screen handle = %d\n", nScreen);
  fprintf(stdout, "display width = %d\n", cDisplayWidth);
  fprintf(stdout, "display height = %d\n", cDisplayHeight);
#endif

  /* load the font */
  if ((xfont = XLoadQueryFont(display, szFontName)) == NULL)
  {
    fprintf(stdout, "ERROR:  cannot open %s font\n", szFontName);
    exit(2);
  }
  cTextWidth = TextWidth;

#if defined(DEBUG)
  fprintf(stdout, "loaded %s font\n", szFontName);
#endif
}


/******************************************************************************
 * Description:  This function sets the current viewport for graphics output.
 */
void GROpenGraphWin ( char *pszWindowTitle, char *pszIconTitle,
                      int left, int top, int right, int bottom,
                      PFREDRAW pf )
{
  int  i;                                /* loop index */
  char ColorName[8];

  /* Origin is at top left corner */
  nXorigin = left;
  nYorigin = top;
  cXresolution = right - left;
  cYresolution = bottom - top;

  /* ---------------------- create window ------------------------- */
  dwForeground = WhitePixel(display, nScreen);
  dwBackground = BlackPixel(display, nScreen);
  window = XCreateSimpleWindow( display, RootWindow(display, nScreen),
                                nXorigin, nYorigin, cXresolution, cYresolution, cBorderWidth,
                                dwForeground, dwBackground);

#if defined(DEBUG)
  fprintf(stdout, "graphics screen coordinates (%d,%d) to (%d,%d)\n",
          left, top, right, bottom);
#endif

  /* load window colors */
  cColorDepth = DisplayPlanes(display, nScreen);
  visual = DefaultVisual(display, nScreen);

  /* Macro returns the default colormap ID */
  colormap = DefaultColormap(display, nScreen);

#if defined(DEBUG)
  fprintf(stdout, "display visual = %d\n", visual);
  fprintf(stdout, "display planes = %d\n", cColorDepth);
#endif

  /* -------------- Set LUT (Look Up Table) for colors -------------- */
  /* monochrome screen : Hardware doesn't support colors */
  if (cColorDepth == 1)
  {
    cMaxColor = 1;
    gc = (GC **) malloc(2 * sizeof(GC *));
    for(i = 0; i < 2; ++i) gc[i] = (GC *) malloc(sizeof(GC));
    GetGC(gc[0], "black");              /* black         */
    GetGC(gc[1], "white");              /* white         */
  }

  /* color screen : Hardware can show different sytles */
  else
    switch (cMaxColor)
    {
    case  99: gc = (GC **) malloc(100 * sizeof(GC *));
      for(i = 0; i < 100; ++i)
      {
        gc[i] = (GC *) malloc(sizeof(GC));
        sprintf(ColorName, "Gray%d", i);
        GetGC(gc[i], ColorName);  /* grey level         */
      }
      break;
    default:  cMaxColor = 15;
      gc = (GC **) malloc((16+Rainbow_Color_NUM) * sizeof(GC *));
      for(i = 0; i < 16+Rainbow_Color_NUM; ++i)
        gc[i] = (GC *) malloc(sizeof(GC));
      GetGC(gc[0],  "black");             /* black         */
      GetGC(gc[1],  "blue");              /* blue          */
      GetGC(gc[2],  "green");             /* green         */
      GetGC(gc[3],  "cyan");              /* cyan          */
      GetGC(gc[4],  "red");               /* red           */
      GetGC(gc[5],  "magenta");           /* magenta       */
      GetGC(gc[6],  "brown");             /* brown         */
      GetGC(gc[7],  "gray");              /* light gray    */
      GetGC(gc[8],  "dark slate gray");   /* dark gray     */
      GetGC(gc[9],  "light blue");        /* light blue    */
      GetGC(gc[10], "lime green");        /* light green   */
      GetGC(gc[11], "slate blue");        /* light cyan    */
      GetGC(gc[12], "medium violet red"); /* light red     */
      GetGC(gc[13], "orange red");        /* light magenta */
      GetGC(gc[14], "yellow");            /* yellow        */
      GetGC(gc[15], "white");             /* white         */
      for(i = 0; i <Rainbow_Color_NUM; ++i)
        GetGCRainbow (gc[15+i+1], rc[i][0],rc[i][1],rc[i][2]);

      break;
    }

  nFillColor    = cMaxColor;
  nCurrentColor = cMaxColor;

#if defined(DEBUG)
  fprintf(stdout, "loaded colors, maximum number of colors = %d\n", cMaxColor);
#endif

  /* create an icon */
  pixmapIcon = XCreateBitmapFromData(display, window, szIconBitmap,
                                     ICON_BITMAP_WIDTH, ICON_BITMAP_HEIGHT);

#if defined(DEBUG)
  fprintf(stdout, "created pixel map\n");
#endif

  /* ----------------- initialize window properties --------------------- */
  sizehints.flags  = USPosition | USSize | PPosition | PSize | PMinSize;
  sizehints.x      = nXorigin;
  sizehints.y      = nYorigin;
  sizehints.width  = cXresolution;
  sizehints.height = cYresolution;
  sizehints.min_width  = cXresolution;
  sizehints.min_height = cYresolution;
  XSetStandardProperties(display, window, pszWindowTitle, pszIconTitle,
                         pixmapIcon, apszArg, cArg, &sizehints);
  XSelectInput(display, window, ExposureMask | KeyPressMask | EnterWindowMask |
               ButtonPressMask | Button1MotionMask | StructureNotifyMask);
  XMapWindow(display, window);
  cTextHeight = xfont->ascent + xfont->descent;

#if defined(DEBUG)
  fprintf(stdout, "window name = %s\n", pszWindowTitle);
  fprintf(stdout, "icon name = %s\n", pszIconTitle);
  fprintf(stdout, "number of command line arguments = %d\n", cArg);
#endif
}


/******************************************************************************
 * Description:  This function clears the current viewport
 */
void GRClearGraphWin ( void )
{
  XClearWindow(display, window);
}


/******************************************************************************
 * Description:  This function shuts down the graphics system
 */
void GRFreeGraphics ( void )
{
  int i;                                /* loop index */

  /* unload the font */
  //XUnloadFont(display, xfont->fid);
  XFreeFont(display, xfont);
  /* unload the colors */
  for(i = 0; i <= cMaxColor; ++i)
  {
    XFreeGC(display, *gc[i]);
    free(gc[i]);
  }
  if(cMaxColor==15)
    for(i = 0; i < Rainbow_Color_NUM; ++i)
    {
      XFreeGC(display, *gc[cMaxColor+1+i]);
      free(gc[cMaxColor+1+i]);
    }
  free(gc);
  /* close the display */
  XCloseDisplay(display);
}


/******************************************************************************
 * Description:  This function sends all queued requests to the server.
 */
void flushdisplay ( void )
{
  XFlush(display);

#if defined(DEBUG)
  fprintf(stdout, "sent all queued requests to the server\n");
#endif
}


/******************************************************************************
 * Description:  This function returns an event.
 */
int getevent ( short int *keycode, int *width, int *height, int *x, int *y)
{
  char buffer[10];
  KeySym keysym;
  XComposeStatus compose;
  int count;

  /* get the event from the X server */
  XWindowEvent(display, window, ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask |
               EnterWindowMask | LeaveWindowMask | Button1MotionMask, &xevent);

  switch (xevent.type)
  {
  case Expose:
    while (XCheckTypedEvent(display, Expose, &xevent));
    return (XCheckWindowEvent(display, window, ButtonPress, &xevent)) ?
           BUTTONPRESS : EXPOSE;

  case EnterNotify:
  case SelectionNotify:
    XRaiseWindow(display, window);
    XSetInputFocus(display, window, RevertToNone, CurrentTime);
    break;

  case LeaveNotify:
    break;

  case ConfigureNotify:
    *width = xevent.xconfigure.width;
    *height = xevent.xconfigure.height;
    return(CONFIGURENOTIFY);

  case KeyPress:
    count = XLookupString(&xevent.xkey, buffer, sizeof(buffer),
                          &keysym, &compose);
    buffer[count] = 0;
    *keycode = (char) 0;
    do
    {
      if       (keysym == XK_Left)        *keycode = LEFTARROW;
      else if  (keysym == XK_Right)       *keycode = RIGHTARROW;
      else if  (keysym == XK_Up)          *keycode = UPARROW;
      else if  (keysym == XK_Down)        *keycode = DOWNARROW;
      else if  (keysym == XK_Escape)      *keycode = ESCAPE;
      else if  (keysym == XK_BackSpace)   *keycode = BACKSPACE;
      else if ((keysym >= XK_space) &&
               (keysym <= XK_asciitilde)) *keycode = buffer[0];
      else                                *keycode = xevent.xkey.keycode;
    }
    while (*keycode == 0);
    return KEYPRESS;

  case ButtonPress:
    *x=xevent.xbutton.x;
    *y=xevent.xbutton.y;
    return BUTTONPRESS;
  case MotionNotify:
    *x=xevent.xmotion.x;
    *y=xevent.xmotion.y;
    return BUTTONMOTION;
  case ButtonRelease:
    *x=xevent.xbutton.x;
    *y=xevent.xbutton.y;
    return BUTTONRELEASE;
  }
  return OTHEREVENT;
}


/******************************************************************************
 * Description:  This function returns the maximum x screen coordinate
 */
int GRGetMaxX ( void )
{
  return (cXresolution <= 0) ? cDisplayWidth : cXresolution;
}


/******************************************************************************
 * Description:  This function returns the maximum y screen coordinate
 */
int GRGetMaxY ( void )
{
  return (cYresolution <= 0) ? cDisplayHeight : cYresolution;
}


/******************************************************************************
 * Description:  This function sets the current drawing color.
 */
void GRSetColor ( int color )
{
  nCurrentColor = color;
  if (filePS) fprintf(filePS, "%d setlc\n", color);
}


/******************************************************************************
 * Description:  This function sets the fill pattern and color.
 */
void GRSetFillColor ( int color, int pattern )
{
  /* X window fill options are FillSolid, FillTiled, */
  /* FillStippled and FillOpaqueStipplied            */

  nFillColor = color;

  if (filePS) fprintf(filePS, "%d setfc\n", color);
  else
  {
    switch (pattern)
    {
    case PS_SOLID:
      XSetFillStyle(display, *gc[color], FillSolid);
      break;

    default:
      XSetFillStyle(display, *gc[color], FillTiled);
      break;
    }
  }
}


/******************************************************************************
 * Description:  This function plots a pixel at a specified point.
 */
void GRPutPixel ( int x, int y, int color )
{
  if (filePS) fprintf(filePS, "%d %d putpixel\n", x, y);
  else XDrawPoint(display, window, *gc[color], x, y);
}


/******************************************************************************
 * Description:  This function draws a line between two specified points.
 */
void GRLine ( int x1, int y1, int x2, int y2 )
{
  if (filePS) fprintf(filePS, "%d %d %d %d line\n", x1, y1, x2, y2);
  else XDrawLine(display, window, *gc[nCurrentColor], x1, y1, x2, y2);
}


/******************************************************************************
 * Description:  This function draws a circle of the given radius with its
 *               center at (x,y).
 */
void GRCircle ( int x, int y, int r )
{
  if (filePS) fprintf(filePS, "%d %d %d circle\n", x, y, r);
  else XDrawArc(display, window, *gc[nCurrentColor], x - r, y - r,
                  2 * r, 2 * r, 0, 23040);
}


/******************************************************************************
 * Description:  This function draws the outline of a polygon.
 */
void GRDrawPoly ( int numpoints, int *polypoints )
{
  int i;                                /* loop index */

  if (filePS)
  {
    for (i = 0; i < numpoints; ++i)
      fprintf(filePS, "%d %d ", polypoints[2*i], polypoints[2*i+1]);
    fprintf(filePS, "%d drawpoly\n", numpoints);
  }
  else
  {
    for (i = 0; i < numpoints - 1; ++i)
      GRLine(polypoints[2*i]  , polypoints[2*i+1],
             polypoints[2*i+2], polypoints[2*i+3]);
  }
}


/******************************************************************************
 * Description:  This function draws and fills a polygon.
 */
void GRFillPoly ( int numpoints, int *polypoints )
{
  XPoint *pxpoint;                      /* pointer to x point structure */
  int i,j;                              /* loop indices */


  if (filePS)
  {
    for (i = 0; i< numpoints; ++i)
      fprintf(filePS, "%d %d ", polypoints[2*i], polypoints[2*i+1]);
    fprintf(filePS, "%d fillpoly\n", numpoints);
  }
  else
  {
    for (i = 0, j = 0, pxpoint = axpoint; i < numpoints; ++i, ++pxpoint)
    {
      pxpoint->x = polypoints[j++];
      pxpoint->y = polypoints[j++];
    }

    /* connect end pooints, if not */
    if ((axpoint[0].x != axpoint[numpoints-1].x) ||
        (axpoint[0].y != axpoint[numpoints-1].y))
    {
      axpoint[numpoints].x = polypoints[0];
      axpoint[numpoints].y = polypoints[1];
      ++numpoints;
    }

    /* x fill mode may be complex, nonconvex or convex */
    XFillPolygon(display, window, *gc[nFillColor], axpoint, numpoints - 1,
                 Nonconvex, CoordModeOrigin);
    XDrawLines(display, window, *gc[nCurrentColor], axpoint, numpoints,
               CoordModeOrigin);
  }
}


/******************************************************************************
 * Description:  This function displays a string at a specified location.
 */
void GRText ( int x, int y, char *textstring )
{
  if (filePS) fprintf(filePS, "%d %d (%s) outtextxy\n", x, y + cTextHeight, textstring);
  else XDrawString(display, window, *gc[nCurrentColor], x, y + cTextHeight,
                     textstring, strlen(textstring));

#if defined(DEBUG)
  fprintf(stdout, "displayed \"%s\" at (%d,%d)\n", textstring, x, y);
#endif
}


/******************************************************************************
 * Description:  This function setup the color LUT (Look Up Table) style.
 *               And return the maxmum color value.
 */
int GRSetLUT ( int monitorstyle )
{
  if (monitorstyle) cMaxColor = 99;             /* graylevel */
  else              cMaxColor = 15;             /* color     */

  return cMaxColor;
}



/*****************************************************************************/
/**************************  PRIVATE FUNCTIONS  ******************************/
/*****************************************************************************/


/******************************************************************************
 * Description:  This function allocates a color.
 */
void GetGC ( GC *gc, char *pszColorName )
{
  DWORD dwValueMask;
  DWORD dwForeground;
  XGCValues values;
  XColor xcolorDef;
  XColor xcolorRGB;

  XAllocNamedColor(display, colormap, pszColorName, &xcolorDef, &xcolorRGB);
  dwForeground = xcolorDef.pixel;
  dwValueMask  = GCForeground | GCBackground;
  values.foreground = dwForeground;
  values.background = dwBackground;
  *gc = XCreateGC(display, window, dwValueMask, &values);
  XSetFont(display, *gc, xfont->fid);
}

/******************************************************************************
 * Description:  This function allocates a rainbow color.
 */
void GetGCRainbow (GC *gc, unsigned short R, unsigned short G,unsigned short B)
{
  DWORD dwValueMask;
  DWORD dwForeground;
  XGCValues values;
  XColor xcolorDef;
  XColor xcolorRGB;
  xcolorDef.red = R;
  xcolorDef.green = G;
  xcolorDef.blue = B;
  xcolorDef.flags = 0;
  XAllocColor(display, colormap, &xcolorDef);
  dwForeground = xcolorDef.pixel;
  dwValueMask  = GCForeground | GCBackground;
  values.foreground = dwForeground;
  values.background = dwBackground;
  *gc = XCreateGC(display, window, dwValueMask, &values);
  XSetFont(display, *gc, xfont->fid);
}



/*****************************************************************************/
/**************************  DUMP SCREEN FUNCTIONS  **************************/
/*****************************************************************************/
#ifdef   HAVE_TIFF
#include <tiffio.h>
int   GRSaveScreen (const char * tiff_file, int width, int height)
{

  TIFF *fout= TIFFOpen(tiff_file, "w");
  if(!fout) return 1;

  int sampleperpixel = 3;    //RGB channel
  XImage *ximage = XGetImage(display, window, 0, 0, width, height, AllPlanes, XYPixmap);
  char *image=new char [width*height*sampleperpixel];

  //it seems the XQueryColors has the limit of query size of 65535 colors.
  XColor  *color=new XColor[65535];
  int lines=65535/width;
  int currentlines;
  int count = 0;

  for (int i = 0; i < height; i+=lines)
  {
    if(i+lines>=height) currentlines=height-i-1;
    else currentlines=lines;
    for (int j = 0; j < currentlines; j++)
      for (int k = 0; k < width; k++)
      {
        unsigned long c;
        c = XGetPixel(ximage, k, i+j);        // X pixel value 
        color[k+j*width].pixel = c;
      }
    XQueryColors(display, colormap, color, currentlines*width);
    for (int j = 0; j < currentlines; j++)
      for (int k = 0; k < width; k++)
      {
        char     r, g, b;
        r = color[k+j*width].red >> 8;
        g = color[k+j*width].green >> 8;
        b = color[k+j*width].blue >> 8;

        image[count++] = r;
        image[count++] = g;
        image[count++] = b;
      }
  }
  delete [] color;

  TIFFSetField(fout, TIFFTAG_IMAGEWIDTH, width);  // set the width of the image
  TIFFSetField(fout, TIFFTAG_IMAGELENGTH, height);    // set the height of the image
  TIFFSetField(fout, TIFFTAG_SAMPLESPERPIXEL, sampleperpixel);   // set number of channels per pixel
  TIFFSetField(fout, TIFFTAG_BITSPERSAMPLE, 8);    // set the size of the channels
  TIFFSetField(fout, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
  TIFFSetField(fout, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(fout, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
  TIFFSetField(fout, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);

  tsize_t linebytes = sampleperpixel*width;     // length in memory of one row of pixel in the image.
  // We set the strip size of the file to be size of one row of pixels
  TIFFSetField(fout, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(fout, width*sampleperpixel));
  //Now writing image to the file one strip at a time
  for (int row = 0; row < height; row++)
  {
    // check the index here, and figure out why not using h*linebytes
    if(TIFFWriteScanline(fout, &image[row*linebytes], row, 0) < 0)
      break;
  }
  TIFFClose(fout);
  delete [] image;
  XDestroyImage(ximage);

  return 0;
}
#else
int   GRSaveScreen (const char * tiff_file, int width, int height)
{
  return 0;
}
#endif

#endif
