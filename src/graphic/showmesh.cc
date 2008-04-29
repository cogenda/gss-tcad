// This file is borrowed from showmesh 1.0
// shipped with easymesh project
// Author:   Bojan NICENO
// Email:    niceno@univ.trieste.it

#include "config.h"
#ifdef HAVE_X11

#include "bsolver.h"
#include "mathfunc.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifndef max
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif
#define SMALL 1e-30
#define GREAT 1e+30

#define MAIN_WDTH     750
#define MAIN_HGHT     500
#define BUTTON_HEIGHT  18
#define BUTTON_WIDTH   80
#define ON              0
#define OFF            -1

#define DELAUNAY  0
#define VORONOI   1
#define MATERIALS 2
#define BOUNDARY  3
#define ZOOM      4
#define MOVE      5
#define FIT       6
#define QUIT      7
#define NBUTTONS  8

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

static Display *display;
static Window  main_win;
static Window  draw_win;
static XTextProperty winname, iconame;
static Pixmap ico_pixm, back_pixm;
static int draw_wdth, draw_hght;
static int scr_num;
static char *prog_name;
static GC          gc_BoW, gc_WoB, gc_XOR, gc_THICK, gc_DASHED, gc_numb;
static GC          gc_color[16+20];
static XFontStruct *text_font, *numb_font;
static Colormap colormap;
static int main_wdth=MAIN_WDTH;
static int main_hght=MAIN_HGHT;
static double  scl, X0, Y0;
static double  xmin,xmax,ymin,ymax;
#define Xmap(x) (int( ((x)-(xmin+xmax)/2)*scl+X0))
#define Ymap(y) (int(-((y)-(ymin+ymax)/2)*scl+Y0))


class Doping2D
{
public:
  int n;
  short int X,Y;
  double doping;
  Doping2D(int _n, short int _XX, short int _YY, double _doping):n(_n),X(_XX),Y(_YY),doping(_doping)
  {}
}
;

inline void swap(Doping2D &A, Doping2D &B)
{
  Doping2D tmp=A;
  A=B;
  B=tmp;
}

struct butt_data
{
  int    x0, y0, hght, wdth, border;
  const char   *caption;
  int    pressed;
}

static butt_data[NBUTTONS]=
  {
    { MAIN_WDTH-BUTTON_WIDTH-20,  10, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Delaunay",  ON},
    { MAIN_WDTH-BUTTON_WIDTH-20,  40, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Voronoi",   OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20,  70, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Materials", ON},
    { MAIN_WDTH-BUTTON_WIDTH-20, 100, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Boundary",  OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 130, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Zoom",      OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 160, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Move",      OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 190, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Fit",       OFF},
    { MAIN_WDTH-BUTTON_WIDTH-20, 250, BUTTON_HEIGHT, BUTTON_WIDTH, 2, "Exit",      OFF},
  };

static Window button[NBUTTONS];

/******************************************************************************
 * Description:  This function allocates a color.
 */
void GetGCbyColorName (Window win, GC *gc, char *pszColorName )
{
  unsigned long dwValueMask;
  unsigned long dwForeground;
  XGCValues values;
  XColor xcolorDef;
  XColor xcolorRGB;

  XAllocNamedColor(display, colormap, pszColorName, &xcolorDef, &xcolorRGB);
  dwForeground = xcolorDef.pixel;
  dwValueMask  = GCForeground | GCBackground;
  values.foreground = dwForeground;
  values.background = BlackPixel(display, scr_num);
  *gc = XCreateGC(display, win, dwValueMask, &values);
  XSetFont(display, *gc, text_font->fid);
}

/******************************************************************************
 * Description:  This function allocates a rainbow color.
 */
void GetGCbyRainbow (Window win, GC *gc, unsigned short R, unsigned short G,unsigned short B)
{
  unsigned long dwValueMask;
  unsigned long dwForeground;
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
  values.background = BlackPixel(display, scr_num);
  *gc = XCreateGC(display, win, dwValueMask, &values);
  XSetFont(display, *gc, text_font->fid);
}

/*========================================================================*/
void getGC ( Window win )
{
  unsigned long valuemask=0;
  XGCValues     values;

  int         dash_offset=0;
  static char dash_list[2]={1,1};
  int         list_len=2;

  /* Normal, thin lines */
  gc_BoW   = XCreateGC ( display, win, valuemask, &values );
  gc_WoB   = XCreateGC ( display, win, valuemask, &values );

  XSetFont ( display, gc_BoW,  text_font->fid );
  XSetFont ( display, gc_WoB,  text_font->fid );

  XSetForeground ( display, gc_BoW, WhitePixel ( display, scr_num ) );
  XSetForeground ( display, gc_WoB, BlackPixel ( display, scr_num ) );

  XSetLineAttributes ( display, gc_BoW, 0, LineSolid, CapRound, JoinRound );
  XSetLineAttributes ( display, gc_WoB, 0, LineSolid, CapRound, JoinRound );

  /* Thick lines */
  gc_THICK = XCreateGC ( display, win, valuemask, &values );
  XSetForeground ( display, gc_THICK, WhitePixel ( display, scr_num ) );
  XSetLineAttributes ( display, gc_THICK, 3, LineSolid, CapRound, JoinRound );

  /* Dashed lines */
  gc_DASHED = XCreateGC ( display, win, valuemask, &values );
  XSetForeground ( display, gc_DASHED, WhitePixel ( display, scr_num ) );
  XSetLineAttributes ( display, gc_DASHED, 0, LineOnOffDash, CapRound, JoinRound );
  XSetDashes ( display, gc_DASHED, dash_offset, dash_list, list_len );

  /* numbers */
  gc_numb = XCreateGC ( display, win, valuemask, &values );
  XSetFont ( display, gc_numb, numb_font->fid );
  XSetForeground ( display, gc_numb, WhitePixel ( display, scr_num ) );

  /* Invisible lines */
  gc_XOR = XCreateGC ( display, win, 0, NULL );
  XSetFunction ( display, gc_XOR, GXxor );
  XSetForeground ( display, gc_XOR, WhitePixel ( display, scr_num ) );

  /* Color GC */
  GetGCbyColorName(win,&gc_color[0],  "black");             /* black         */
  GetGCbyColorName(win,&gc_color[1],  "blue");              /* blue          */
  GetGCbyColorName(win,&gc_color[2],  "green");             /* green         */
  GetGCbyColorName(win,&gc_color[3],  "cyan");              /* cyan          */
  GetGCbyColorName(win,&gc_color[4],  "red");               /* red           */
  GetGCbyColorName(win,&gc_color[5],  "magenta");           /* magenta       */
  GetGCbyColorName(win,&gc_color[6],  "brown");             /* brown         */
  GetGCbyColorName(win,&gc_color[7],  "gray");              /* light gray    */
  GetGCbyColorName(win,&gc_color[8],  "dark slate gray");   /* dark gray     */
  GetGCbyColorName(win,&gc_color[9],  "light blue");        /* light blue    */
  GetGCbyColorName(win,&gc_color[10], "lime green");        /* light green   */
  GetGCbyColorName(win,&gc_color[11], "slate blue");        /* light cyan    */
  GetGCbyColorName(win,&gc_color[12], "medium violet red"); /* light red     */
  GetGCbyColorName(win,&gc_color[13], "orange red");        /* light magenta */
  GetGCbyColorName(win,&gc_color[14], "yellow");            /* yellow        */
  GetGCbyColorName(win,&gc_color[15], "white");             /* white         */
  for(int i = 0; i <Rainbow_Color_NUM; i++)
    GetGCbyRainbow (win,&gc_color[16+i], rc[i][0],rc[i][1],rc[i][2]);

  for(int i=0; i<16+Rainbow_Color_NUM; i++)
    XSetLineAttributes ( display, gc_color[i], 3, LineSolid, CapRound, JoinRound );
}



/*========================================================================*/
void load_fonts()
{
  if ( ( text_font = XLoadQueryFont ( display, "-*-helvetica-bold-r-normal--12-*" ) ) == NULL )
  {
    if ( ( text_font = XLoadQueryFont ( display, "7x14" ) ) == NULL )
    {
      ( void ) fprintf ( stderr, "%s: Cannot open font\n", prog_name );
      exit ( -1 );
    }
  }

  if ( ( numb_font = XLoadQueryFont ( display, "-*-helvetica-bold-r-normal--10-*" ) ) == NULL )
  {
    if ( ( numb_font = XLoadQueryFont ( display, "7x14" ) ) == NULL )
    {
      ( void ) fprintf ( stderr, "%s: Cannot open font\n", prog_name );
      exit ( -1 );
    }
  }
}

/*========================================================================*/
void draw ( Window win, GC gc, int win_x_dim, int win_y_dim )
{
  int x_0, y_0, x_dim, y_dim;

  x_0 = win_x_dim/5;
  y_0 = win_y_dim/5;

  x_dim = 3*win_x_dim/5;
  y_dim = 3*win_y_dim/5;

  XDrawRectangle ( display, win, gc, x_0, y_0, x_dim, y_dim );
}

/*========================================================================*/
void place_text ( Window win, GC gc, XFontStruct *text_font,
                  int win_x_dim, int win_y_dim, char *string )
{
  int width, height; /* string height and width */

  width  = XTextWidth ( text_font, string, strlen ( string ) );
  height = text_font->ascent + text_font->descent;

  XDrawString ( display, win, gc,
                ( win_x_dim-width ) /2, ( win_y_dim+height ) /2,
                string, strlen ( string ) );
}


/*========================================================================*/
void create_buttons ( Window parent )
{
  int b;

  unsigned long        valuemask = CWWinGravity;
  unsigned long        border, background;
  XSetWindowAttributes attr;

  attr.win_gravity = NorthEastGravity;

  for ( b=0; b<NBUTTONS; b++ )
  {
    if ( butt_data[b].pressed==OFF )
    {
      border    =BlackPixel ( display, scr_num );
      background=WhitePixel ( display, scr_num );
    }
    else
    {
      border    =WhitePixel ( display, scr_num );
      background=BlackPixel ( display, scr_num );
    }

    button[b] = XCreateSimpleWindow ( display, parent,
                                      butt_data[b].x0,   butt_data[b].y0,
                                      butt_data[b].wdth, butt_data[b].hght,
                                      butt_data[b].border,
                                      border, background );

    XSelectInput ( display, button[b], ExposureMask | ButtonPressMask );
    XChangeWindowAttributes ( display, button[b], valuemask, &attr );
    XMapWindow ( display, button[b] );
  }
}

/*========================================================================*/
void write_on_button ( int b )
{
  int width, height; /* string height and width */

  unsigned long        valuemask = CWBackPixel | CWBorderPixel;
  XSetWindowAttributes attr;
  GC gc;

  height = text_font->ascent + text_font->descent;
  width  = XTextWidth ( text_font, butt_data[b].caption, strlen ( butt_data[b].caption ) );

  if ( butt_data[b].pressed==ON )
  {
    attr.background_pixel = WhitePixel ( display, scr_num );
    attr.border_pixel     = BlackPixel ( display, scr_num );
  }

  if ( butt_data[b].pressed==OFF )
  {
    attr.background_pixel = BlackPixel ( display, scr_num );
    attr.border_pixel     = WhitePixel ( display, scr_num );
  }

  XChangeWindowAttributes ( display, button[b], valuemask, &attr );
  XClearWindow ( display, button[b] );
  /* XFlush(display); */

  if ( butt_data[b].pressed==ON )  gc=gc_WoB;
  if ( butt_data[b].pressed==OFF ) gc=gc_BoW;

  XDrawString ( display, button[b], gc,
                ( butt_data[b].wdth-width ) /2, ( butt_data[b].hght+height ) /2,
                butt_data[b].caption, strlen ( butt_data[b].caption ) );
}

/*=========================================================================*/
void BSolver ::show_mesh_redwraw()
{
  XPoint axpoint[10];

  /***********************
  *  Draw Delaunay Mesh  *
  ***********************/
  //draw triangle from zone level
  for ( int z = 0; z < zone_num; z++ )
    for ( int itri = 0; itri < zone[z].datri.size(); itri++ )
    {
      Tri ptri = zone[z].datri[itri];
      int Bcolor=GR_WHITE;
      if(zonedata[z]->material_type == Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast<SMCZone *> (zonedata[z]);
        
	/*
	double dopingA = (pzonedata->aux[ptri.node[0]].Nd - pzonedata->aux[ptri.node[0]].Na)*1e18;
        double dopingB = (pzonedata->aux[ptri.node[1]].Nd - pzonedata->aux[ptri.node[1]].Na)*1e18;
        double dopingC = (pzonedata->aux[ptri.node[2]].Nd - pzonedata->aux[ptri.node[2]].Na)*1e18;
        Doping2D A(ptri.node[0], axpoint[0].x,axpoint[0].y,dsign(dopingA)*log(1+fabs(dopingA)));
	Doping2D B(ptri.node[1], axpoint[1].x,axpoint[1].y,dsign(dopingB)*log(1+fabs(dopingB)));
	Doping2D C(ptri.node[2], axpoint[2].x,axpoint[2].y,dsign(dopingC)*log(1+fabs(dopingC)));
	*/
	double dopingA = (pzonedata->aux[ptri.node[0]].Nd - pzonedata->aux[ptri.node[0]].Na);
        double dopingB = (pzonedata->aux[ptri.node[1]].Nd - pzonedata->aux[ptri.node[1]].Na);
        double dopingC = (pzonedata->aux[ptri.node[2]].Nd - pzonedata->aux[ptri.node[2]].Na);
        Doping2D A(ptri.node[0], Xmap(zone[z].danode[ptri.node[0]].x),Ymap(zone[z].danode[ptri.node[0]].y),(dopingA));
	Doping2D B(ptri.node[1], Xmap(zone[z].danode[ptri.node[1]].x),Ymap(zone[z].danode[ptri.node[1]].y),(dopingB));
	Doping2D C(ptri.node[2], Xmap(zone[z].danode[ptri.node[2]].x),Ymap(zone[z].danode[ptri.node[2]].y),(dopingC));
	
        // Sort the three points A B C in increasing doping components;
        if (A.doping > B.doping) swap(A, B);
        if (B.doping > C.doping)
        {
          if(C.doping < A.doping) swap(A, C);
          swap(B, C);
        }

        if(A.doping>=0 || C.doping<=0)             
        {
         if(A.doping>=0)   Bcolor= GR_LIGHTGREEN;  //the triangle is N-Type
         if(C.doping<=0)   Bcolor= Rainbow_15;     //the triangle is P-Type
          
          axpoint[0].x=Xmap(zone[z].danode[ptri.node[0]].x);
          axpoint[0].y=Ymap(zone[z].danode[ptri.node[0]].y);
          axpoint[1].x=Xmap(zone[z].danode[ptri.node[1]].x);
          axpoint[1].y=Ymap(zone[z].danode[ptri.node[1]].y);
          axpoint[2].x=Xmap(zone[z].danode[ptri.node[2]].x);
          axpoint[2].y=Ymap(zone[z].danode[ptri.node[2]].y);
          axpoint[3].x=axpoint[0].x;
          axpoint[3].y=axpoint[0].y;
          if ( butt_data[MATERIALS].pressed==ON )
            XFillPolygon(display, draw_win, gc_color[Bcolor], axpoint, 4 - 1, Nonconvex, CoordModeOrigin);
        }
        //we must draw the metallurgic juction
        if(A.doping<0 && B.doping>=0 && C.doping>=0)
        {
          // B and C above the A
          // A<0<=B<=C
          double factor = (0 - A.doping) / ( B.doping - A.doping);
          int temp1x = int(A.X + (B.X - A.X) * factor);
          int temp1y = int(A.Y + (B.Y - A.Y) * factor);

          factor = (0 - A.doping) / (C.doping - A.doping);
          int temp2x = int(A.X + (C.X - A.X) * factor);
          int temp2y = int(A.Y + (C.Y - A.Y) * factor);
          
          Bcolor= Rainbow_15;
          axpoint[0].x=A.X;
          axpoint[0].y=A.Y;
          axpoint[1].x=temp1x;
          axpoint[1].y=temp1y;
          axpoint[2].x=temp2x;
          axpoint[2].y=temp2y;
          axpoint[3].x=axpoint[0].x;
          axpoint[3].y=axpoint[0].y;
          if ( butt_data[MATERIALS].pressed==ON )
            XFillPolygon(display, draw_win, gc_color[Bcolor], axpoint, 4 - 1, Nonconvex, CoordModeOrigin);
          //if ( butt_data[DELAUNAY].pressed==ON )
          //  XDrawLines(display, draw_win, gc_BoW, axpoint, 4,  CoordModeOrigin);
          
          Bcolor= GR_LIGHTGREEN;
          axpoint[0].x=temp1x;
          axpoint[0].y=temp1y;
          axpoint[1].x=temp2x;
          axpoint[1].y=temp2y;
          axpoint[2].x=C.X;
          axpoint[2].y=C.Y;
          axpoint[3].x=B.X;
          axpoint[3].y=B.Y;
          axpoint[4].x=axpoint[0].x;
          axpoint[4].y=axpoint[0].y;  
          if ( butt_data[MATERIALS].pressed==ON )
            XFillPolygon(display, draw_win, gc_color[Bcolor], axpoint, 5 - 1, Nonconvex, CoordModeOrigin);
          //if ( butt_data[DELAUNAY].pressed==ON )
          //  XDrawLines(display, draw_win, gc_BoW, axpoint, 5,  CoordModeOrigin);  
        }
        if(A.doping<0 && B.doping<0 && C.doping>=0)
        {
          // A and B are below C
          // a<b<0<=c

         double factor = (0 - A.doping) / (C.doping - A.doping);
         int temp1x = int(A.X + (C.X - A.X) * factor);
         int temp1y = int(A.Y + (C.Y - A.Y) * factor);

         factor = (0 - B.doping) / (C.doping - B.doping);
         int temp2x = int(B.X + (C.X - B.X) * factor);
         int temp2y = int(B.Y + (C.Y - B.Y) * factor);
	
	 Bcolor= Rainbow_15;
         axpoint[0].x=temp1x;
         axpoint[0].y=temp1y;
         axpoint[1].x=temp2x;
         axpoint[1].y=temp2y;
         axpoint[2].x=B.X;
         axpoint[2].y=B.Y;
         axpoint[3].x=A.X;
         axpoint[3].y=A.Y;
         axpoint[4].x=axpoint[0].x;
         axpoint[4].y=axpoint[0].y;  
         if ( butt_data[MATERIALS].pressed==ON )
           XFillPolygon(display, draw_win, gc_color[Bcolor], axpoint, 5 - 1, Nonconvex, CoordModeOrigin);
         //if ( butt_data[DELAUNAY].pressed==ON )
         //  XDrawLines(display, draw_win, gc_BoW, axpoint, 5,  CoordModeOrigin);  	
         
          Bcolor= GR_LIGHTGREEN;
          axpoint[0].x=C.X;
          axpoint[0].y=C.Y;
          axpoint[1].x=temp1x;
          axpoint[1].y=temp1y;
          axpoint[2].x=temp2x;
          axpoint[2].y=temp2y;
          axpoint[3].x=axpoint[0].x;
          axpoint[3].y=axpoint[0].y;
          if ( butt_data[MATERIALS].pressed==ON )
            XFillPolygon(display, draw_win, gc_color[Bcolor], axpoint, 4 - 1, Nonconvex, CoordModeOrigin);
          //if ( butt_data[DELAUNAY].pressed==ON )
          //  XDrawLines(display, draw_win, gc_BoW, axpoint, 4,  CoordModeOrigin);
        	
        }	
          axpoint[0].x=Xmap(zone[z].danode[ptri.node[0]].x);
          axpoint[0].y=Ymap(zone[z].danode[ptri.node[0]].y);
          axpoint[1].x=Xmap(zone[z].danode[ptri.node[1]].x);
          axpoint[1].y=Ymap(zone[z].danode[ptri.node[1]].y);
          axpoint[2].x=Xmap(zone[z].danode[ptri.node[2]].x);
          axpoint[2].y=Ymap(zone[z].danode[ptri.node[2]].y);
          axpoint[3].x=axpoint[0].x;
          axpoint[3].y=axpoint[0].y;
        if ( butt_data[DELAUNAY].pressed==ON )
          XDrawLines(display, draw_win, gc_BoW, axpoint, 4,  CoordModeOrigin);
      }
      else //other zones
      {
        axpoint[0].x=Xmap(zone[z].danode[ptri.node[0]].x);
        axpoint[0].y=Ymap(zone[z].danode[ptri.node[0]].y);
        axpoint[1].x=Xmap(zone[z].danode[ptri.node[1]].x);
        axpoint[1].y=Ymap(zone[z].danode[ptri.node[1]].y);
        axpoint[2].x=Xmap(zone[z].danode[ptri.node[2]].x);
        axpoint[2].y=Ymap(zone[z].danode[ptri.node[2]].y);
        axpoint[3].x=axpoint[0].x;
        axpoint[3].y=axpoint[0].y;
	
	switch (zonedata[z]->material_type)
        {
        case Insulator:
          {
            if (zonedata[z]->material==SiO2)
              Bcolor= GR_CYAN;//Rainbow_3;
            else if (zonedata[z]->material==Nitride)
              Bcolor= GR_MAGENTA;
            else   Bcolor= Rainbow_7;
            break;
          }
        case Conductor:
          Bcolor= GR_BLUE;
          break;
        case Vacuum:
          Bcolor= GR_LIGHTCYAN;
          break;
        case PML:
          Bcolor= Rainbow_3;
          break;
        }
        
        if ( butt_data[MATERIALS].pressed==ON )
          XFillPolygon(display, draw_win, gc_color[Bcolor], axpoint, 4 - 1, Nonconvex, CoordModeOrigin);
        if ( butt_data[DELAUNAY].pressed==ON )
          XDrawLines(display, draw_win, gc_BoW, axpoint, 4,  CoordModeOrigin);
      }
    }

  /**********************
  *  Draw Voronoi Mesh  *
  **********************/
  if ( butt_data[VORONOI].pressed==ON )
    for(int z = 0; z < zone_num; z++)
    {
      for(int i = 0; i < zone[z].davedge.size(); i++)
      {
        VEdge Edge = zone[z].davedge[i];
        int sx1=Xmap(Edge.x1);
        int sy1=Ymap(Edge.y1);
        int sx2=Xmap(Edge.x2);
        int sy2=Ymap(Edge.y2);
        XDrawLine ( display, draw_win, gc_DASHED, sx1, sy1, sx2, sy2);
      }
    }

  /*****************************
  *  Draw Boundary Conditions  *
  *****************************/
  if ( butt_data[BOUNDARY].pressed==ON )
    for(int z = 0; z < zone_num; z++)
    {
      for(int i=0;i<zone[z].dasegment.size();i++)
      {
        int color=GR_WHITE;
        int bc_index = zone[z].dasegment[i].bc_index-1;
        if(bc[bc_index].BCType & ContactMask)   color=GR_BLUE;
        if(bc[bc_index].BCType & InterfaceMask) color=GR_RED;

        for(int k=0;k<zone[z].dasegment[i].edge_num;k++)
        {
          int p1=zone[z].dasegment[i].edge_array[k].p1;
          int p2=zone[z].dasegment[i].edge_array[k].p2;
          int sx1=Xmap(zone[z].danode[p1].x);
          int sy1=Ymap(zone[z].danode[p1].y);
          int sx2=Xmap(zone[z].danode[p2].x);
          int sy2=Ymap(zone[z].danode[p2].y);
          //XDrawLine ( display, draw_win, gc_THICK, sx1, sy1, sx2, sy2);
          XDrawLine ( display, draw_win, gc_color[color], sx1, sy1, sx2, sy2);
        }
      }
    }
  XFlush(display);
}
/*-show_mesh_redwraw--------------------------------------------------------------*/




/*========================================================================*/
void show_mesh_init()
{
  XWindowAttributes main_win_attr;

  /*----------------------+
  |  Connect to X server  |
  +----------------------*/
  {
    char *disp_name=NULL; /* ako nije definirano od korisnika mora biti NULL */

    if ( ( display=XOpenDisplay ( disp_name ) ) == NULL )
    {
      fprintf ( stderr, "%s: cannot connect to X server %s\n",
                prog_name, XDisplayName ( disp_name ) );
      exit ( -1 );
    }
  }

  scr_num = DefaultScreen ( display );

  /*---------------------------------------+
  |  Creating the main application window  |
  +---------------------------------------*/
  main_win = XCreateSimpleWindow ( display, RootWindow ( display, scr_num ),
                                   0, 0, MAIN_WDTH, MAIN_HGHT, 4,
                                   BlackPixel ( display, scr_num ),
                                   WhitePixel ( display, scr_num ) );
  XGetWindowAttributes ( display, main_win, &main_win_attr );
  main_wdth=main_win_attr.width;
  main_hght=main_win_attr.height;
  colormap = DefaultColormap(display, scr_num);
  /* load window colors */
  int ColorDepth = DisplayPlanes(display, scr_num);
  Visual *visual = DefaultVisual(display, scr_num);
  /*------------------------------------+
  |  Preparing an icon an a background  |
  +------------------------------------*/
  {
#define icon_width 57
#define icon_height 57
    static  char icon_bits[] =
      {
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0xfc, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x00,
        0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0xe4, 0xff, 0xd7, 0xff,
        0xaf, 0xff, 0x4f, 0x00, 0xd4, 0xff, 0xbb, 0xff, 0x77, 0xff, 0x57, 0x00,
        0xb4, 0xff, 0xbb, 0xff, 0x7b, 0xff, 0x5b, 0x00, 0x74, 0xff, 0x7d, 0xff,
        0xfb, 0xfe, 0x5d, 0x00, 0xf4, 0xfc, 0xfd, 0xfe, 0xfd, 0xfe, 0x5e, 0x00,
        0xf4, 0xfb, 0xfe, 0xfe, 0xfe, 0x7d, 0x5f, 0x00, 0xf4, 0xf7, 0xfe, 0x7d,
        0xff, 0xbd, 0x5f, 0x00, 0xf4, 0x6f, 0xff, 0x7b, 0xff, 0xdb, 0x5f, 0x00,
        0xf4, 0x5f, 0xff, 0xbb, 0xff, 0xeb, 0x5f, 0x00, 0xf4, 0x1f, 0xfc, 0xd7,
        0xff, 0xf0, 0x5f, 0x00, 0xf4, 0xa7, 0x03, 0xd6, 0x00, 0xc7, 0x5f, 0x00,
        0xf4, 0x79, 0xff, 0x01, 0xff, 0xbb, 0x5f, 0x00, 0x74, 0x7e, 0xff, 0xd7,
        0xff, 0x7b, 0x5e, 0x00, 0x94, 0x7f, 0xff, 0xbb, 0xff, 0xfb, 0x59, 0x00,
        0xe4, 0x7f, 0xff, 0xbb, 0xff, 0xfd, 0x57, 0x00, 0xf4, 0xff, 0xfe, 0x7d,
        0xff, 0xfd, 0x4f, 0x00, 0xe4, 0xff, 0xfe, 0xfe, 0xfe, 0xfd, 0x5f, 0x00,
        0x94, 0xff, 0x7e, 0xff, 0xfd, 0xfe, 0x47, 0x00, 0x74, 0xfe, 0x7e, 0xff,
        0xfd, 0xfe, 0x59, 0x00, 0xf4, 0xfd, 0xbd, 0xff, 0xfb, 0x7e, 0x5e, 0x00,
        0xf4, 0xf3, 0xdd, 0xff, 0x77, 0x9f, 0x5f, 0x00, 0xf4, 0xcf, 0xed, 0xff,
        0x6f, 0xe7, 0x5f, 0x00, 0xf4, 0x3f, 0xed, 0xff, 0x6f, 0xf9, 0x5f, 0x00,
        0xf4, 0xff, 0xf0, 0xff, 0x1f, 0xfe, 0x5f, 0x00, 0xf4, 0xff, 0x03, 0x00,
        0x80, 0xff, 0x5f, 0x00, 0xf4, 0xff, 0xf0, 0xff, 0x1f, 0xfe, 0x5f, 0x00,
        0xf4, 0x3f, 0xf5, 0xff, 0x5f, 0xf9, 0x5f, 0x00, 0xf4, 0xcf, 0xed, 0xff,
        0x6f, 0xe7, 0x5f, 0x00, 0xf4, 0xf3, 0xdd, 0xff, 0x77, 0xdf, 0x5f, 0x00,
        0xf4, 0xfd, 0xbd, 0xff, 0x7b, 0x3f, 0x5f, 0x00, 0x74, 0xfe, 0xbe, 0xff,
        0xfb, 0xfe, 0x5c, 0x00, 0x94, 0xff, 0x7e, 0xff, 0xfd, 0xfe, 0x53, 0x00,
        0xe4, 0xff, 0xfe, 0xfe, 0xfe, 0xfe, 0x4f, 0x00, 0xf4, 0xff, 0xfe, 0xfe,
        0xfe, 0xfe, 0x5f, 0x00, 0xe4, 0x7f, 0xff, 0x7d, 0xff, 0xfd, 0x47, 0x00,
        0x94, 0x7f, 0xff, 0xbb, 0xff, 0xfd, 0x59, 0x00, 0x74, 0x7e, 0xff, 0xd7,
        0xff, 0x7d, 0x5e, 0x00, 0xf4, 0x79, 0xff, 0xd7, 0xff, 0x9d, 0x5f, 0x00,
        0xf4, 0xa7, 0x1f, 0x00, 0xf8, 0xe3, 0x5f, 0x00, 0xf4, 0x1f, 0xe0, 0xd7,
        0x07, 0xf8, 0x5f, 0x00, 0xf4, 0x9f, 0xff, 0xbb, 0xff, 0xf3, 0x5f, 0x00,
        0xf4, 0x6f, 0xff, 0xbb, 0xff, 0xed, 0x5f, 0x00, 0xf4, 0x77, 0xff, 0x7d,
        0xff, 0xdd, 0x5f, 0x00, 0xf4, 0xfb, 0xfe, 0xfe, 0xfe, 0xbe, 0x5f, 0x00,
        0xf4, 0xfc, 0x7e, 0xff, 0xfd, 0x7e, 0x5f, 0x00, 0x74, 0xff, 0x7d, 0xff,
        0x7d, 0xff, 0x5c, 0x00, 0xb4, 0xff, 0xbd, 0xff, 0x7b, 0xff, 0x5b, 0x00,
        0xd4, 0xff, 0xdb, 0xff, 0xb7, 0xff, 0x57, 0x00, 0xe4, 0xff, 0xeb, 0xff,
        0xaf, 0xff, 0x4f, 0x00, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00,
        0xfc, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
      };

#define back_width 4
#define back_height 4
    static  char back_bits[] = {0x05, 0x0a, 0x05, 0x0a};

    ico_pixm = XCreateBitmapFromData ( display, main_win,
                                       icon_bits, icon_width, icon_height );

    back_pixm = XCreatePixmapFromBitmapData ( display, main_win,
                back_bits, back_width, back_height,
                BlackPixel ( display, scr_num ), WhitePixel ( display, scr_num ),
                main_win_attr.depth );
  }

  /*----------------------------+
  |  Set the window background  |
  +----------------------------*/
  XSetWindowBackgroundPixmap ( display, main_win, back_pixm );

  /*------------------------+
  |  Set window properties  |
  +------------------------*/
  {
    XSizeHints    size_hints;
    XWMHints      wm_hints;
    XClassHint    class_hints;

    char *win_name="Mesh Topology";
    char *ico_name="ShowMesh";

    size_hints.flags      = PPosition | PSize | PMinSize;
    size_hints.min_width  = 50; /* 300 */
    size_hints.min_height = NBUTTONS* ( BUTTON_HEIGHT+12 ) +10;

    wm_hints.initial_state = NormalState;  /* Normal or Iconified     */
    wm_hints.input         = True;         /* It needs keyboard input */
    wm_hints.icon_pixmap   = ico_pixm;
    wm_hints.flags         = StateHint | IconPixmapHint | InputHint;

    class_hints.res_name  = prog_name;
    class_hints.res_class = "ShowMesh";

    if ( XStringListToTextProperty ( &win_name,1,&winname ) == 0 )
    {
      fprintf ( stderr, "%s: structure allocation for window name failed\n",
                prog_name );
      exit ( -1 );
    }

    if ( XStringListToTextProperty ( &ico_name,1,&iconame ) == 0 )
    {
      fprintf ( stderr, "%s: structure allocation for icon name failed\n",
                prog_name );
      exit ( -1 );
    }

    XSetWMProperties ( display, main_win,
                       &winname, &iconame,
                       0, NULL,
                       &size_hints, &wm_hints, &class_hints );

  }

  /*---------------------+
  |  Select input types  |
  +---------------------*/
  XSelectInput ( display, main_win, ExposureMask | StructureNotifyMask | KeyPressMask);

  load_fonts();
  getGC ( main_win );

  /*---------------------+
  |  Display the window  |
  +---------------------*/
  XMapWindow ( display, main_win );

  /*--------------------+
  |  Create the button  |
  +--------------------*/
  create_buttons ( main_win );

  /*----------------------------+
  |  Create the drawing window  |
  +----------------------------*/
  draw_wdth = main_wdth-BUTTON_WIDTH-40;
  draw_hght = main_hght-20;
  draw_win = XCreateSimpleWindow ( display, main_win,
                                   10, 10, draw_wdth, draw_hght, 3,
                                   WhitePixel ( display, scr_num ),
                                   BlackPixel ( display, scr_num ) );

  XSelectInput ( display, draw_win, ExposureMask | PointerMotionMask | ButtonPressMask );
  XMapWindow ( display, draw_win );

}
/*------------------------------------------------------------------------*/

void show_mesh_close()
{
  XFree(winname.value);
  XFree(iconame.value);
  XFreePixmap(display,ico_pixm);
  XFreePixmap(display,back_pixm);
  XUnloadFont(display, text_font->fid);
  XUnloadFont(display, numb_font->fid);
  XFreeGC(display, gc_WoB);
  XFreeGC(display, gc_BoW);
  XFreeGC(display, gc_XOR);
  XFreeGC(display, gc_THICK);
  XFreeGC(display, gc_DASHED);
  XFreeGC(display, gc_numb);
  for(int i=0;i<16+Rainbow_Color_NUM;i++)
    XFreeGC(display, gc_color[i]);
  XCloseDisplay(display);
}

#ifdef   HAVE_TIFF
#include <tiffio.h>
int show_mesh_save_screen (const char * tiff_file, Window window, int width, int height)
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
        c = XGetPixel(ximage, k, i+j);        /* X pixel value */
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
/* do nothing */
int show_mesh_save_screen (const char * tiff_file, Window window, int width, int height)
{
  return 0;
}
#endif

/*========================================================================*/
void BSolver::show_mesh(PlotDefine &plot)
{
  XEvent      report;

  /*-------------------------------+
  |  Copy the name of the program  |
  +-------------------------------*/
  prog_name="ShowMesh";
  show_mesh_init();
  xmin=xMin;xmax=xMax;ymin=yMin;ymax=yMax;
  scl =min ( ( 0.9* ( double ) draw_hght ) / ( yMax-yMin+SMALL ),
             ( 0.9* ( double ) draw_wdth ) / ( xMax-xMin+SMALL ) );
  X0 = ( int ) (0.5*draw_wdth);
  Y0 = ( int ) (0.5*draw_hght);

  /*============#
  #  MAIN LOOP  #
  #============*/
  {
    Window        root, child; /* for XQuerryPointer */
    unsigned int  mouse_butt;
    int           b, x0=OFF, y0, x_root, y_root;
    int           x_new, y_new;
    int           x_old=int(0.45*draw_wdth), y_old=int(0.45*draw_hght);
    double        x0_fiz, y0_fiz, scl_new;
    char buffer[10];
    int count;
    KeySym keysym;
    XComposeStatus compose;
    short int keycode;

    while ( 1 )
    {
      XNextEvent ( display, &report );
      switch ( report.type )
      {
        /******************
        *  Expose Window  *
        ******************/
      case Expose:
        while (XCheckTypedWindowEvent(display, report.xany.window, Expose, &report));
        if ( report.xany.window == draw_win )
          show_mesh_redwraw();
        for ( b=0; b<NBUTTONS; b++ )
          if ( report.xany.window==button[b] )
            write_on_button ( b );
        break;

        /******************
        *  Resize Window  *
        ******************/
      case ConfigureNotify:
        main_wdth = report.xconfigure.width;
        main_hght = report.xconfigure.height;
        draw_wdth = main_wdth-BUTTON_WIDTH-40;
        draw_hght = main_hght-20;
        XResizeWindow ( display, draw_win, draw_wdth, draw_hght );
        break;

        /*****************
        *  Button Press  *
        *****************/
      case ButtonPress:
        if ( report.xany.window==button[QUIT] )
        {
          if(plot.GenerateTIFF)
            show_mesh_save_screen(plot.TIFFFileName.c_str(),draw_win, draw_wdth, draw_hght);
          show_mesh_close();
          return;
        }

        for ( b=0; b<NBUTTONS; b++ )
        {
          if ( report.xany.window==button[b] )
          {
            if ( butt_data[b].pressed==ON )
            {butt_data[b].pressed=OFF; write_on_button ( b );}
            else
            {butt_data[b].pressed=ON;  write_on_button ( b );}
          }
        }

        if ( report.xany.window==button[MATERIALS] && butt_data[MATERIALS].pressed==ON )
        { write_on_button ( MATERIALS );}

        if ( report.xany.window==button[FIT] )
        {
          butt_data[FIT].pressed=OFF;
          scl =min ( ( 0.9* ( double ) draw_hght ) / ( yMax-yMin+SMALL ),
                     ( 0.9* ( double ) draw_wdth ) / ( xMax-xMin+SMALL ) );
          X0 = ( int ) (0.5*draw_wdth);
          Y0 = ( int ) (0.5*draw_hght);
          XClearWindow ( display, draw_win );
          show_mesh_redwraw();
          write_on_button ( FIT );
        }

        if ( report.xany.window>=button[DELAUNAY] && report.xany.window<=button[BOUNDARY] )
        {
          XClearWindow ( display, draw_win );
          show_mesh_redwraw();
        }

        if ( report.xany.window==draw_win )
        {
          if ( butt_data[MOVE].pressed==ON )
          {
            if ( x0==OFF ) {x0=report.xmotion.x; y0=report.xmotion.y;}
            else
            {
              butt_data[MOVE].pressed=OFF;
              write_on_button ( MOVE );
              X0+= ( x_new-x0 ); Y0+= ( y_new-y0 );   x0=OFF;
              XClearWindow ( display, draw_win ); show_mesh_redwraw();
            }
          }
          if ( butt_data[ZOOM].pressed==ON )
          {
            if ( x0==OFF ) {x0=report.xmotion.x; y0=report.xmotion.y;}
            else
            {
              butt_data[ZOOM].pressed=OFF;
              write_on_button ( ZOOM );
              x0_fiz = ( min ( x0, x_new )-X0 ) /scl;
              y0_fiz = ( min ( y0, y_new )-Y0 ) /scl;
              if ( x0!=x_new && y0!=y_new )
                scl_new= ( min ( ( double ) draw_wdth/abs ( x0-x_new ),
                                 ( double ) draw_hght/abs ( y0-y_new ) ) ) *scl;
              if ( max ( scl_new*xMax, scl_new*yMax ) < 32768 )
              {
                scl=scl_new;
                X0=-x0_fiz*scl;
                Y0=-y0_fiz*scl;
              }
              x0=OFF;
              XClearWindow ( display, draw_win ); show_mesh_redwraw();
            }
          }
        }
        break;

        /*****************
        *  Mouse Motion  *
        *****************/
      case MotionNotify:
        if ( report.xany.window==draw_win )
        {
          if ( butt_data[MOVE].pressed==ON )
          {
            x_old=x_new; y_old=y_new;
            XQueryPointer ( display, report.xmotion.window,
                            &root, &child, &x_root, &y_root, &x_new, &y_new,
                            &mouse_butt );
            if ( x0!=OFF )
            {
              XDrawLine ( display, draw_win, gc_XOR, x0, y0, x_old, y_old );
              XDrawLine ( display, draw_win, gc_XOR, x0, y0, x_new, y_new );
            }
          }

          if ( butt_data[ZOOM].pressed==ON )
          {
            x_old=x_new; y_old=y_new;
            XQueryPointer ( display, report.xmotion.window,
                            &root, &child, &x_root, &y_root, &x_new, &y_new,
                            &mouse_butt );
            if ( x0!=OFF )
            {
              XDrawRectangle ( display, draw_win, gc_XOR,
                               min ( x0, x_old ), min ( y0, y_old ), abs ( x_old-x0 ), abs ( y_old-y0 ) );
              XDrawRectangle ( display, draw_win, gc_XOR,
                               min ( x0, x_new ), min ( y0, y_new ), abs ( x_new-x0 ), abs ( y_new-y0 ) );
            }
          }
        }

        break;
        /*****************
        *  Key Press     *
        *****************/
      case KeyPress:
        count = XLookupString(&report.xkey, buffer, sizeof(buffer),
                              &keysym, &compose);
        buffer[count] = 0;
        keycode = (char) 0;
        do
        {
          switch (keysym)
          {
          case XK_Left:
            X0--;XClearWindow ( display, draw_win ); show_mesh_redwraw(); break;
          case XK_Right:
            X0++;XClearWindow ( display, draw_win ); show_mesh_redwraw(); break;
          case XK_Up:
            Y0--;XClearWindow ( display, draw_win ); show_mesh_redwraw(); break;
          case XK_Down:
            Y0++;XClearWindow ( display, draw_win ); show_mesh_redwraw(); break;
          case XK_R:
          case XK_r:
            scl =min ( ( 0.9* ( double ) draw_hght ) / ( yMax-yMin+SMALL ),
                       ( 0.9* ( double ) draw_wdth ) / ( xMax-xMin+SMALL ) );
            X0 = ( int ) (0.5*draw_wdth);
            Y0 = ( int ) (0.5*draw_hght);
            XClearWindow ( display, draw_win );
            show_mesh_redwraw();
            break;
          case XK_Z:
            x0_fiz = ( 0.5*draw_wdth-X0 ) /scl;
            y0_fiz = ( 0.5*draw_hght-Y0 ) /scl;
            scl*=2;
            X0+=x0_fiz*scl;
            Y0+=y0_fiz*scl;
            XClearWindow ( display, draw_win ); show_mesh_redwraw();
            break;
          case XK_z:
            x0_fiz = ( 0.5*draw_wdth-X0 ) /scl;
            y0_fiz = ( 0.5*draw_hght-Y0 ) /scl;
            scl/=2;
            X0+=x0_fiz*scl;
            Y0+=y0_fiz*scl;
            XClearWindow ( display, draw_win ); show_mesh_redwraw();
            break;
          case XK_Q:
          case XK_q:
          case XK_Escape:
            if(plot.GenerateTIFF)
              show_mesh_save_screen(plot.TIFFFileName.c_str(),draw_win, draw_wdth, draw_hght);
            show_mesh_close();
            return;
          case XK_BackSpace: break;
          default: break;
          }
          keycode = report.xkey.keycode;
        }
        while (keycode == 0);

        break;
      } /* end switch */
    } /* end while */
  } /* end of block */

}

#endif
