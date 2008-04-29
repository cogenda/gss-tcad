#include "config.h"

#ifdef HAVE_WIN32 

#include <windows.h>
#include <string.h>
#include <stdio.h>
#include "wgraph.h"


static char szClassName[] = "GraphicsWindow";
static FILE *filePS = NULL;
static int  cTextHeight;

static COLORREF acolorref[36] =
{ RGB(0xFF,0xFF,0xFF),    /* white */
  RGB(0x00,0x00,0xFF),    /* blue */
  RGB(0x00,0xFF,0x00),    /* green */
  RGB(0x00,0xFF,0xFF),    /* cyan */
  RGB(0xFF,0x00,0x00),    /* red */
  RGB(0xFF,0xFF,0x43),    /* magenta */
  RGB(0x80,0x66,0x4C),    /* brown */
  RGB(0xB2,0xB2,0xB2),    /* light gray */
  RGB(0x4C,0x4C,0x4C),    /* dark gray */
  RGB(0x43,0x43,0xFF),    /* light blue */
  RGB(0x43,0xFF,0x43),    /* light green */
  RGB(0x80,0xFF,0xFF),    /* light cyan */
  RGB(0xFF,0x43,0x43),    /* light red */
  RGB(0xFF,0xFF,0x80),    /* light magenta */
  RGB(0xFF,0xFF,0x00),    /* yellow */
  RGB(0x00,0x00,0x00),    /* black */
  RGB(0,     0, 255),
  RGB(0,   114, 255),
  RGB(0,   161, 255),
  RGB(0,   197, 255),
  RGB(0,   228, 255),
  RGB(0,   255, 255),
  RGB(0,   255, 221),
  RGB(0,   255, 179),
  RGB(0,   255, 127),
  RGB(0,   255,   0),
  RGB(125, 255,   0),
  RGB(180, 255,   0),
  RGB(220, 255,   0),
  RGB(255, 255,   0),
  RGB(255, 233,   0),
  RGB(255, 207,   0),
  RGB(255, 178,   0),
  RGB(255, 148,   0),
  RGB(255, 103,   0),
  RGB(255,   0,   0)

};

static HDC hDC;
static HPEN ahpen[36];
static HBRUSH ahbrush[36];
static POINT apoint[256];
PFREDRAW pfRedrawFunction;

LRESULT CALLBACK _export WndProc (HWND hWnd, UINT iMessage, WPARAM wParam, LPARAM lParam)
{
  PAINTSTRUCT ps;
  HPEN hpen;
  HBRUSH hbrush;

  switch (iMessage)
  { case WM_CREATE:
      break;

    case WM_PAINT:
      hDC = BeginPaint(hWnd, &ps);
      hpen = (HPEN) SelectObject(hDC, ahpen[GR_WHITE]);
      hbrush = (HBRUSH) SelectObject(hDC, ahbrush[GR_BLACK]);
      (*pfRedrawFunction)();
      SelectObject(hDC, hpen);
      SelectObject(hDC, hbrush);
      EndPaint(hWnd, &ps);
      break;

    case WM_DESTROY:
                PostQuitMessage(0);
                break;
    default:
                return DefWindowProc(hWnd, iMessage, wParam, lParam);
  }
  return 0;
}

int GRPrintOpen (const char *pv, int nResFactor)
{
  char pszFileName[128];
  sprintf(pszFileName, "%s.ps",  pv);
  filePS = fopen(pszFileName, "wt");
  if(!filePS) return 1;
  
  fprintf(filePS,
         "%%!PS-Adobe-1.0\n"
         "%%%%Creator: SGFramework\n"
         "%%%%Title:  %s\n"
         "%%%%Pages: 1\n"
         "%%%%BoundingBox: 72 162 540 640 \n"
         "%%%%EndComments\n\n", pszFileName);
  fprintf(filePS,
         "%%! SGFramework prolog Version 1.0\n"
         "%%%%\n"
         "\n"
         "/setlc {/linecolor exch def} def\n"
         "/setfc {/fillcolor exch def} def\n"
         "\n"
         "/colors\n"
         "[ [1.0 1.0 1.0] %% white\n"
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

  cTextHeight *= nResFactor;
  return 0;
}


void GRPrintClose ()
{
  fprintf(filePS, "showpage\n");
  fclose(filePS);
  filePS = NULL;
}


void GRInitGraphics (int TextWidth, int TextHeight)
{
  int i;
  cTextHeight = TextHeight;
  WNDCLASS wndclass;

  wndclass.style         = CS_HREDRAW | CS_VREDRAW;
  wndclass.lpfnWndProc   = WndProc;
  wndclass.cbClsExtra    = 0;
  wndclass.cbWndExtra    = 0;
  wndclass.hInstance     = 0;
  wndclass.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
  wndclass.hCursor       = LoadCursor(NULL, IDC_ARROW);
  wndclass.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
  wndclass.lpszMenuName  = NULL;
  wndclass.lpszClassName = szClassName;
  if (!RegisterClass(&wndclass)) exit(FALSE);

  for(i = 0; i < 36; ++i)
  { ahpen[i] = CreatePen(PS_SOLID, 1, acolorref[i]);
    ahbrush[i] = CreateSolidBrush(acolorref[i]);
  }
}


void GRFreeGraphics ()
{
  int i;

  for(i = 0; i < 36; ++i)
  { DeleteObject(ahpen[i]);
    DeleteObject(ahbrush[i]);
  }
}


void GROpenGraphWin (char *pszWindowTitle, char *pszIconTitle, int x, int y, int w, int h, PFREDRAW pf)
{
  MSG msg;
  HWND hWnd;

  if (x < 0) { x = CW_USEDEFAULT; y = 0; }
  pfRedrawFunction = pf;
  hWnd = CreateWindow(szClassName, pszWindowTitle, WS_TILED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX, x, y, w, h, NULL, NULL, 0, NULL);
  if (!hWnd) exit(FALSE);

  ShowWindow(hWnd, SW_SHOW);
  UpdateWindow(hWnd);

  while (GetMessage(&msg, NULL, 0, 0))
  { TranslateMessage(&msg);
    DispatchMessage(&msg);
  }
}


void GRSetColor (int nColor)
{
  if (filePS) fprintf(filePS, "%d setlc\n", nColor);
  else SelectObject(hDC, ahpen[nColor]);
}


void GRSetFillColor (int nColor, int nPattern)
{
  if (filePS) fprintf(filePS, "%d setfc\n", nColor);
  else SelectObject(hDC, ahbrush[nColor]);
}


void GRPutPixel (int nColor, int x, int y)
{
  GRSetColor(nColor);
  if (filePS) fprintf(filePS, "%d %d putpixel\n");
  else GRCircle(x, y, 1);
}


void GRLine (int x1, int y1, int x2, int y2)
{
  if (filePS) fprintf(filePS, "%d %d %d %d line\n", x1, y1, x2, y2);
  else
  {     MoveToEx(hDC, x1, y1, NULL);
          LineTo(hDC, x2, y2);
  }
}


void GRCircle (int x, int y, int r)
{
  if (filePS) fprintf(filePS, "%d %d %d circle\n", x, y, r);
  else AngleArc(hDC, x, y, r, 0, 360);
}


void GRDrawPoly (int c, int *an)
{
  int i, cT;
  int *pn = an;
  POINT* ppoint;
        if (filePS)
  { cT = c << 1;
        for(i = 0; i < cT; ++i)
        fprintf(filePS, "%d ", *pn++);
    fprintf(filePS, "%d drawpoly\n", c);
  }
  else
  {     ppoint = apoint;
          for(i = 0; i < c; ++i)
          { ppoint->x = *pn++;
            ppoint->y = *pn++;
            ++ppoint;
          }
          Polyline(hDC, apoint, c);
  }
}


void GRFillPoly (int c, int *an)
{
  int i, cT;
  int *pn = an;
  POINT* ppoint;
        if (filePS)
  { cT = c << 1;
        for(i = 0; i < cT; ++i)
        fprintf(filePS, "%d ", *pn++);
    fprintf(filePS, "%d fillpoly\n", c);
  }
  else
  {     ppoint = apoint;
          for(i = 0; i < c; ++i)
          { ppoint->x = *pn++;
            ppoint->y = *pn++;
            ++ppoint;
          }
          Polygon(hDC, apoint, c);
  }
}


void GRText (int x, int y, char *psz)
{
  if (filePS) fprintf(filePS, "%d %d (%s) outtextxy\n", x, y + cTextHeight, psz);
  else TextOut(hDC, x, y, psz, strlen(psz));
}

int GRSetLUT(int fGrayLevel)
{
  return 15;
}


#endif




