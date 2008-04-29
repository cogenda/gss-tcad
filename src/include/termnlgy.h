// ----------------------------------------------------------------------------
// PROJECT:    SIMULATION GENERATION FRAMEWORK (SIMGEN)
//
// AUTHOR:     Kevin M. Kramer
//             J. Joshua Feng
//
// REVISION HISTORY
// DATE      INITIALS   DESCRIPTION
// --------  --------   ------------------------------------------------------
// 09/15/94    KMK      initial implementation
// ??/??/96    JJF      seperate it from the headers
//
// DESCRIPTION: Terminology / Nomenclature in programming
//              Type definitions and macro definitions, if handled by the
//              precompiler, should not increase the size of the executable,
//              nor the objective code.
// ----------------------------------------------------------------------------


#ifndef __termnlgy_h
#define __termnlgy_h

///////////////////////////////////////////////////////////////////////////////
////////////////////////////  TYPE DEFINITIONS  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// boolean constants
#define TRUE           1
#define FALSE          0

// global type definitions
typedef int            BOOL;                    // f
typedef unsigned char  BYTE;                    // b
typedef unsigned int   UINT;                    // u
typedef char *         PSTR;                    // psz
typedef int  *         PINT;                    // pn
typedef long *         PLONG;                   // pl

typedef void *         PVOID;                   // pv
typedef int *          PBOOL;                   // pf

typedef int            I3ARRAY[3];              // i3a x1, x2, x3 coordinates 
typedef void          (*PFV) (BOOL);

#if   defined(SPREC)
typedef float          REAL;                                      // r
typedef float *        PREAL;                                     // pr
typedef float         (*PFEQU) (I3ARRAY);                 
typedef void          (*PFPD)  (float *, UINT *, I3ARRAY);
#elif defined(LPREC)
typedef long double    REAL;                                      // r
typedef long double *  PREAL;                                     // pr
typedef long double   (*PFEQU) (I3ARRAY);
typedef void          (*PFPD)  (long double *, UINT *, I3ARRAY);
#else
typedef double         REAL;                                      // r
typedef double *       PREAL;                                     // pr
typedef double        (*PFEQU) (I3ARRAY);
typedef void          (*PFPD)  (double *, UINT *, I3ARRAY);
#endif

typedef double        (*PRF) ();                                       // prf
typedef double        (*PRF1)(double);                                 // prf1
typedef double        (*PRF2)(double, double);                         // prf2
typedef double        (*PRF3)(double, double, double);                 // prf3
typedef double        (*PRF4)(double, double, double, double);         // prf4
typedef double        (*PRF5)(double, double, double, double, double); // prf5

typedef unsigned char  UCHAR;            // uch
typedef unsigned long  ULONG;            // ul
typedef unsigned long  DWORD;            // dw
typedef unsigned long *PDWORD;           // pdw

#if defined(UNIX)
typedef unsigned int   WORD;             // w
typedef unsigned int  *PWORD;            // pw
#define O_BINARY 0x0000
#endif

#endif

