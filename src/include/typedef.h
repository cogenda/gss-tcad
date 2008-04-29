/*****************************************************************************/
/*   	        8888888         88888888         88888888                    */
/*  	      8                8                8                            */
/* 	     8                 8                8                            */
/*  	     8                  88888888         88888888                    */
/* 	     8      8888                8                8                   */
/* 	      8       8                 8                8                   */
/* 	        888888         888888888        888888888                    */
/*                                                                           */
/*       A Two-Dimensional General Purpose Semiconductor Simulator.          */
/*                                                                           */
/*  GSS 0.4x                                                                 */
/*  Last update: Sep 5, 2005                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

/*                                                                           */
/*                       global typedef header                               */
/*                                                                           */

#ifndef _typedef_h_
#define _typedef_h_


/******************************************************************************/
/****************************  TYPE DEFINITIONS	*******************************/
/******************************************************************************/

typedef unsigned int UINT;

const int  ERROR = 1;
/*
       the default precision for float point arithmetic is double
*/

#define DPREC
typedef double	     REAL;

// PETSC_SCALAR is double precision
#define DP_PETSC_SCALAR

// PETSC_SCALAR is long double precision
#define LP_PETSC_SCALAR

#endif
