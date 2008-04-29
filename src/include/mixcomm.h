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
/*  Last update: May 09, 2007                                                */
/*                                                                           */
/*  Gong Ding   Xuan Chun                                                    */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#ifndef _mixcomm_h_
#define _mixcomm_h_

/* network function */
#include <errno.h>
#include <netinet/in.h> /* IPv4 socket address structres. */
#include <netdb.h> /* Access to DNS lookups. */
#include <arpa/inet.h> /* inet_ntop function. */
#include <sys/socket.h> /* Socket functions. */
#include <setjmp.h>
#include <signal.h>


#ifdef CYGWIN
#define MSG_FLAG 0
#else
#define MSG_FLAG MSG_WAITALL
#endif

typedef struct
{
  PetscScalar pdI_pdV;
  Vec         pdI_pdw;
  Vec         pdF_pdV;
  Vec         pdw_pdV;
} sPINcond;  
extern sigjmp_buf net_buf;
extern int NDEVServer(int port,int &listener, int &client);
#endif
