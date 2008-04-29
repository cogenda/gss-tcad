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
/*  Last update: Sep 5, 2006                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

/*  We need two file streams for message output: console and log file.       */
/*  Output Device Select                                                     */

#include "log.h"
//use C++ stream system
GSS_LOG_STREAM gss_log;

//support the old c fprintf front end.
char log_buf[128];
void GSS_LOG()
{
  gss_log.string_buf()<<log_buf;
  gss_log.record();
}
