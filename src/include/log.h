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
/*  Last update: Sep 5, 2006                                                 */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

/*  We need two streams for message output: screen and log file.             */
/*  Output Device Select                                                     */

#ifndef _log_h_
#define _log_h_

#include <iostream>  
#include <fstream> 
#include <sstream>
#include <string>

using namespace std;

extern char log_buf[128];
extern void GSS_LOG();

class GSS_LOG_STREAM
{
  private:
    stringstream str_buf;    //the string buf for out stream.
    ofstream     log_file;   //the file stream
  public:  		
  void GSS_LOG_BEGIN(char *logfile)   { log_file.open(logfile,ofstream::out | ofstream::trunc); }
  void GSS_LOG_END()                  { log_file.close(); }
  stringstream & string_buf()         { return str_buf;}
  
  void record()
  {
     cout<<str_buf.str();
     cout.flush();
     log_file<<str_buf.str();
     str_buf.str("");
  }

  GSS_LOG_STREAM()
  {
    str_buf.setf(ios::scientific);
  }	
};

extern GSS_LOG_STREAM gss_log;

#endif
