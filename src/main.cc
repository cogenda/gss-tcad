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
/*  Last update: Oct 15, 2007                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/

#include "petsc.h"
#include "petsctime.h"
#include "log.h"
#include "ctrl.h"


int main(int argc, char ** args)
{
  PetscInitialize(&argc,&args,PETSC_NULL,PETSC_NULL);
  
  PetscLogDouble v1,v2,elapsed_time;
  PetscGetTime(&v1);
    
  cout<<"****************************************************************************"<<endl;
  cout<<"*              8888888         88888888         88888888                   *"<<endl;
  cout<<"*            8                8                8                           *"<<endl;
  cout<<"*           8                 8                8                           *"<<endl;
  cout<<"*           8                  88888888         88888888                   *"<<endl;
  cout<<"*           8      8888                8                8                  *"<<endl;
  cout<<"*            8       8                 8                8                  *"<<endl;
  cout<<"*              888888         888888888        888888888                   *"<<endl;
  cout<<"*                                                                          *"<<endl;
  cout<<"*       A Two-Dimensional General-purpose Semiconductor Simulator.         *"<<endl;
  cout<<"*                                                                          *"<<endl;
  if(sizeof(PetscScalar)==sizeof(double))
  cout<<"* This is GSS version 0.46.08 with double precision.                       *"<<endl;
  if(sizeof(PetscScalar)==sizeof(long double))
  cout<<"* This is GSS version 0.46.08 with long double precision.                  *"<<endl;
  cout<<"* Copyright (C) 2005-2008 by Gong Ding and Zhang XiangHua.                 *"<<endl;
  cout<<"* Please send ANY comments, corrections or suggestions to gdiso@ustc.edu.  *"<<endl;
  cout<<"****************************************************************************"<<endl<<endl;

  if(argc<2)
  {
    cout<<"usage: gss card_file [petsc_options]"<<endl;
    cout<<"For more information, please refer to GSS User's Guide."<<endl<<endl;
    exit(0);
  }

  CmdBuf *cmdbuf= new CmdBuf;
  if(cmdbuf->parse_file(args[1]))
  {
    delete cmdbuf;
    exit(0);
  }

  char *logfile = new char[strlen(args[1])+16];
  sprintf(logfile,"%s.log",args[1]);
  gss_log.GSS_LOG_BEGIN(logfile);

  SolveControl *solver= new SolveControl;
  if(solver->static_config(cmdbuf)) exit(0);
  if(solver->run_control(cmdbuf->cmd_list)) exit(0);

  delete solver;
  delete cmdbuf;
  delete [] logfile;

  gss_log.GSS_LOG_END();
  
  PetscGetTime(&v2);
  elapsed_time = v2 - v1; 
  cout<<"GSS finished. Total time is "<<elapsed_time<<" seconds. Good bye!"<<endl;
  
  PetscFinalize();
  return 0;
}
