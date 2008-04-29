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
/*  Last update: Feb 19, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#include "em_fem.h"
#include "zonedata.h"
#include "log.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/* ----------------------------------------------------------------------------
 * EM_FEM_Solver----------------by zhangxih 
 */
/* ----------------------------------------------------------------------------
 * EM_FEM_Solver::constructor
 */
EM_FEM_Solver::EM_FEM_Solver()
{
  N=0;
  its  = 0;
  norm = 0.0;
}

/* ----------------------------------------------------------------------------
 * EM_FEM_Solver::init_solver:  This function do initial setup for FEM solver
 */
int EM_FEM_Solver::init_solver(SolveDefine &sv)
{
  //search for PHOTOGEN statement
  int pg_flag=0;
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();pcmdbuf->goto_next_cmd())
  {
    if(pcmdbuf->is_current_cmd("PHOTOGEN"))
    {
      list<Cmd>::iterator    pcmd = pcmdbuf->get_current_cmd();

      //get incidence light's wave length from input file, else its value is 0.532 micron
      lamda = pcmd->get_number("wavelen",0,0.532)*scale_unit.s_micron;

      //get incident power of incidence light from input file, else its value is 0.0
      power = pcmd->get_number("intensity",0,1.0)*(scale_unit.s_joule/scale_unit.s_second)
              /(scale_unit.s_centimeter*scale_unit.s_centimeter);

      // get weight of TE mode in incidence light
      WeightOfTE=pcmd->get_number("wte",0,0.0);
      //get weight of TM mode in incidence light
      WeightOfTM=pcmd->get_number("wtm",0,1.0);
      if(fabs((WeightOfTE+WeightOfTM)-1.0)>=1e-3)
      {
        gss_log.string_buf()<<"PHOTOGEN: Error! The sum of TM mode Weight and TE mode Weight must be 1.0\n";
        gss_log.record();
        return 1;
      }

      //get incident angle of incidence light from input file, else its value is 90.0
      OptAngle = pcmd->get_number("angle",0,90.0);
      if (OptAngle<0||OptAngle>180)
      {
        gss_log.string_buf()<<"PHOTOGEN: Error! Incident angle exceeding allowed range 0<theta<180\n";
        gss_log.record();
        return 1;
      }
      OptTheta=OptAngle*PI/180;
      //if the optG should be added to original value instead of overwrite it, for multi-freq wave calculation
      append = pcmd->get_bool("append",0,false);

      //get quantum efficiency of semiconductor material
      eta = pcmd->get_number("quan.eff",0,1.0);
      if (eta<=0||eta>1)
      {
        gss_log.string_buf()<<"PHOTOGEN: Error! Quantum efficiency exceeding allowed range 0<=eff<=1.0\n";
        gss_log.record();
        return 1;
      }

      //get phase differnce between TE wave and TM wave
      InitPhi = pcmd->get_number("phase.diff",0,0.0);
      if(fabs(InitPhi)>360)
      {
        gss_log.string_buf()<<"PHOTOGEN: Error! Phase difference exceeding allowed range 0<phi<360 \n";
        gss_log.record();
        return 1;
      }
      PhiDiff = InitPhi*PI/180;
      
      //this PHOTOGEN statement is parsed, delete it
      pg_flag=1;
      pcmdbuf->delete_current_cmd();
      break;
    }
  }

  if(!pg_flag)
  {
    gss_log.string_buf()<<"NO suitable PHOTOGEN statement found! Each EMFEM solver requires unique one.\n";
    gss_log.record();
    return 1;
  }
  
  for(int i=0;i<3;i++)//initializing element matrix
  {
    for(int j=0;j<3;j++)
      eleMatrix[i][j]=0;
    eleV[i]=0;
    eleNode[i]=0;
  }

  //fprintf(SCREEN,"FEM solver init...\n");
  gss_log.string_buf()<<"FEM solver init...\n";
  gss_log.record();
  gss_log.string_buf()<<"Linear Solver "<<sv.LS<<" is loaded.\n\n";
  gss_log.record(); 
  // compute the scale of problem
  N=gnode.size();

  //ierr=PetscOptionsSetValue("-options_left",PETSC_NULL);
  ierr=VecCreateSeq(PETSC_COMM_SELF,2*N,&x);//create vector and matrix in solution this problem
  ierr=VecDuplicate(x,&b);
  MatCreate(PETSC_COMM_SELF,&A);
  MatSetSizes(A,2*N,2*N,2*N,2*N);
  if(sv.LS=="superlu")
  {
#ifdef PETSC_HAVE_SUPERLU  
    MatSetType(A,MATSUPERLU);
#else
    MatSetType(A,MATSEQAIJ);    
#endif    
  }
  else if(sv.LS=="umfpack")
  {
#ifdef PETSC_HAVE_UMFPACK
    MatSetType(A,MATUMFPACK);
#else
    MatSetType(A,MATSEQAIJ);    
#endif    
  }
  else
  {
    MatSetType(A,MATSEQAIJ);
  }
  MatSeqAIJSetPreallocation(A,18,0);
  
  ierr=KSPCreate(PETSC_COMM_SELF,&solver);
  ierr=KSPSetType(solver,KSPPREONLY);
  //ierr=KSPSetInitialGuessNonzero(solver,PETSC_TRUE);
  ierr=KSPGetPC(solver,&prec);
  if(sv.LS=="lu"||sv.LS=="superlu"||sv.LS=="umfpack")
  {
    KSPSetType(solver,KSPPREONLY);
    PCSetType(prec,PCLU);
    PCFactorSetPivoting(prec,1.0);
    PCFactorReorderForNonzeroDiagonal(prec,1e-14);
    PCFactorSetShiftNonzero(prec,1e-12);
  }
  else
  {
    KSPSetType(solver,sv.LS.c_str());
    if(sv.LS=="gmres")  KSPGMRESSetRestart(solver,150);
    KSPSetTolerances(solver,1e-10*N,1e-20*N,PETSC_DEFAULT,N/10);
    PCSetType(prec,PCILU);
    PCFactorSetLevels(prec,3);
    PCFactorSetShiftNonzero(prec,1e-14);
    PCFactorReorderForNonzeroDiagonal(prec,1e-12);
  }
  ierr=KSPSetFromOptions(solver);

  TMEZ.resize(2*N);
  TMHX.resize(2*N);
  TMHY.resize(2*N);
  TEHZ.resize(2*N);
  TEEX.resize(2*N);
  TEEY.resize(2*N);

  //initializing Maximum and Minimum coordinate values
  XLeft   =  1.0e9;
  XRight  = -1.0e9;
  YTop    = -1.0e9;
  YBottom =  1.0e9;

  //Find Maximum and Minimum coordinate values
  for(int WholeNode=0;WholeNode<gnode.size();WholeNode++)
  {
    if (gnode[WholeNode].x<=XLeft)  XLeft = gnode[WholeNode].x;
    else if (gnode[WholeNode].x>=XRight)  XRight = gnode[WholeNode].x;
    if (gnode[WholeNode].y>=YTop)  YTop = gnode[WholeNode].y;
    else if (gnode[WholeNode].y<=YBottom)  YBottom = gnode[WholeNode].y;
  }

  //initializing 4 corners of vacuum zone:
  VacXMax=-1.0e9;   //give enough small value to VacXMax and VacYMax
  VacYMax=-1.0e9;
  VacXMin=1.0e9;  //give enough large value to VacXMin and VacYMin
  VacYMin=1.0e9;
  //find 4 corners of vacuum zone
  for(int i=0;i<gtri.size();i++)
  {
    int ZoneSerial=gtri[i].zone_index;
    if(zonedata[ZoneSerial]->material_type==Vacuum)//for vacuum zone
    {
      for (int CorFind=0;CorFind<3;CorFind++)
      {
        PetscScalar TriNodeX=gnode[gtri[i].g_node[CorFind]].x;
        PetscScalar TriNodeY=gnode[gtri[i].g_node[CorFind]].y;
        if (TriNodeX>=VacXMax)  VacXMax=TriNodeX;
        if (TriNodeX<=VacXMin)  VacXMin=TriNodeX;
        if (TriNodeY>=VacYMax)  VacYMax=TriNodeY;
        if (TriNodeY<=VacYMin)  VacYMin=TriNodeY;
      }
    }
  }
  
  //Just for Source Calculating
  EpsInVac = 8.854187818e-12*scale_unit.s_coulomb/scale_unit.s_volt/scale_unit.s_meter;
  MuInVac  = 12.56637061e-7*pow(scale_unit.s_second,2)/scale_unit.s_coulomb*scale_unit.s_volt/scale_unit.s_meter;
  WaveVector0 = 2*PI/lamda; //lamda input by users
  SourceTM = pow(WaveVector0,2)*sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac));
  SourceTE = pow(WaveVector0,2)*sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac));

  // constant for UPML
  m=2;
  R=1.0e-3;

  return 0;
}

/* ----------------------------------------------------------------------------
 * EM_FEM_Solver::do_solve:  This function solve the problem
 */
int EM_FEM_Solver::do_solve(SolveDefine &sv)
{
  if(WeightOfTM)
  {
    gss_log.string_buf()<<"Begin assembling TM mode FEM matrix...\n";
    gss_log.record();

    //assembling the stiff matrix---for TM mode
    for(int i=0;i<gtri.size();i++)
    {
      int ZoneSerial=gtri[i].zone_index;

      eleB[0]=gnode[gtri[i].g_node[1]].y-gnode[gtri[i].g_node[2]].y;
      eleB[1]=gnode[gtri[i].g_node[2]].y-gnode[gtri[i].g_node[0]].y;
      eleB[2]=gnode[gtri[i].g_node[0]].y-gnode[gtri[i].g_node[1]].y;
      eleC[0]=gnode[gtri[i].g_node[2]].x-gnode[gtri[i].g_node[1]].x;
      eleC[1]=gnode[gtri[i].g_node[0]].x-gnode[gtri[i].g_node[2]].x;
      eleC[2]=gnode[gtri[i].g_node[1]].x-gnode[gtri[i].g_node[0]].x;

      if(zonedata[ZoneSerial]->material_type==Semiconductor)//semiconductor zone
      {
        SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[ZoneSerial]);

        PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
        PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
        RealPermi = n*n-k*k;
        ImagePermi = 2*n*k;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta));

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/6;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/12;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=-gtri[i].area*SourceTM*((1-RealPermi)*cos(SourceAlpha)-ImagePermi*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=-gtri[i].area*SourceTM*(ImagePermi*cos(SourceAlpha)+(1-RealPermi)*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==Insulator)//insulator zone
      {
        ISZone *pzonedata= dynamic_cast< ISZone * >(zonedata[ZoneSerial]);

        PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
        PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
        RealPermi = n*n-k*k;
        ImagePermi = 2*n*k;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta));

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/6;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/12;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=-gtri[i].area*SourceTM*((1-RealPermi)*cos(SourceAlpha)-ImagePermi*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=-gtri[i].area*SourceTM*(ImagePermi*cos(SourceAlpha)+(1-RealPermi)*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==Conductor)//Conductor zone
      {
        ElZone *pzonedata = dynamic_cast<ElZone *>(zonedata[ZoneSerial]);

        PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
        PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
        RealPermi = n*n-k*k;
        ImagePermi = 2*n*k;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta));

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/6;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/12;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=-gtri[i].area*SourceTM*((1-RealPermi)*cos(SourceAlpha)-ImagePermi*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=-gtri[i].area*SourceTM*(ImagePermi*cos(SourceAlpha)+(1-RealPermi)*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==Vacuum)//for vacuum zone
      {
        VacuumZone *pzonedata= dynamic_cast< VacuumZone * >(zonedata[ZoneSerial]);

        RealPermi = 1.0;
        ImagePermi = 0.0;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta));

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/6;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*RealPermi*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=gtri[i].area*ImagePermi*pow(WaveVector0,2)/12;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=-gtri[i].area*SourceTM*((1-RealPermi)*cos(SourceAlpha)-ImagePermi*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=-gtri[i].area*SourceTM*(ImagePermi*cos(SourceAlpha)+(1-RealPermi)*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==PML)//for PML zone
      {
        PMLZone *pzonedata= dynamic_cast< PMLZone * >(zonedata[ZoneSerial]);

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          //transition variable
          TransVar = ((m+1)*log(1/R)*lamda)/(4*PI*(VacXMin-XLeft));
          CurrentX = gnode[MatJi].x;
          CurrentY = gnode[MatJi].y;
          // give Variable values in different region
          if (gtri[i].xc<=VacXMin&&gtri[i].yc>=VacYMax)  //left-top corner
          {
            PetscScalar CurrentP = pow((VacXMin-CurrentX)/(VacXMin-XLeft),m);
            PetscScalar CurrentQ = pow((CurrentY-VacYMax)/(YTop-VacYMax),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
            //printf("Value of TransVar in function 1: %e\n",TransVar); //this place is ok
          }
          else if (gtri[i].xc<VacXMin&&gtri[i].yc<VacYMax&&gtri[i].yc>VacYMin)  //left part
          {
            PetscScalar CurrentP = pow((VacXMin-CurrentX)/(VacXMin-XLeft),m);
            //PetscScalar CurrentQ = ;
            RealAlphaX = 1/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = 1;
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = TransVar*CurrentP/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = -TransVar*CurrentP;
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentP*TransVar);
            //printf("Value of TransVar in function 2: %e\n",TransVar);  //this place is ok
          }
          else if (gtri[i].xc<=VacXMin&&gtri[i].yc<=VacYMin)  //left-bottom corner
          {
            PetscScalar CurrentP = pow((VacXMin-CurrentX)/(VacXMin-XLeft),m);
            PetscScalar CurrentQ = pow((VacYMin-CurrentY)/(VacYMin-YBottom),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
          }
          else if (gtri[i].xc>VacXMin&&gtri[i].xc<VacXMax&&gtri[i].yc<VacYMin)  //bottom part
          {
            //PetscScalar CurrentP = ;
            PetscScalar CurrentQ = pow((VacYMin-CurrentY)/(VacYMin-YBottom),m);
            RealAlphaX = 1;
            RealAlphaY = 1/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = -TransVar*CurrentQ;
            ImageAlphaY = TransVar*CurrentQ/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentQ*TransVar);
          }
          else if (gtri[i].xc>=VacXMax&&gtri[i].yc<=VacYMin)  //right-bottom corner
          {
            PetscScalar CurrentP = pow((CurrentX-VacXMax)/(XRight-VacXMax),m);
            PetscScalar CurrentQ = pow((VacYMin-CurrentY)/(VacYMin-YBottom),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
          }
          else if (gtri[i].xc>VacXMax && gtri[i].yc>VacYMin && gtri[i].yc<VacYMax)  //right part
          {
            PetscScalar CurrentP = pow((CurrentX-VacXMax)/(XRight-VacXMax),m);
            //PetscScalar CurrentQ = ;
            RealAlphaX = 1/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = 1;
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = TransVar*CurrentP/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = -TransVar*CurrentP;
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentP*TransVar);
          }
          else if (gtri[i].xc>=VacXMax&&gtri[i].yc>=VacYMax)  //right-top corner
          {
            PetscScalar CurrentP = pow((CurrentX-VacXMax)/(XRight-VacXMax),m);
            PetscScalar CurrentQ = pow((CurrentY-VacYMax)/(YTop-VacYMax),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
          }
          //else if (gtri[i].xc>VacXMin&&gtri[i].xc<VacXMax&&gtri[i].yc>VacYMax)  //top part
          else                                                                    //top part
          {
            //PetscScalar CurrentP = ;
            PetscScalar CurrentQ = pow((CurrentY-VacYMax)/(YTop-VacYMax),m);
            RealAlphaX = 1;
            RealAlphaY = 1/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = -TransVar*CurrentQ;
            ImageAlphaY = TransVar*CurrentQ/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentQ*TransVar);
          }

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(RealAlphaX*eleB[j]*eleB[k]+RealAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*RealBeta/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImageAlphaX*eleB[j]*eleB[k]+ImageAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*ImageBeta/6;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(RealAlphaX*eleB[j]*eleB[k]+RealAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*RealBeta/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImageAlphaX*eleB[j]*eleB[k]+ImageAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*ImageBeta/12;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          v=0.0;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }
    }

    //Boundary outside UPML is metal boundary, which is Dirichlet Boundary for TM mode and should be deal seperately
    for(int i=0;i<gnode.size();i++)
    {
      PetscScalar ErroLimit=1e-10;

      if(fabs(gnode[i].x-XLeft)<=ErroLimit||fabs(gnode[i].x-XRight)<=ErroLimit||fabs(gnode[i].y-YBottom)<=ErroLimit||fabs(gnode[i].y-YTop)<=ErroLimit)
      {
        PetscScalar UnLimit=1e50; //define infinite value
        int MatJiPosOne=2*i;
        int MatJiPosTwo=2*i+1;

        v=UnLimit; //for stiff matrix
        ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatJiPosOne,&v,ADD_VALUES);
        ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatJiPosTwo,&v,ADD_VALUES);

        v=0.0*UnLimit; //for right hand vector
        ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
        ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
      }
    }

    ierr=MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr=MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    ierr=VecAssemblyBegin(b);
    ierr=VecAssemblyEnd(b);

    gss_log.string_buf()<<"TM mode FEM matrix assembling finished.\n";
    gss_log.record();
    //MatView(A,PETSC_VIEWER_DRAW_WORLD);
    //getchar();
    
    ierr=KSPSetOperators(solver,A,A,DIFFERENT_NONZERO_PATTERN);
    ierr=KSPSetUp(solver);

    ierr=KSPSolve(solver,b,x);
    ierr=KSPGetConvergedReason(solver,&reason);

    if (reason<0)
    {
      gss_log.string_buf()<<"Solution divergence! You may need to change your linear solver.\n";
      gss_log.record();
      return 1;
    }
    else
    {
      ierr=KSPGetIterationNumber(solver,&its);
      gss_log.string_buf()<<"Solution Converged in "<<its<<" iterations.\n";
      gss_log.record();
    }

    //output the result
    VecGetArray(x,&TM);

    gss_log.string_buf()<<"TM mode solve finished.\n\n";
    gss_log.record();

    //update the solution for plot
    TM_solution_update(sv);

    //zero the matrix and vector to use them again before destroy
    MatZeroEntries(A);
    VecZeroEntries(b);
    VecZeroEntries(x);
  }

  //for TE mode
  if(WeightOfTE)
  {
    gss_log.string_buf()<<"Begin assembling TE mode FEM matrix...\n";
    gss_log.record();

    //assembling the stiff matrix---for TE mode
    for(int i=0;i<gtri.size();i++)
    {
      int ZoneSerial=gtri[i].zone_index;

      eleB[0]=gnode[gtri[i].g_node[1]].y-gnode[gtri[i].g_node[2]].y;
      eleB[1]=gnode[gtri[i].g_node[2]].y-gnode[gtri[i].g_node[0]].y;
      eleB[2]=gnode[gtri[i].g_node[0]].y-gnode[gtri[i].g_node[1]].y;
      eleC[0]=gnode[gtri[i].g_node[2]].x-gnode[gtri[i].g_node[1]].x;
      eleC[1]=gnode[gtri[i].g_node[0]].x-gnode[gtri[i].g_node[2]].x;
      eleC[2]=gnode[gtri[i].g_node[1]].x-gnode[gtri[i].g_node[0]].x;

      if(zonedata[ZoneSerial]->material_type==Semiconductor)//semiconductor zone
      {
        SMCZone *pzonedata= dynamic_cast< SMCZone * >(zonedata[ZoneSerial]);
        PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
        PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
        RealPermi = n*n-k*k;
        ImagePermi = 2*n*k;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta))-PhiDiff;

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=gtri[i].area*SourceTE*((1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*cos(SourceAlpha)
                                   +ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=gtri[i].area*SourceTE*(-ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*cos(SourceAlpha)
                                   +(1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==Insulator)//insulator zone
      {
        ISZone *pzonedata= dynamic_cast< ISZone * >(zonedata[ZoneSerial]);

        PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
        PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
        RealPermi = n*n-k*k;
        ImagePermi = 2*n*k;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta))-PhiDiff;

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=gtri[i].area*SourceTE*((1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*cos(SourceAlpha)
                                   +ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=gtri[i].area*SourceTE*(-ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*cos(SourceAlpha)
                                   +(1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==Conductor)//Conductor zone
      {
        ElZone *pzonedata = dynamic_cast<ElZone *>(zonedata[ZoneSerial]);

        PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
        PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
        RealPermi = n*n-k*k;
        ImagePermi = 2*n*k;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta))-PhiDiff;

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=gtri[i].area*SourceTE*((1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*cos(SourceAlpha)
                                   +ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=gtri[i].area*SourceTE*(-ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*cos(SourceAlpha)
                                   +(1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==Vacuum)//for vacuum zone
      {
        VacuumZone *pzonedata= dynamic_cast< VacuumZone * >(zonedata[ZoneSerial]);

        RealPermi = 1.0;
        ImagePermi = 0.0;

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          SourceAlpha = -WaveVector0*(gnode[MatJi].x*cos(OptTheta)-gnode[MatJi].y*sin(OptTheta))-PhiDiff;

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area)
                              -gtri[i].area*pow(WaveVector0,2)/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*(eleB[j]*eleB[k]+eleC[j]*eleC[k])/(4*gtri[i].area);
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          //real part
          v=gtri[i].area*SourceTE*((1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*cos(SourceAlpha)
                                   +ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          //image part
          v=gtri[i].area*SourceTE*(-ImagePermi/(pow(RealPermi,2)+pow(ImagePermi,2))*cos(SourceAlpha)
                                   +(1-RealPermi/(pow(RealPermi,2)+pow(ImagePermi,2)))*sin(SourceAlpha))/3;
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }

      else if(zonedata[ZoneSerial]->material_type==PML)//for PML zone
      {
        PMLZone *pzonedata= dynamic_cast< PMLZone * >(zonedata[ZoneSerial]);

        for(int j=0;j<3;j++)
        {
          int MatJi=gtri[i].g_node[j];
          int MatJiPosOne=2*MatJi;
          int MatJiPosTwo=2*MatJi+1;

          //transition variable
          TransVar = ((m+1)*log(1/R)*lamda)/(4*PI*(VacXMin-XLeft));
          CurrentX = gnode[MatJi].x;
          CurrentY = gnode[MatJi].y;
          // give Variable values in different region
          if (gtri[i].xc<=VacXMin&&gtri[i].yc>=VacYMax)  //left-top corner
          {
            PetscScalar CurrentP = pow((VacXMin-CurrentX)/(VacXMin-XLeft),m);
            PetscScalar CurrentQ = pow((CurrentY-VacYMax)/(YTop-VacYMax),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
          }
          else if (gtri[i].xc<VacXMin&&gtri[i].yc<VacYMax&&gtri[i].yc>VacYMin)  //left part
          {
            PetscScalar CurrentP = pow((VacXMin-CurrentX)/(VacXMin-XLeft),m);
            //PetscScalar CurrentQ = ;
            RealAlphaX = 1/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = 1;
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = TransVar*CurrentP/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = -TransVar*CurrentP;
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentP*TransVar);
          }
          else if (gtri[i].xc<=VacXMin&&gtri[i].yc<=VacYMin)  //left-bottom corner
          {
            PetscScalar CurrentP = pow((VacXMin-CurrentX)/(VacXMin-XLeft),m);
            PetscScalar CurrentQ = pow((VacYMin-CurrentY)/(VacYMin-YBottom),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
          }
          else if (gtri[i].xc>VacXMin&&gtri[i].xc<VacXMax&&gtri[i].yc<VacYMin)  //bottom part
          {
            //PetscScalar CurrentP = ;
            PetscScalar CurrentQ = pow((VacYMin-CurrentY)/(VacYMin-YBottom),m);
            RealAlphaX = 1;
            RealAlphaY = 1/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = -TransVar*CurrentQ;
            ImageAlphaY = TransVar*CurrentQ/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentQ*TransVar);
          }
          else if (gtri[i].xc>=VacXMax&&gtri[i].yc<=VacYMin)  //right-bottom corner
          {
            PetscScalar CurrentP = pow((CurrentX-VacXMax)/(XRight-VacXMax),m);
            PetscScalar CurrentQ = pow((VacYMin-CurrentY)/(VacYMin-YBottom),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
          }
          else if (gtri[i].xc>VacXMax && gtri[i].yc>VacYMin && gtri[i].yc<VacYMax)  //right part
          {
            PetscScalar CurrentP = pow((CurrentX-VacXMax)/(XRight-VacXMax),m);
            //PetscScalar CurrentQ = ;
            RealAlphaX = 1/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = 1;
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = TransVar*CurrentP/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = -TransVar*CurrentP;
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentP*TransVar);
          }
          else if (gtri[i].xc>=VacXMax&&gtri[i].yc>=VacYMax)  //right-top corner
          {
            PetscScalar CurrentP = pow((CurrentX-VacXMax)/(XRight-VacXMax),m);
            PetscScalar CurrentQ = pow((CurrentY-VacYMax)/(YTop-VacYMax),m);
            RealAlphaX = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentP,2));
            RealAlphaY = (1+pow(TransVar,2)*CurrentP*CurrentQ)/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*(1-CurrentP*CurrentQ*pow(TransVar,2));
            ImageAlphaX = ((CurrentP-CurrentQ)*TransVar)/(1+pow(TransVar*CurrentP,2));
            ImageAlphaY = ((CurrentQ-CurrentP)*TransVar)/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-(CurrentP+CurrentQ)*TransVar);
          }
          else                                                                    //top part
          {
            //PetscScalar CurrentP = ;
            PetscScalar CurrentQ = pow((CurrentY-VacYMax)/(YTop-VacYMax),m);
            RealAlphaX = 1;
            RealAlphaY = 1/(1+pow(TransVar*CurrentQ,2));
            RealBeta   = (-pow((2*PI/lamda),2))*1;
            ImageAlphaX = -TransVar*CurrentQ;
            ImageAlphaY = TransVar*CurrentQ/(1+pow(TransVar*CurrentQ,2));
            ImageBeta   = (-pow((2*PI/lamda),2))*(-CurrentQ*TransVar);
          }

          for(int k=0;k<3;k++)
          {
            int MatKi=gtri[i].g_node[k];
            int MatKiPosOne=2*MatKi;
            int MatKiPosTwo=2*MatKi+1;

            if(k==j)
            {
              eleMatrix[j][k]=(RealAlphaX*eleB[j]*eleB[k]+RealAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*RealBeta/6;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImageAlphaX*eleB[j]*eleB[k]+ImageAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*ImageBeta/6;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
            else
            {
              eleMatrix[j][k]=(RealAlphaX*eleB[j]*eleB[k]+RealAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*RealBeta/12;
              v=eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosOne,&v,ADD_VALUES);
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosTwo,&v,ADD_VALUES);

              eleMatrix[j][k]=(ImageAlphaX*eleB[j]*eleB[k]+ImageAlphaY*eleC[j]*eleC[k])/(4*gtri[i].area)
                              +gtri[i].area*ImageBeta/12;
              v=-eleMatrix[j][k];
              ierr=MatSetValues(A,1,&MatJiPosOne,1,&MatKiPosTwo,&v,ADD_VALUES);
              v=-v;
              ierr=MatSetValues(A,1,&MatJiPosTwo,1,&MatKiPosOne,&v,ADD_VALUES);
            }
          }
          v=0.0;
          ierr=VecSetValues(b,1,&MatJiPosOne,&v,ADD_VALUES);
          ierr=VecSetValues(b,1,&MatJiPosTwo,&v,ADD_VALUES);
        }
      }
    }

    ierr=MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    ierr=MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    ierr=VecAssemblyBegin(b);
    ierr=VecAssemblyEnd(b);
    //MatView(A,PETSC_VIEWER_DRAW_WORLD);
    //getchar();
    gss_log.string_buf()<<"TE mode FEM matrix assembling finished.\n";
    gss_log.record();
    
    ierr=KSPSetOperators(solver,A,A,DIFFERENT_NONZERO_PATTERN);
    ierr=KSPSetUp(solver);
    ierr=KSPSolve(solver,b,x);
    ierr=KSPGetConvergedReason(solver,&reason);
    
    if (reason<0)
    {
      gss_log.string_buf()<<"Solution divergence! You may need to change your linear solver.\n";
      gss_log.record();
      return 1;
    }
    else
    {
      ierr=KSPGetIterationNumber(solver,&its);
      gss_log.string_buf()<<"Solution Converged in "<<its<<" iterations.\n";
      gss_log.record();
    }

    //output the result
    VecGetArray(x,&TE);

    gss_log.string_buf()<<"TE mode solve finished.\n\n";
    gss_log.record();

    //update the solution for plot
    TE_solution_update(sv);
  }

  //calculating other field components in both modes
  field_calculation(sv);

  //calculating optical generation rate
  generation_cal(sv);

  gss_log.string_buf()<<"Carrier generation ratio calculating finished.\n\n";
  gss_log.record();

  return 0;
}


/* ------------------------------------------------------------------------------------
 * EM_FEM_Solver::TM_solution_update:  This function do update the result of TM mode
 */
int EM_FEM_Solver::TM_solution_update(SolveDefine &sv)
{
  int offset=0;
  for(int i=0;i<gtri.size();i++)//update by triangle
  {
    int ZoneSerial=gtri[i].zone_index;
    if(zonedata[ZoneSerial]->material_type==Semiconductor)//semiconductor zone
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        PetscScalar imageInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        pzonedata->aux[localIndex].OpEz = sqrt(pow((TM[2*offset]+realInc),2)+pow((TM[2*offset+1]+imageInc),2));
        TMEZ[2*offset] = TM[2*offset]+realInc;
        TMEZ[2*offset+1] = TM[2*offset+1]+imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==Insulator)//insulator zone
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        PetscScalar imageInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        pzonedata->aux[localIndex].OpEz = sqrt(pow((TM[2*offset]+realInc),2)+pow((TM[2*offset+1]+imageInc),2));
        TMEZ[2*offset] = TM[2*offset]+realInc;
        TMEZ[2*offset+1] = TM[2*offset+1]+imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==Conductor)//conductor zone
    {

      ElZone *pzonedata =dynamic_cast< ElZone *>(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        PetscScalar imageInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        pzonedata->aux[localIndex].OpEz = sqrt(pow((TM[2*offset]+realInc),2)+pow((TM[2*offset+1]+imageInc),2));
        TMEZ[2*offset] = TM[2*offset]+realInc;
        TMEZ[2*offset+1] = TM[2*offset+1]+imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==Vacuum)//vacuum zone
    {
      VacuumZone *pzonedata = dynamic_cast< VacuumZone *>(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        PetscScalar imageInc=sqrt(WeightOfTM*power*sqrt(MuInVac/EpsInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)-gnode[offset].y*sin(OptTheta)));
        pzonedata->aux[localIndex].OpEz = sqrt(pow((TM[2*offset]+realInc),2)+pow((TM[2*offset+1]+imageInc),2));
        TMEZ[2*offset] = TM[2*offset]+realInc;
        TMEZ[2*offset+1] = TM[2*offset+1]+imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==PML)//PML zone
    {
      PMLZone *pzonedata = dynamic_cast< PMLZone *>(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        pzonedata->aux[localIndex].OpEz = sqrt(pow(TM[2*offset],2)+pow(TM[2*offset+1],2));
        TMEZ[2*offset] = TM[2*offset];
        TMEZ[2*offset+1] = TM[2*offset+1];
      }
    }
  }
  return 0;
}


/* ----------------------------------------------------------------------------------
 * EM_FEM_Solver::TE_solution_update:  This function do update the result of TE mode
 */
int EM_FEM_Solver::TE_solution_update(SolveDefine &sv)
{
  int offset=0;
  for(int i=0;i<gtri.size();i++)//update by triangle
  {
    int ZoneSerial=gtri[i].zone_index;
    if(zonedata[ZoneSerial]->material_type==Semiconductor)//semiconductor zone
    {
      SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                            -gnode[offset].y*sin(OptTheta))-PhiDiff);
        PetscScalar imageInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                             -gnode[offset].y*sin(OptTheta))-PhiDiff);
        pzonedata->aux[localIndex].OpHz = sqrt(pow((TE[2*offset]+realInc),2)+pow((TE[2*offset+1]+imageInc),2));
        TEHZ[2*offset] = TE[2*offset]+realInc;
        TEHZ[2*offset+1] = TE[2*offset+1] +imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==Insulator)//insulator zone
    {
      ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                            -gnode[offset].y*sin(OptTheta))-PhiDiff);
        PetscScalar imageInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                             -gnode[offset].y*sin(OptTheta))-PhiDiff);
        pzonedata->aux[localIndex].OpHz = sqrt(pow((TE[2*offset]+realInc),2)+pow((TE[2*offset+1]+imageInc),2));
        TEHZ[2*offset] = TE[2*offset]+realInc;
        TEHZ[2*offset+1] = TE[2*offset+1] +imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==Conductor)//conductor zone
    {

      ElZone *pzonedata =dynamic_cast< ElZone *>(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                            -gnode[offset].y*sin(OptTheta))-PhiDiff);
        PetscScalar imageInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                             -gnode[offset].y*sin(OptTheta))-PhiDiff);
        pzonedata->aux[localIndex].OpHz = sqrt(pow((TE[2*offset]+realInc),2)+pow((TE[2*offset+1]+imageInc),2));
        TEHZ[2*offset] = TE[2*offset]+realInc;
        TEHZ[2*offset+1] = TE[2*offset+1] +imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==Vacuum)//vacuum zone
    {
      VacuumZone *pzonedata = dynamic_cast< VacuumZone *>(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        //total field=scattering field + incidence field
        PetscScalar realInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*cos(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                            -gnode[offset].y*sin(OptTheta))-PhiDiff);
        PetscScalar imageInc=sqrt(WeightOfTE*power*sqrt(EpsInVac/MuInVac))*sin(-WaveVector0*(gnode[offset].x*cos(OptTheta)
                             -gnode[offset].y*sin(OptTheta))-PhiDiff);
        pzonedata->aux[localIndex].OpHz = sqrt(pow((TE[2*offset]+realInc),2)+pow((TE[2*offset+1]+imageInc),2));
        TEHZ[2*offset] = TE[2*offset]+realInc;
        TEHZ[2*offset+1] = TE[2*offset+1] +imageInc;
      }
    }
    if(zonedata[ZoneSerial]->material_type==PML)//PML zone
    {
      PMLZone *pzonedata = dynamic_cast< PMLZone *>(zonedata[ZoneSerial]);
      for(int j=0;j<3;j++)
      {
        offset=gtri[i].g_node[j];
        int localIndex=gtri[i].node[j];
        pzonedata->aux[localIndex].OpHz = sqrt(pow(TE[2*offset],2)+pow(TE[2*offset+1],2));
        TEHZ[2*offset] = TE[2*offset];
        TEHZ[2*offset+1] = TE[2*offset+1];
      }
    }
  }
  return 0;
}


/*----------------------------------------------------------------------------------------------
*EM_FEM_Solver::generation_cal: this function calculates carrier generation rates
*/
int EM_FEM_Solver::generation_cal(SolveDefine &sv)
{ 
  for(int i=0;i<gtri.size();i++) //calculating generation rates
  {
    int ZoneSerial = gtri[i].zone_index;
    if(zonedata[ZoneSerial]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *>(zonedata[ZoneSerial]);
      PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
      PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
      RealPermi = n*n-k*k;
      ImagePermi = 2*n*k;
      PetscScalar PrivateEps0 = pzonedata->mt->eps0;
      PetscScalar PrivateMu0  = pzonedata->mt->mu0;
      PetscScalar PlankConst = pzonedata->mt->h;
      PetscScalar PrivateOmega = 2*PI/(lamda*sqrt(PrivateEps0*PrivateMu0));
      PetscScalar CoeffOfG  = eta*lamda*sqrt(PrivateEps0*PrivateMu0)/PlankConst;      
      
      for(int j=0;j<3;j++)
      {
        int LocalIndex = gtri[i].node[j];       
        int offset     = gtri[i].g_node[j];
        
        PetscScalar AreaJ = pzonedata->pzone->davcell.GetPointer(LocalIndex)->area;
        pzonedata->aux[LocalIndex].OptG += CoeffOfG*(PrivateOmega*PrivateEps0*ImagePermi*
                                          (pow(TMEZ[2*offset],2)+pow(TMEZ[2*offset+1],2)
                                          +pow(TEEX[2*offset],2)+pow(TEEX[2*offset+1],2)
                                          +pow(TEEY[2*offset],2)+pow(TEEY[2*offset+1],2))/2)*gtri[i].s[j]/AreaJ;                       
        } 
    }
  }
  return 0;
}

/*-------------------------------------------------------------------------------------------------
 * EM_FEM_Solver::fiel_calculation:  This function calculates other field components in both modes
 */
int EM_FEM_Solver::field_calculation(SolveDefine &sv)
{
  for(int i=0;i<gtri.size();i++)  //calculating by triangle element
  {
    int ZoneSerial=gtri[i].zone_index;
    if(zonedata[ZoneSerial]->material_type==Semiconductor) //semiconductor zone
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *>(zonedata[ZoneSerial]);
      PetscScalar n = pzonedata->mt->optical->RefractionIndex(lamda).real();
      PetscScalar k = pzonedata->mt->optical->RefractionIndex(lamda).imag();
      RealPermi = n*n-k*k;
      ImagePermi = 2*n*k;
      PetscScalar PrivateEps0 = pzonedata->mt->eps0;
      PetscScalar PrivateMu0  = pzonedata->mt->mu0;
      PetscScalar PrivateOmega = 2*PI/(lamda*sqrt(PrivateEps0*PrivateMu0));

      eleB[0] = gnode[gtri[i].g_node[1]].y-gnode[gtri[i].g_node[2]].y;
      eleB[1] = gnode[gtri[i].g_node[2]].y-gnode[gtri[i].g_node[0]].y;
      eleB[2] = gnode[gtri[i].g_node[0]].y-gnode[gtri[i].g_node[1]].y;
      eleC[0] = gnode[gtri[i].g_node[2]].x-gnode[gtri[i].g_node[1]].x;
      eleC[1] = gnode[gtri[i].g_node[0]].x-gnode[gtri[i].g_node[2]].x;
      eleC[2] = gnode[gtri[i].g_node[1]].x-gnode[gtri[i].g_node[0]].x;

      for(int j=0;j<3;j++)
      {
        int LocalIndex = gtri[i].node[j];
        int goffset    = gtri[i].g_node[j];

        PetscScalar RealHx  = 0.0;
        PetscScalar ImageHx = 0.0;
        PetscScalar RealHy  = 0.0;
        PetscScalar ImageHy = 0.0;
        PetscScalar RealEx  = 0.0;
        PetscScalar ImageEx = 0.0;
        PetscScalar RealEy  = 0.0;
        PetscScalar ImageEy = 0.0;

        if(WeightOfTM) //for TM mode
        {
          for(int k=0;k<3;k++)
          {
            int offset = gtri[i].g_node[k];
            RealHx = 1/(PrivateOmega*PrivateMu0)*eleC[k]*TMEZ[2*offset+1]/(2*gtri[i].area)+RealHx;
            ImageHx = -1/(PrivateOmega*PrivateMu0)*eleC[k]*TMEZ[2*offset]/(2*gtri[i].area)+ImageHx;
            RealHy = -1/(PrivateOmega*PrivateMu0)*eleB[k]*TMEZ[2*offset+1]/(2*gtri[i].area)+RealHy;
            ImageHy = 1/(PrivateOmega*PrivateMu0)*eleB[k]*TMEZ[2*offset]/(2*gtri[i].area)+ImageHy;
          }
          PetscScalar AreaJ = pzonedata->pzone->davcell.GetPointer(LocalIndex)->area;

          pzonedata->aux[LocalIndex].OpHx = pzonedata->aux[LocalIndex].OpHx+sqrt(pow(RealHx,2)+pow(ImageHx,2))*gtri[i].s[j]/AreaJ;
          pzonedata->aux[LocalIndex].OpHy = pzonedata->aux[LocalIndex].OpHy+sqrt(pow(RealHy,2)+pow(ImageHy,2))*gtri[i].s[j]/AreaJ;

          TMHX[2*goffset] = TMHX[2*goffset]+RealHx*gtri[i].s[j]/AreaJ;
          TMHX[2*goffset+1] = TMHX[2*goffset+1]+ImageHx*gtri[i].s[j]/AreaJ;
          TMHY[2*goffset] = TMHY[2*goffset]+RealHy*gtri[i].s[j]/AreaJ;
          TMHY[2*goffset+1] = TMHY[2*goffset+1]+ImageHy*gtri[i].s[j]/AreaJ;
        }

        if(WeightOfTE)  //for TE mode
        {
          for(int k=0;k<3;k++)
          {
            int offset = gtri[i].g_node[k];
            RealEx = -1/(PrivateOmega*PrivateEps0*(pow(RealPermi,2)+pow(ImagePermi,2)))*(RealPermi*eleC[k]*TEHZ[2*offset+1]
                     +ImagePermi*eleC[k]*TEHZ[2*offset])/(2*gtri[i].area)+RealEx;
            ImageEx = 1/(PrivateOmega*PrivateEps0*(pow(RealPermi,2)+pow(ImagePermi,2)))*(RealPermi*eleC[k]*TEHZ[2*offset]
                      -ImagePermi*eleC[k]*TEHZ[2*offset+1])/(2*gtri[i].area)+ImageEx;
            RealEy = 1/(PrivateOmega*PrivateEps0*(pow(RealPermi,2)+pow(ImagePermi,2)))*(RealPermi*eleB[k]*TEHZ[2*offset+1]
                     +ImagePermi*eleB[k]*TEHZ[2*offset])/(2*gtri[i].area)+RealEy;
            ImageEy = -1/(PrivateOmega*PrivateEps0*(pow(RealPermi,2)+pow(ImagePermi,2)))*(RealPermi*eleB[k]*TEHZ[2*offset]
                      -ImagePermi*eleB[k]*TEHZ[2*offset+1])/(2*gtri[i].area)+ImageEx;
          }
          PetscScalar AreaJ = pzonedata->pzone->davcell.GetPointer(LocalIndex)->area;

          pzonedata->aux[LocalIndex].OpEx = pzonedata->aux[LocalIndex].OpEx+sqrt(pow(RealEx,2)+pow(ImageEx,2))*gtri[i].s[j]/AreaJ;
          pzonedata->aux[LocalIndex].OpEy = pzonedata->aux[LocalIndex].OpEy+sqrt(pow(RealEy,2)+pow(ImageEy,2))*gtri[i].s[j]/AreaJ;

          TEEX[2*goffset] = TEEX[2*goffset]+RealEx*gtri[i].s[j]/AreaJ;
          TEEX[2*goffset+1] = TEEX[2*goffset+1]+ImageEx*gtri[i].s[j]/AreaJ;
          TEEY[2*goffset] = TEEY[2*goffset]+RealEy*gtri[i].s[j]/AreaJ;
          TEEY[2*goffset+1] = TEEY[2*goffset+1]+ImageEy*gtri[i].s[j]/AreaJ;
        }
      }
    }
  }
  return 0;
}


/* ----------------------------------------------------------------------------
 * EM_FEM_Solver::destroy_solver:  This function do destroy the FEM solver
 */
int EM_FEM_Solver::destroy_solver(SolveDefine &sv)
{
  ierr=VecDestroy(x);
  ierr=VecDestroy(b);
  ierr=MatDestroy(A);
  ierr=KSPDestroy(solver);
  return 0;
}
