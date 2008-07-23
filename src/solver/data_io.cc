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
/*  Last update: Nov 28, 2005                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#include "bsolver.h"
#include "zonedata.h"
#include "log.h"
#include <math.h>
#include <assert.h>
#include <cgnslib.h>

//simple error handle function
inline int cgns_error_handle(const char *s,int fn)
{
  gss_log.string_buf()<<s;
  gss_log.record();
  cg_close(fn);
  return 1;
}



int SMCZone::import_doping(char *cgnsfile,PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,BS=1,nfields;
  int  isize[3],nsol;
  char zonelabel[32],solname[32];
  GridLocation_t location;

  int node_num = pzone->davcell.size();

  double *Na, *Nd, *P, *As, *Sb, *B;
  Na = new double[node_num];
  Nd = new double[node_num];
  P  = new double[node_num];
  As = new double[node_num];
  Sb = new double[node_num];
  B  = new double[node_num];

  cg_open(cgnsfile,MODE_READ,&fn);
  cg_zone_read(fn,BS,Z,zonelabel,isize);

  gss_log.string_buf()<<"read Doping of Region "<<pzone->zonename<<"... ";
  //---------------------------------------------------------------------------------------
  //how many fieldsol node?
  assert(!cg_nsols(fn,BS,Z,&nsol));
  if(nsol == 0)
    return cgns_error_handle("no fieldsol node,not a supported cgns file\n",fn);
  //search for the "doping" FlowSolution_t
  int flag = 0;
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,BS,Z,i,solname,&location));
    if(!strcmp(solname,"Doping"))
    {
      flag = 1;//we found doping information
      int imin=1;

      cg_nfields(fn, BS, Z, i, &nfields);
      if(nfields <= 2)
      {
        assert(!cg_field_read(fn,BS,Z,i,"Nd",RealDouble,&imin,&isize[0],Nd));
        assert(!cg_field_read(fn,BS,Z,i,"Na",RealDouble,&imin,&isize[0],Na));
      }
      if(nfields == 6)
      {
        assert(!cg_field_read(fn,BS,Z,i,"Nd",RealDouble,&imin,&isize[0],Nd));
        assert(!cg_field_read(fn,BS,Z,i,"Na",RealDouble,&imin,&isize[0],Na));
	assert(!cg_field_read(fn,BS,Z,i,"P" ,RealDouble,&imin,&isize[0],P ));
        assert(!cg_field_read(fn,BS,Z,i,"As",RealDouble,&imin,&isize[0],As));
        assert(!cg_field_read(fn,BS,Z,i,"Sb",RealDouble,&imin,&isize[0],Sb));
        assert(!cg_field_read(fn,BS,Z,i,"B" ,RealDouble,&imin,&isize[0],B ));
      }
      gss_log.string_buf()<<"done\n";
      gss_log.record();
    }
  }

  if(flag == 0 )
    return cgns_error_handle("no doping found",fn);
  else
    cg_close(fn);

  if(nfields <= 2)
  {
    for(int i=0;i<node_num;i++)
    {     //also,doping will be scaled. nagtive value convert to positive
      aux[i].Nd = fabs(Nd[i])/pow(scale->s_centimeter,3);
      aux[i].Na = fabs(Na[i])/pow(scale->s_centimeter,3);
    }
  }
  if(nfields == 6)
  {
    for(int i=0;i<node_num;i++)
    {     //also,doping will be scaled. nagtive value convert to positive
      aux[i].Nd = fabs(Nd[i])/pow(scale->s_centimeter,3);
      aux[i].Na = fabs(Na[i])/pow(scale->s_centimeter,3);
      aux[i].P  = fabs(P[i])/pow(scale->s_centimeter,3);
      aux[i].As = fabs(As[i])/pow(scale->s_centimeter,3);
      aux[i].Sb = fabs(Sb[i])/pow(scale->s_centimeter,3);
      aux[i].B  = fabs(B[i])/pow(scale->s_centimeter,3);
    }
  }
  
  delete [] Na;
  delete [] Nd;
  delete [] P ;
  delete [] As;
  delete [] Sb;
  delete [] B ;
  return 0;

}


int SMCZone::export_doping(char *cgnsfile,PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,BS=1,SOL,F;
  int  isize[3],nsol;
  char zonelabel[32],solname[32];
  GridLocation_t location;

  int node_num = pzone->davcell.size();

  double *Na, *Nd, *P, *As, *Sb, *B;
  Na = new double[node_num];
  Nd = new double[node_num];
  P  = new double[node_num];
  As = new double[node_num];
  Sb = new double[node_num];
  B  = new double[node_num];

  for(int i=0;i<node_num;i++)
  {
    Nd[i] = double(aux[i].Nd*pow(scale->s_centimeter,3));
    Na[i] = double(aux[i].Na*pow(scale->s_centimeter,3));
    P[i]  = double(aux[i].P *pow(scale->s_centimeter,3));
    As[i] = double(aux[i].As*pow(scale->s_centimeter,3));
    Sb[i] = double(aux[i].Sb*pow(scale->s_centimeter,3));
    B[i]  = double(aux[i].B *pow(scale->s_centimeter,3));
  }

  cg_open(cgnsfile,MODE_MODIFY,&fn);
  cg_zone_read(fn,BS,Z,zonelabel,isize);

  //---------------------------------------------------------------------------------------
  //search for FlowSolution_t if the "Doping" is exit
  assert(!cg_nsols(fn,BS,Z,&nsol));
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,BS,Z,i,solname,&location));
    if(!strcmp(solname,"Doping"))
    {
      gss_log.string_buf()<<"waring !!! doping exit,delete it.\n";
      gss_log.record();
      assert(!cg_goto(fn,BS,"Zone_t",Z,"end"));
      assert(!cg_delete_node("Doping"));
    }
  }

  //write field solution

  assert(!cg_sol_write(fn,BS,Z,"Doping",Vertex,&SOL));
  assert(!cg_field_write(fn,BS,Z,SOL,RealDouble,"Na",Na,&F));
  assert(!cg_field_write(fn,BS,Z,SOL,RealDouble,"Nd",Nd,&F));
  assert(!cg_field_write(fn,BS,Z,SOL,RealDouble,"P", P, &F));
  assert(!cg_field_write(fn,BS,Z,SOL,RealDouble,"As",As,&F));
  assert(!cg_field_write(fn,BS,Z,SOL,RealDouble,"Sb",Sb,&F));
  assert(!cg_field_write(fn,BS,Z,SOL,RealDouble,"B", B, &F));

  cg_close(fn);

  delete [] Na;
  delete [] Nd;
  delete [] P ;
  delete [] As;
  delete [] Sb;
  delete [] B ;
  return 0;

}


int SMCZone::export_accptor_ascii(char *dopingfile,PhysicalUnitScale *scale)
{
  FILE *fp=fopen(dopingfile,"a");//opne in append model, for multi-region
  if (fp == NULL) return 1;
  for(int i=0;i<node_num;i++)
  {
    double x = pzone->danode[i].x/scale->s_micron;
    double y = pzone->danode[i].y/scale->s_micron;
    fprintf(fp,"%e %e %e\n",x,y,double((aux[i].Na+ aux[i].B) *pow(scale->s_centimeter,3)));
  }
  fclose(fp);
  return 0;
}


int SMCZone::export_donor_ascii(char *dopingfile,PhysicalUnitScale *scale)
{
  FILE *fp=fopen(dopingfile,"a");//opne in append model, for multi-region
  if (fp == NULL) return 1;
  for(int i=0;i<node_num;i++)
  {
    double x = pzone->danode[i].x/scale->s_micron;
    double y = pzone->danode[i].y/scale->s_micron;
    fprintf(fp,"%e %e %e\n",x,y,double((aux[i].Nd+aux[i].P+aux[i].As+aux[i].Sb)*pow(scale->s_centimeter,3)));
  }
  fclose(fp);
  return 0;
}


int SMCZone::import_mole(char *cgnsfile,int comp)
{
  int  Z = zone_index+1,fn,B=1;
  int  isize[3],nsol;
  char zonelabel[32],solname[32];
  GridLocation_t location;

  int node_num = pzone->davcell.size();

  double *mole_x,*mole_y;
  mole_x = new double[node_num];
  if(comp==2) mole_y = new double[node_num];

  cg_open(cgnsfile,MODE_READ,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);

  gss_log.string_buf()<<"read mole function of Region "<<pzone->zonename<<"... ";

  //---------------------------------------------------------------------------------------
  //how many fieldsol node?
  assert(!cg_nsols(fn,B,Z,&nsol));
  if(nsol == 0)
    return cgns_error_handle("no fieldsol node,not a supported cgns file\n",fn);
  //search for the "doping" FlowSolution_t
  int flag = 0;
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,solname,&location));
    if(!strcmp(solname,"Mole"))
    {
      int imin = 1;
      flag = 1;
      gss_log.string_buf()<<"done\n";
      gss_log.record();
      assert(!cg_field_read(fn,B,Z,i,"mole_x",RealDouble,&imin,&isize[0],mole_x));
      if(comp==2) assert(!cg_field_read(fn,B,Z,i,"mole_y",RealDouble,&imin,&isize[0],mole_y));
    }
  }

  if(flag == 0 )
    return cgns_error_handle("no mole data found",fn);
  else
    cg_close(fn);

  for(int i=0;i<node_num;i++)
  {     //also,doping will be scaled. nagtive value convert to positive
    aux[i].mole_x = mole_x[i];
    if(comp==2) aux[i].mole_y = mole_y[i];
  }

  delete [] mole_x;
  if(comp==2) delete [] mole_y;
  return 0;

}


int SMCZone::export_mole(char *cgnsfile,int comp)
{
  int  Z = zone_index+1,fn,B=1,SOL,F;
  int  isize[3],nsol;
  char zonelabel[32],solname[32];
  GridLocation_t location;

  int node_num = pzone->davcell.size();

  double *mole_x,*mole_y;
  mole_x = new double[node_num];
  if(comp==2) mole_y = new double[node_num];

  for(int i=0;i<node_num;i++)
  {
    mole_x[i] = double(aux[i].mole_x);
    if(comp==2) mole_y[i] = double(aux[i].mole_y);
  }

  cg_open(cgnsfile,MODE_MODIFY,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);

  //---------------------------------------------------------------------------------------
  //search for FlowSolution_t if the "Doping" is exit
  assert(!cg_nsols(fn,B,Z,&nsol));
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,solname,&location));
    if(!strcmp(solname,"Mole"))
    {
      gss_log.string_buf()<<"waring !!! mole data exit,delete it.\n";
      gss_log.record();
      assert(!cg_goto(fn,B,"Zone_t",Z,"end"));
      assert(!cg_delete_node("Mole"));
    }
  }

  //write field solution

  assert(!cg_sol_write(fn,B,Z,"Mole",Vertex,&SOL));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"mole_x",mole_x,&F));
  if(comp==2) assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"mole_y",mole_y,&F));

  cg_close(fn);

  delete [] mole_x;
  if(comp==2) delete [] mole_y;
  return 0;

}


int SMCZone::export_solution(char *cgnsfile,char *solname, DABC &bc, PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,B=1,SOL,F,BCNum;
  char zonelabel[32],sol_exist[32],boconame[32];
  int  isize[3],nsol,npnts,normalindex,normallistflag,ndataset;
  GridLocation_t location;
  DataType_t     normaldatatype;
  PointSetType_t ptset;
  BCType_t       bocotype;

  int node_num = pzone->davcell.size();

  double *n,*p,*P,*T,*Tn,*Tp,*Eqc,*Eqv;

  n  = new double[node_num];
  p  = new double[node_num];
  P  = new double[node_num];
  T  = new double[node_num];
  Tn  = new double[node_num];
  Tp  = new double[node_num];
  Eqc = new double[node_num];
  Eqv = new double[node_num];
  for(int i=0;i<node_num;i++)
  {
    n[i]   = double(fs[i].n*pow(scale->s_centimeter,3));        //density unit cm-3
    p[i]   = double(fs[i].p*pow(scale->s_centimeter,3));
    P[i]   = double(fs[i].P/scale->s_volt);
    T[i]   = double(fs[i].T/scale->s_kelvin);
    Tn[i]  = double(fs[i].Tn/scale->s_kelvin);
    Tp[i]  = double(fs[i].Tp/scale->s_kelvin);
    Eqc[i] = double(fs[i].Eqc/scale->s_eV);
    Eqv[i] = double(fs[i].Eqv/scale->s_eV);
  }

  cg_open(cgnsfile,MODE_MODIFY,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);

  //------------------------------------------------------------
  //search for FlowSolution_t if the "solution" is exit
  assert(!cg_nsols(fn,B,Z,&nsol));
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,sol_exist,&location));
    if(!strcmp(sol_exist,solname))
    {
      gss_log.string_buf()<<"waring !!! Solution exit,delete it.\n";
      gss_log.record();
      assert(!cg_goto(fn,B,"Zone_t",Z,"end"));
      assert(!cg_delete_node(solname));
    }
  }
  assert(!cg_sol_write(fn,B,Z,solname,Vertex,&SOL));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"elec_density",n,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"hole_density",p,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"potential",P,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"lattice_temperature",T,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"elec_temperature",Tn,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"hole_temperature",Tp,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"elec_quantum_potential",Eqc,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"hole_quantum_potential",Eqv,&F));

  //-------------------------------------------------------------
  //write electrode potential data
  assert(!cg_nbocos(fn,B,Z,&BCNum));
  for(int BCIndex=1;BCIndex<=BCNum;BCIndex++)
  {
    assert(!cg_boco_info(fn,B,Z,BCIndex,boconame,&bocotype,&ptset,&npnts,
                         &normalindex,&normallistflag,&normaldatatype,&ndataset));
    int bc_index=bc.Get_bc_index(zone_index,boconame);
    if(bc.is_electrode(bc_index))
    {
      double potential=bc.Get_pointer(bc_index)->Get_Potential();
      int    DimensionVector=1;
      cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BCIndex, "end");
      cg_user_data_write ("Extra_data_for_electrode");
      cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BCIndex, "UserDefinedData_t", 1, "end");
      cg_array_write("electrode_potential", RealDouble, 1, &DimensionVector, &potential);
    }
  }

  cg_close(fn);

  //gss_log.string_buf()<<"Export Data for Semiconductor Zone.\n";
  //gss_log.record();

  delete [] n;
  delete [] p;
  delete [] P;
  delete [] T;
  delete [] Tn;
  delete [] Tp;
  delete [] Eqc;
  delete [] Eqv;
  return 0;

}

int ElZone::export_solution(char *cgnsfile,char *solname,DABC &bc,PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,B=1,SOL,F;
  int  isize[3],nsol;
  char zonelabel[32],sol_exist[32];
  GridLocation_t location;

  int node_num = pzone->davcell.size();

  double *P,*T;

  P  = new double[node_num];
  T  = new double[node_num];

  for(int i=0;i<node_num;i++)
  {
    P[i]   = double(fs[i].P/scale->s_volt);
    T[i]   = double(fs[i].T/scale->s_kelvin);
  }

  cg_open(cgnsfile,MODE_MODIFY,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);
  //------------------------------------------------------------
  //search for FlowSolution_t if the "solution" is exit
  assert(!cg_nsols(fn,B,Z,&nsol));
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,sol_exist,&location));
    if(!strcmp(sol_exist,solname))
    {
      gss_log.string_buf()<<"waring !!! Solution exit,delete it.\n";
      gss_log.record();
      assert(!cg_goto(fn,B,"Zone_t",Z,"end"));
      assert(!cg_delete_node(solname));
    }
  }
  //-----------------------------------------------------------
  assert(!cg_sol_write(fn,B,Z,solname,Vertex,&SOL));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"potential",P,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"lattice_temperature",T,&F));

  cg_close(fn);

  delete [] P;
  delete [] T;
  return 0;

}

int ISZone::export_solution(char *cgnsfile,char *solname,DABC &bc,PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,B=1,SOL,F,BCNum;
  char zonelabel[32],sol_exist[32],boconame[32];
  int  isize[3],nsol,npnts,normalindex,normallistflag,ndataset;
  GridLocation_t location;
  DataType_t     normaldatatype;
  PointSetType_t ptset;
  BCType_t       bocotype;

  int node_num = pzone->davcell.size();
  double *P,*T;

  P  = new double[node_num];
  T  = new double[node_num];

  for(int i=0;i<node_num;i++)
  {
    P[i]   = double(fs[i].P/scale->s_volt);
    T[i]   = double(fs[i].T/scale->s_kelvin);
  }

  cg_open(cgnsfile,MODE_MODIFY,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);
  //------------------------------------------------------------
  //search for FlowSolution_t if the "solution" is exit
  assert(!cg_nsols(fn,B,Z,&nsol));
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,sol_exist,&location));
    if(!strcmp(sol_exist,solname))
    {
      gss_log.string_buf()<<"waring !!! Solution exit,delete it.\n";
      gss_log.record();
      assert(!cg_goto(fn,B,"Zone_t",Z,"end"));
      assert(!cg_delete_node(solname));
    }
  }
  assert(!cg_sol_write(fn,B,Z,solname,Vertex,&SOL));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"potential",P,&F));
  assert(!cg_field_write(fn,B,Z,SOL,RealDouble,"lattice_temperature",T,&F));


  //-------------------------------------------------------------
  //write electrode potential data
  assert(!cg_nbocos(fn,B,Z,&BCNum));
  for(int BCIndex=1;BCIndex<=BCNum;BCIndex++)
  {
    assert(!cg_boco_info(fn,B,Z,BCIndex,boconame,&bocotype,&ptset,&npnts,
                         &normalindex,&normallistflag,&normaldatatype,&ndataset));
    int bc_index=bc.Get_bc_index(zone_index,boconame);
    if(bc.is_electrode(bc_index))
    {
      double potential=bc.Get_pointer(bc_index)->Get_Potential();
      int    DimensionVector=1;
      cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BCIndex, "end");
      cg_user_data_write ("Extra_data_for_electrode");
      cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BCIndex, "UserDefinedData_t", 1, "end");
      cg_array_write("electrode_potential", RealDouble, 1, &DimensionVector, &potential);
    }
  }

  cg_close(fn);

  delete [] P;
  delete [] T;
  return 0;

}



int SMCZone::import_solution(char *cgnsfile,char *solname,DABC &bc,PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,B=1,BCNum;
  int  isize[3],nsol,npnts,normalindex,normallistflag,ndataset;
  char zonelabel[32],sol[32],boconame[32];
  GridLocation_t location;
  DataType_t     normaldatatype;
  PointSetType_t ptset;
  BCType_t       bocotype;
  int narrays=0;

  char ArrayName[32]="";
  DataType_t DataType=RealDouble;
  int  DataDimension=1,DimensionVector=1;
  int nfields;
  int node_num = pzone->davcell.size();

  double *n,*p,*P,*T,*Tn,*Tp,*Eqc,*Eqv;

  n  = new double[node_num];
  p  = new double[node_num];
  P  = new double[node_num];
  T  = new double[node_num];
  Tn = new double[node_num];
  Tp = new double[node_num];
  Eqc = new double[node_num];
  Eqv = new double[node_num];

  cg_open(cgnsfile,MODE_READ,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);

  gss_log.string_buf()<<"read Solution of Region "<<pzone->zonename<<"... ";
  //---------------------------------------------------------------------------------------
  //how many fieldsol node?
  assert(!cg_nsols(fn,B,Z,&nsol));
  if(nsol == 0)
    return cgns_error_handle("no fieldsol node,not a supported CoreFile\n",fn);
  //search for the "solution" FlowSolution_t
  int flag = 0;
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,sol,&location));
    if(!strcmp(sol,solname))
    {
      int imin = 1;
      flag = 1;
      gss_log.string_buf()<<"done\n";
      gss_log.record();
      cg_nfields(fn,B,Z,i,&nfields);

      assert(!cg_field_read(fn,B,Z,i,"potential",RealDouble,&imin,&isize[0],P));
      assert(!cg_field_read(fn,B,Z,i,"elec_density",RealDouble,&imin,&isize[0],n));
      assert(!cg_field_read(fn,B,Z,i,"hole_density",RealDouble,&imin,&isize[0],p));
      assert(!cg_field_read(fn,B,Z,i,"lattice_temperature",RealDouble,&imin,&isize[0],T));
      if(nfields>4)
      {
        assert(!cg_field_read(fn,B,Z,i,"elec_temperature",RealDouble,&imin,&isize[0],Tn));
        assert(!cg_field_read(fn,B,Z,i,"hole_temperature",RealDouble,&imin,&isize[0],Tp));
      }
      if(nfields>6)
      {
        assert(!cg_field_read(fn,B,Z,i,"elec_quantum_potential",RealDouble,&imin,&isize[0],Eqc));
        assert(!cg_field_read(fn,B,Z,i,"hole_quantum_potential",RealDouble,&imin,&isize[0],Eqv));
      }
    }
  }
  for(int i=0;i<node_num;i++)
  {
    fs[i].n = n[i]/pow(scale->s_centimeter,3);
    fs[i].p = p[i]/pow(scale->s_centimeter,3);
    fs[i].P = P[i]*scale->s_volt;
    fs[i].T = T[i]*scale->s_kelvin;
    if(nfields>4)
    {
      fs[i].Tn = Tn[i]*scale->s_kelvin;
      fs[i].Tp = Tp[i]*scale->s_kelvin;
    }
    else
    {
      fs[i].Tn = fs[i].T;
      fs[i].Tp = fs[i].T;
    }
    if(nfields>6)
    {
      fs[i].Eqc = Eqc[i]*scale->s_eV;
      fs[i].Eqv = Eqv[i]*scale->s_eV;
    }
    else
    {
      fs[i].Eqc = -scale->s_eV*(fs[i].P + aux[i].affinity);
      fs[i].Eqv = -scale->s_eV*(fs[i].P + aux[i].affinity + aux[i].Eg);
    }
  }

  delete [] n;
  delete [] p;
  delete [] P;
  delete [] T;
  delete [] Tn;
  delete [] Tp;
  delete [] Eqc;
  delete [] Eqv;
  //read bc potential
  assert(!cg_nbocos(fn,B,Z,&BCNum));
  for(int BCIndex=1;BCIndex<=BCNum;BCIndex++)
  {
    assert(!cg_boco_info(fn,B,Z,BCIndex,boconame,&bocotype,&ptset,&npnts,
                         &normalindex,&normallistflag,&normaldatatype,&ndataset));
    int bc_index=bc.Get_bc_index(zone_index,boconame);
    if(bc.is_electrode(bc_index))
    {
      double potential=0;
      cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BCIndex, "UserDefinedData_t", 1, "end");
      assert(!cg_array_read_as(1,RealDouble, &potential));
      bc.Get_pointer(bc_index)->Set_Potential(potential);
    }
  }

  if(flag == 0 )
    return cgns_error_handle("no Core Data found\n",fn);
  else
    cg_close(fn);

  return 0;
}


int ISZone::import_solution(char *cgnsfile,char *solname,DABC &bc,PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,B=1,BCNum;
  int  isize[3],nsol,npnts,normalindex,normallistflag,ndataset;
  char zonelabel[32],sol[32],boconame[32];
  GridLocation_t location;
  DataType_t     normaldatatype;
  PointSetType_t ptset;
  BCType_t       bocotype;
  int node_num = pzone->davcell.size();

  double *n,*p,*P,*T;

  P  = new double[node_num];
  T  = new double[node_num];


  cg_open(cgnsfile,MODE_READ,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);

  gss_log.string_buf()<<"read Solution of Region "<<pzone->zonename<<"... ";
  //---------------------------------------------------------------------------------------
  //how many fieldsol node?
  assert(!cg_nsols(fn,B,Z,&nsol));
  if(nsol == 0)
    return cgns_error_handle("no fieldsol node,not a supported CoreFile\n",fn);
  //search for the "solution" FlowSolution_t
  int flag = 0;
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,sol,&location));
    if(!strcmp(sol,solname))
    {
      int imin = 1;
      flag = 1;
      gss_log.string_buf()<<"done\n";
      gss_log.record();
      assert(!cg_field_read(fn,B,Z,i,"lattice_temperature",RealDouble,&imin,&isize[0],T));
      assert(!cg_field_read(fn,B,Z,i,"potential",RealDouble,&imin,&isize[0],P));
    }
  }
  for(int i=0;i<node_num;i++)
  {
    fs[i].P = P[i]*scale->s_volt;
    fs[i].T = T[i]*scale->s_kelvin;
  }

  delete [] P;
  delete [] T;

  //read bc potential
  assert(!cg_nbocos(fn,B,Z,&BCNum));
  for(int BCIndex=1;BCIndex<=BCNum;BCIndex++)
  {
    assert(!cg_boco_info(fn,B,Z,BCIndex,boconame,&bocotype,&ptset,&npnts,
                         &normalindex,&normallistflag,&normaldatatype,&ndataset));
    int bc_index=bc.Get_bc_index(zone_index,boconame);
    if(bc.is_electrode(bc_index))
    {
      double potential=0;
      cg_goto(fn, B, "Zone_t", Z, "ZoneBC_t", 1, "BC_t", BCIndex, "UserDefinedData_t", 1, "end");
      assert(!cg_array_read_as(1,RealDouble, &potential));
      bc.Get_pointer(bc_index)->Set_Potential(potential);
    }
  }

  if(flag == 0 )
    return cgns_error_handle("no Core Data found\n",fn);
  else
    cg_close(fn);

  return 0;
}

int ElZone::import_solution(char *cgnsfile,char *solname,DABC &bc,PhysicalUnitScale *scale)
{
  int  Z = zone_index+1,fn,B=1;
  int  isize[3],nsol;
  char zonelabel[32],sol[32];
  GridLocation_t location;
  int node_num = pzone->davcell.size();

  double *n,*p,*P,*T;

  P  = new double[node_num];
  T  = new double[node_num];


  cg_open(cgnsfile,MODE_READ,&fn);
  cg_zone_read(fn,B,Z,zonelabel,isize);

  gss_log.string_buf()<<"read Solution of Region "<<pzone->zonename<<"... ";
  //---------------------------------------------------------------------------------------
  //how many fieldsol node?
  assert(!cg_nsols(fn,B,Z,&nsol));
  if(nsol == 0)
    return cgns_error_handle("no fieldsol node,not a supported CoreFile\n",fn);
  //search for the "solution" FlowSolution_t
  int flag = 0;
  for(int i=1;i<=nsol;i++)
  {
    assert(!cg_sol_info(fn,B,Z,i,sol,&location));
    if(!strcmp(sol,solname))
    {
      int imin = 1;
      flag = 1;
      gss_log.string_buf()<<"done\n";
      gss_log.record();
      assert(!cg_field_read(fn,B,Z,i,"lattice_temperature",RealDouble,&imin,&isize[0],T));
      assert(!cg_field_read(fn,B,Z,i,"potential",RealDouble,&imin,&isize[0],P));
    }
  }
  for(int i=0;i<node_num;i++)
  {
    fs[i].P = P[i]*scale->s_volt;
    fs[i].T = T[i]*scale->s_kelvin;
  }

  delete [] P;
  delete [] T;

  if(flag == 0 )
    return cgns_error_handle("no Core Data found\n",fn);
  else
    cg_close(fn);

  return 0;
}

int BSolver::extract_ascii(char *filename)
{

  FILE *fp;
  if((fp=fopen(filename,"w"))==NULL) return 1;

  zone_to_field();

  fprintf(fp,"h TIF V1.2.1 generated by GSS 0.4x\n");

  //output nodes
  for(int i=0; i<gnode.size();i++)
  {
    fprintf(fp,"c\t%d\t%e\t%e\t%e\n",
            i+1,
            gnode[i].x/scale_unit.s_centimeter*1e4,
            -gnode[i].y/scale_unit.s_centimeter*1e4,
            0.0);
  }

  //output edges
  int edge_index = 1;
  for(int i=0; i<gsegment.size();i++)
    for(int j=0;j<gsegment[i].edge_array.size();j++,edge_index++)
    {
      fprintf(fp,"e\t%d\t%d\t%d\t%d\n",
              edge_index,
              gsegment[i].edge_array[j].p1+1,
              gsegment[i].edge_array[j].p2+1,
              0);
    }

  //output region, electrode and boundary type
  int elec_index = 1;
  for(int z=0;z<zone_num;z++)
  {
    // regions
    fprintf(fp,"r\t%d\t%s\t%s\n",
            z+1,
            zone[z].zonelabel,
            zone[z].zonename
           );
    // output boundarys
    for(int s=0; s<zone[z].dasegment.size();s++)
    {
      int g_index = zone[z].dasegment[s].g_index;
      edge_index = 1;
      for(int i=0; i<g_index;i++) edge_index+=gsegment[i].edge_num;
      for(int j=0;j<gsegment[g_index].edge_array.size()/2;j++,edge_index++)
        fprintf(fp,"\tb %d\n",edge_index);
    }
    //skip electrode region
    if(!strcmp(zone[z].zonelabel,"Elec"))
    {
      //fprintf(fp,"i %d %s %s %d\n",elec_index++,"Elec",zone[z].zonename,1);
    }
    //electrode boundary
    else
    {
      edge_index = 1;
      for(int i=0; i<gsegment.size();i++)
      {
        int bc_index = gsegment[i].bc_index - 1;
        if(gsegment[i].zone_index==z&&
            gsegment[i].interface==-1&&
            (bc[bc_index].BCType==OhmicContact||
             bc[bc_index].BCType==SchottkyContact))
        {
          fprintf(fp,"i %d %s %s %d\n",elec_index++,"Elec",gsegment[i].label,0);
          for(int j=0;j<gsegment[i].edge_array.size()/2;j++,edge_index++)
            fprintf(fp,"\tj %d\n",edge_index);
        }
        else
          edge_index+=gsegment[i].edge_num;
      }
    }
  }

  //output electrode region here
  for(int z=0;z<zone_num;z++)
  {
    if(!strcmp(zone[z].zonelabel,"Elec"))
    {
      fprintf(fp,"i %d %s %s %d\n",elec_index++,"Elec",zone[z].zonename,1);
    }
  }

  //output triangles
  int tri_index = 1;
  for(int i=0; i<gtri.size();i++)
  {
    fprintf(fp,"t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
            tri_index++,
            gtri[i].zone_index+1,
            gtri[i].g_node[0]+1,
            gtri[i].g_node[1]+1,
            gtri[i].g_node[2]+1,
            -1,
            -1,
            -1);
  }

  // set all the interface node in global node array as semiconductor node
  for(int z=0;z<zone_num;z++)
    if(zonedata[z]->material_type == Semiconductor)
    {
      for(int i=0;i<zone[z].danode.size();i++)
      {
        gnode[zone[z].danode[i].g_index].local_index = i;
        gnode[zone[z].danode[i].g_index].zone_index = z;
      }
    }
  fprintf(fp,"s\t%d Donor Accept n p v Ex Ey T\n",8);

  for(int i=0; i<gnode.size();i++)
  {
    int z = gnode[i].zone_index;
    int j = gnode[i].local_index;
    if(zonedata[z]->material_type == Semiconductor)
    {
      SMCZone *pzonedata = dynamic_cast<SMCZone *>(zonedata[z]);
      fprintf(fp,"n %d %s %e %e %e %e %e %e %e %e\n",
              i+1,
              zone[z].zonelabel,
              double(pzonedata->aux[j].Total_Nd()*pow(scale_unit.s_centimeter,3)),
              double(pzonedata->aux[j].Total_Na()*pow(scale_unit.s_centimeter,3)),
              double(pzonedata->fs[j].n*pow(scale_unit.s_centimeter,3)),
              double(pzonedata->fs[j].p*pow(scale_unit.s_centimeter,3)),
              double(pzonedata->fs[j].P/scale_unit.s_volt),
              double(pzonedata->aux[j].Ex/scale_unit.s_volt*scale_unit.s_centimeter),
              double(pzonedata->aux[j].Ey/scale_unit.s_volt*scale_unit.s_centimeter),
              double(pzonedata->fs[j].T/scale_unit.s_kelvin));
    }

    if(zonedata[z]->material_type == Insulator)
    {
      ISZone *pzonedata = dynamic_cast<ISZone *>(zonedata[z]);
      fprintf(fp,"n %d %s %e %e %e %e %e %e %e %e\n",
              i+1,
              zone[z].zonelabel,
              0.0,
              0.0,
              0.0,
              0.0,
              double(pzonedata->fs[j].P/scale_unit.s_volt),
              double(pzonedata->aux[j].Ex/scale_unit.s_volt*scale_unit.s_centimeter),
              double(pzonedata->aux[j].Ey/scale_unit.s_volt*scale_unit.s_centimeter),
              double(pzonedata->fs[j].T/scale_unit.s_kelvin));
    }

    if(zonedata[z]->material_type == Conductor)
    {
      ElZone *pzonedata = dynamic_cast<ElZone *>(zonedata[z]);
      fprintf(fp,"n %d %s %e %e %e %e %e %e %e %e\n",
              i+1,
              zone[z].zonelabel,
              0.0,
              0.0,
              0.0,
              0.0,
              double(pzonedata->fs[j].P/scale_unit.s_volt),
              double(pzonedata->aux[j].Ex/scale_unit.s_volt*scale_unit.s_centimeter),
              double(pzonedata->aux[j].Ey/scale_unit.s_volt*scale_unit.s_centimeter),
              double(pzonedata->fs[j].T/scale_unit.s_kelvin));
    }
  }

  fclose(fp);

  return 0;
}
