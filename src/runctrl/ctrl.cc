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
/*  Last update: Apr 13, 2006                                                */
/*                                                                           */
/*  Gong Ding                                                                */
/*  gdiso@ustc.edu                                                           */
/*  NINT, No.69 P.O.Box, Xi'an City, China                                   */
/*                                                                           */
/*****************************************************************************/
#include "config.h"
#include "log.h"
#include "ctrl.h"

#include "xgraph.h"
#include "wgraph.h"
#include "grafix3d.h"
#include <ctype.h>


int  SolveControl::run_control(list<Cmd>& cmdlist)
{

  //search the CmdBuf, delete each command card after its execution.
  for(pcmdbuf->cmd_search_begin();!pcmdbuf->cmd_search_end();)
  {
    if(pcmdbuf->is_current_cmd("METHOD"))
    {
      if(set_method(pcmdbuf->get_current_cmd())) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("SOLVE"))
    {
      if(set_solve(pcmdbuf->get_current_cmd())) return 1;
      if(do_solve()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("EXPORT"))
    {
      if(set_export(pcmdbuf->get_current_cmd())) return 1;
      if(do_export()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("IMPORT"))
    {
      if(set_import(pcmdbuf->get_current_cmd())) return 1;
      if(do_import()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("ATTACH"))
    {
      if(set_attach(pcmdbuf->get_current_cmd())) return 1;
      if(do_attach()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("REFINE"))
    {
      if(set_refine(pcmdbuf->get_current_cmd())) return 1;
      if(do_refine()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("PLOT"))
    {
      if(set_plot(pcmdbuf->get_current_cmd())) return 1;
      if(do_plot()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("PLOTMESH"))
    {
      if(set_plotmesh(pcmdbuf->get_current_cmd())) return 1;
      if(do_plotmesh()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("PLOTVTK"))
    {
      if(set_vtk_plot(pcmdbuf->get_current_cmd())) return 1;
      if(do_vtk_plot()) return 1;
      pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("PROBE"))
    {
      if(set_probe(pcmdbuf->get_current_cmd())) return 1;
       pcmdbuf->delete_current_cmd();
    }

    else if(pcmdbuf->is_current_cmd("END"))
    {
      return 0;
    }

    else
      pcmdbuf->goto_next_cmd();
  }

  return 0;

}



int SolveControl::do_solve()
{

  if(solve_define.Solver == DDML1E)
  {
    DDM_Solver_L1E::init_solver(solve_define);
    DDM_Solver_L1E::do_solve(solve_define);
    DDM_Solver_L1E::destroy_solver(solve_define);
  }

  if(solve_define.Solver == DDML1AC)
  {
    DDM_Solver_L1AC::init_solver(solve_define);
    DDM_Solver_L1AC::do_solve(solve_define);
    DDM_Solver_L1AC::destroy_solver(solve_define);
  }

  if(solve_define.Solver == QDDML1E)
  {
    QDDM_Solver_L1E::init_solver(solve_define);
    QDDM_Solver_L1E::do_solve(solve_define);
    QDDM_Solver_L1E::destroy_solver(solve_define);
  }

  if(solve_define.Solver == DDML1MIX)
  {
    DDM_Mix_Solver_L1E::init_solver(solve_define);
    DDM_Mix_Solver_L1E::do_solve(solve_define);
    DDM_Mix_Solver_L1E::destroy_solver(solve_define);
  }

  if(solve_define.Solver == DDML2E)
  {
    DDM_Solver_L2E::init_solver(solve_define);
    DDM_Solver_L2E::do_solve(solve_define);
    DDM_Solver_L2E::destroy_solver(solve_define);
  }

  if(solve_define.Solver == DDML2MIX)
  {
    DDM_Mix_Solver_L2E::init_solver(solve_define);
    DDM_Mix_Solver_L2E::do_solve(solve_define);
    DDM_Mix_Solver_L2E::destroy_solver(solve_define);
  }

  if(solve_define.Solver == EBML3E)
  {
    EBM_Solver_L3E::init_solver(solve_define);
    EBM_Solver_L3E::do_solve(solve_define);
    EBM_Solver_L3E::destroy_solver(solve_define);
  }

  if(solve_define.Solver == EMFEM)
  {
    if(EM_FEM_Solver::init_solver(solve_define))    return 1;
    if(EM_FEM_Solver::do_solve(solve_define))       return 1;
    if(EM_FEM_Solver::destroy_solver(solve_define)) return 1;
  }

  return 0;
}

int  SolveControl::do_refine()
{
  sprintf(log_buf,"Refine mesh...\n");GSS_LOG();
  //delete old solution data
  if(BSolver::refine(refine_define)) return 1;
  reorder();
  zone_to_field();
  build_least_squares();
  return 0;
}

int SolveControl::do_export()
{
  if(export_define.core)
  {
    sprintf(log_buf,"Export Core Data to %s...",export_define.CoreFile);GSS_LOG();
    //export mesh
    export_mesh(export_define.CoreFile);

    //export data
    for(int z=0;z<zone_num;z++)
      if(zonedata[z]->material_type==Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
        pzonedata->export_doping(export_define.CoreFile,&scale_unit);
        if(IsSingleCompSemiconductor(zone[z].zonelabel))
          if(pzonedata->export_mole(export_define.CoreFile,1)) return 1;
        pzonedata->export_solution(export_define.CoreFile,"solution",bc,&scale_unit);
      }
      else if(zonedata[z]->material_type==Insulator)
      {
        ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
        pzonedata->export_solution(export_define.CoreFile,"solution",bc,&scale_unit);
      }
      else if(zonedata[z]->material_type==Conductor)
      {
        ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
        pzonedata->export_solution(export_define.CoreFile,"solution",bc,&scale_unit);
      }
       sprintf(log_buf,"done\n");GSS_LOG();
  }

  if(export_define.ascii)
  {
    sprintf(log_buf,"Export Data to TIF File %s...",export_define.AscFile);GSS_LOG();
    extract_ascii(export_define.AscFile);
    sprintf(log_buf,"done\n");GSS_LOG();
  }

  if(export_define.vtk)
  {
    sprintf(log_buf,"Export Data to VTK File %s...",export_define.VTKFile);GSS_LOG();
    vtk_output_file(export_define.VTKFile);
    sprintf(log_buf,"done\n");GSS_LOG();
  }
  
  if(export_define.acceptor)
  {
    /* remove old file if exist*/
    remove(export_define.AcceptorFile);
    for(int z=0;z<zone_num;z++)
      if(zonedata[z]->material_type==Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
        pzonedata->export_accptor_ascii(export_define.AcceptorFile,&scale_unit);
      }
  }
  
  if(export_define.donor)
  {
    /* remove old file if exist*/
    remove(export_define.DonorFile);
    for(int z=0;z<zone_num;z++)
      if(zonedata[z]->material_type==Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
        pzonedata->export_donor_ascii(export_define.DonorFile,&scale_unit);
      }
  }
  
  return 0;
}


int  SolveControl::do_import()
{

  //delete old solution data
  bc.clear();
  clear_mesh();
  clear_data();
  if(import_define.file_type==MODEL)
  {
    sprintf(log_buf,"Import ModelData from %s...\n",import_define.ModelFile);GSS_LOG();
    if(import_cgns(import_define.ModelFile)) return 1;
    build_zone();
    if(setup_bc())   return 1;
    if(build_zonedata())    return 1;
    if(doping_func.size())
      setup_doping();
    else
      if(import_doping_from_cgns(import_define.ModelFile)) return 1;
    if(import_mole_from_cgns(import_define.ModelFile)) return 1;
    if(setup_init_data()) return 1;
    //reorder the mesh, important for implicit solver
    reorder();
    zone_to_field();
    //for fast call to least-squares arithmetic
    build_least_squares();
  }
  else
  {
    sprintf(log_buf,"Import CoreData from %s...\n",import_define.CoreFile);GSS_LOG();
    //init new solver
    if(import_cgns(import_define.CoreFile)) return 1;
    build_zone();
    if(setup_bc())   return 1;
    if(build_zonedata())    return 1;
    if(doping_func.size())
      setup_doping();
    else
      if(import_doping_from_cgns(import_define.CoreFile)) return 1;
    if(import_mole_from_cgns(import_define.CoreFile)) return 1;
    if(setup_init_data()) return 1;

    //import data
    for(int z=0;z<zone_num;z++)
    {
      if(zonedata[z]->material_type==Semiconductor)
      {
        SMCZone *pzonedata = dynamic_cast< SMCZone * >(zonedata[z]);
        if(pzonedata->import_solution(import_define.CoreFile,"solution",bc,&scale_unit)) return 1;
      }
      else if(zonedata[z]->material_type==Insulator)
      {
        ISZone *pzonedata = dynamic_cast< ISZone * >(zonedata[z]);
        if(pzonedata->import_solution(import_define.CoreFile,"solution",bc,&scale_unit)) return 1;
      }
      else if(zonedata[z]->material_type==Conductor)
      {
        ElZone *pzonedata = dynamic_cast< ElZone * >(zonedata[z]);
        if(pzonedata->import_solution(import_define.CoreFile,"solution",bc,&scale_unit)) return 1;
      }
    }
    //reorder the mesh, important for implicit solver
    reorder();
    zone_to_field();
    build_least_squares();
  }
  sprintf(log_buf,"\n");GSS_LOG();
  return 0;

}


int  SolveControl::do_attach()
{
  if(attach_define.electrode_type==VoltageBC)
  {
    vector<VSource *> bc_vsrc;
    for(int i=0;i<attach_define.vsrc_index.size();i++)
    {
      int vsrc_index=attach_define.vsrc_index[i];
      bc_vsrc.push_back(vsrc[vsrc_index]);
      sprintf(log_buf,"Vsource %s attached to bc %s\n",vsrc[vsrc_index]->label,attach_define.Electrode);
      GSS_LOG();
    }
    bc.Attach_Vapp(attach_define.Electrode,bc_vsrc);
  }
  if(attach_define.electrode_type==CurrentBC)
  {
    vector<ISource *> bc_isrc;
    for(int i=0;i<attach_define.isrc_index.size();i++)
    {
      int isrc_index=attach_define.isrc_index[i];
      bc_isrc.push_back(isrc[isrc_index]);
      sprintf(log_buf,"Isource %s attached to bc %s\n",isrc[isrc_index]->label,attach_define.Electrode);
      GSS_LOG();
    }
    bc.Attach_Iapp(attach_define.Electrode,bc_isrc);
  }
  return 0;
}

int  SolveControl::do_plot()
{
#ifdef HAVE_WIN32
  pid_t pid = fork();

  if (pid == (pid_t) 0)
  {
    // The child process to do plot.
    if(plot_define.MeshOnly)
    {
      if(plot_define.GeneratePS)
        plot_mesh_ps(plot_define);
      else
        plot_mesh_screen(plot_define);
    }
    else
      plot_data(plot_define);
    exit(0);
  }
  else if (pid < (pid_t) 0)
  {
    // The fork failed.
    fprintf (stderr, "Fork failed.\n");
  }
#endif

#ifdef HAVE_X11
  if(plot_define.MeshOnly)
  {
    /*
    plot_mesh_screen(plot_define);
    if(plot_define.GeneratePS)
      plot_mesh_ps(plot_define);
    */
    //for X11 exist, use show mesh instead
    show_mesh(plot_define);  
  }
  else
    plot_data(plot_define);
#endif
  return 0;
}


int  SolveControl::do_plotmesh()
{
#ifdef HAVE_X11
  show_mesh(plot_define);
#endif
  return 0;
}


int  SolveControl::do_vtk_plot()
{
#ifdef HAVE_VTK
  pid_t pid = fork();
  if (pid == (pid_t) 0)
  {
    vtk_display_data(plot_define);
    exit(0);
  }
  else if (pid < (pid_t) 0)
  {
    // The fork failed.
    fprintf (stderr, "Fork failed.\n");
  }
#endif
  return 0;
}

//------------------------------------------------------------------------------------


int  SolveControl::set_method(list<Cmd>::iterator  pcmd)
{
  PetscTruth flg;

  if(!pcmd->allowed_args(27,"type","scheme","ns","damping","ls","serverport","ejmodel","fermi","projection",
                         "highfieldmobility","impactionization","ii.type","bandbandtunneling","qnfactor","qpfactor",
                         "relative.tol","possion.tol","elec.continuty.tol","hole.continuty.tol","latt.temp.tol",
                         "elec.energy.tol","hole.energy.tol","elec.quantum.tol","hole.quantum.tol",
			 "electrode.tol","toler.relax","maxiteration"))
  {
    sprintf(log_buf,"line %d METHOD: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  solve_define.HighFieldMobility = pcmd->get_bool("highfieldmobility", 0, true);
  solve_define.BandBandTunneling = pcmd->get_bool("bandbandtunneling", 0, false);
  solve_define.QNFactor          = pcmd->get_number("qnfactor",0,1.0);
  solve_define.QPFactor          = pcmd->get_number("qpfactor",0,1.0);
  solve_define.EJModel           = pcmd->get_bool("ejmodel", 0, false);
  solve_define.Fermi             = pcmd->get_bool("fermi", 0, false);
  solve_define.Projection        = pcmd->get_bool("projection", 0, false);
  
  solve_define.maxit                    = pcmd->get_integer("maxiteration", 0, 30);
  solve_define.relative_toler           = pcmd->get_number("relative.tol",0,1e-5);
  solve_define.toler_relax              = pcmd->get_number("toler.relax",0,1e4);
  solve_define.possion_abs_toler        = pcmd->get_number("possion.tol",0,1e-29)*scale_unit.s_coulomb/scale_unit.s_micron;
  solve_define.elec_continuty_abs_toler = pcmd->get_number("elec.continuty.tol",0,5e-18)*scale_unit.s_A/scale_unit.s_micron;
  solve_define.hole_continuty_abs_toler = pcmd->get_number("hole.continuty.tol",0,5e-18)*scale_unit.s_A/scale_unit.s_micron;
  solve_define.heat_equation_abs_toler  = pcmd->get_number("latt.temp.tol",0,1e-11)*scale_unit.s_W/scale_unit.s_micron;
  solve_define.elec_energy_abs_toler    = pcmd->get_number("elec.energy.tol",0,1e-18)*scale_unit.s_W/scale_unit.s_micron;
  solve_define.hole_energy_abs_toler    = pcmd->get_number("hole.energy.tol",0,1e-18)*scale_unit.s_W/scale_unit.s_micron;
  solve_define.electrode_abs_toler      = pcmd->get_number("electrode.tol",0,1e-9)*scale_unit.s_volt;
  solve_define.elec_quantum_abs_toler   = pcmd->get_number("elec.quantum.tol",0,1e-29)*scale_unit.s_coulomb/scale_unit.s_micron;
  solve_define.hole_quantum_abs_toler   = pcmd->get_number("hole.quantum.tol",0,1e-29)*scale_unit.s_coulomb/scale_unit.s_micron;
  
  PetscOptionsGetInt(NULL, "-serverport",&solve_define.port,&flg);
  if(!flg) solve_define.port = pcmd->get_integer("serverport", 0, 1611);

  if(pcmd->is_arg_exist("type"))
  {
    if (pcmd->is_arg_value("type", "DDML1"))                 solve_define.Solver = DDML1E;
    else if (pcmd->is_arg_value("type", "DDML2"))            solve_define.Solver = DDML2E;
    else if (pcmd->is_arg_value("type", "DDML1E"))           solve_define.Solver = DDML1E;
    else if (pcmd->is_arg_value("type", "DDML2E"))           solve_define.Solver = DDML2E;
    else if (pcmd->is_arg_value("type", "DDML1AC"))          solve_define.Solver = DDML1AC;
    else if (pcmd->is_arg_value("type", "QDDML1E"))          solve_define.Solver = QDDML1E;
    else if (pcmd->is_arg_value("type", "DDML1MIX"))         solve_define.Solver = DDML1MIX;
    else if (pcmd->is_arg_value("type", "DDML2MIX"))         solve_define.Solver = DDML2MIX;
    else if (pcmd->is_arg_value("type", "EBML3"))            solve_define.Solver = EBML3E;
    else if (pcmd->is_arg_value("type", "EBML3E"))           solve_define.Solver = EBML3E;
    else if (pcmd->is_arg_value("type", "EMFEM"))            solve_define.Solver = EMFEM;
    else {sprintf(log_buf,"line %d METHOD: Type syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }
  
  solve_define.ImpactIonization  = pcmd->get_bool("impactionization", 0, false);
  if(solve_define.EJModel)
    solve_define.IIType = EdotJ;
  else if(solve_define.Solver == EBML3E)
    solve_define.IIType = ESide;
  else
    solve_define.IIType = GradQf;
  if(pcmd->is_arg_exist("ii.type"))
  {
    if (pcmd->is_arg_value("ii.type", "EdotJ"))                  solve_define.IIType = EdotJ;
    else if (pcmd->is_arg_value("ii.type", "GradQf"))            solve_define.IIType = GradQf;
    else if (pcmd->is_arg_value("ii.type", "EVector"))           solve_define.IIType = EVector;
    else if (pcmd->is_arg_value("ii.type", "ESide"))             solve_define.IIType = ESide;
    else if (pcmd->is_arg_value("ii.type", "Soft.II"))           solve_define.IIType = SoftII;
    else if (pcmd->is_arg_value("ii.type", "Hard.II"))           solve_define.IIType = HardII;
    else if (pcmd->is_arg_value("ii.type", "Temp.II"))           solve_define.IIType = TempII;
  else {sprintf(log_buf,"line %d METHOD: II.Type syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  

  if(pcmd->is_arg_exist("scheme"))
  {
    if (pcmd->is_arg_value("scheme", "Newton"))              solve_define.Scheme = DDM_Newton;
    else if (pcmd->is_arg_value("scheme", "Gummel"))         solve_define.Scheme = DDM_Gummel;
  else {sprintf(log_buf,"line %d METHOD: Scheme syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("ns"))
  {
    if (pcmd->is_arg_value("ns", "LineSearch"))              solve_define.NS = LineSearch;
    else if (pcmd->is_arg_value("ns", "Basic"))              solve_define.NS = Basic;
    else if (pcmd->is_arg_value("ns", "TrustRegion"))        solve_define.NS = TrustRegion;
    else {sprintf(log_buf,"line %d METHOD: NS syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("ls"))
  {
    char * linear_solver = pcmd->get_string("ls",0,"lu");
    for(int c=0;c<strlen(linear_solver);c++)
      if(isupper(linear_solver[c])) linear_solver[c]=tolower(linear_solver[c]);
    solve_define.LS = linear_solver;
  }

  if(pcmd->is_arg_exist("damping"))
  {
    if (pcmd->is_arg_value("damping", "No"))                    solve_define.Damping = DampingNo;
    else if (pcmd->is_arg_value("damping", "BankRose"))         solve_define.Damping = DampingBankRose;
    else if (pcmd->is_arg_value("damping", "Potential"))        solve_define.Damping = DampingPotential;
    else {sprintf(log_buf,"line %d METHOD: Damping syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }
  
  switch(solve_define.Solver)
  {
  case DDML1:   sprintf(log_buf,"Solver type       : DDM Level 1\n");  GSS_LOG();break;
  case DDML1E:  sprintf(log_buf,"Solver type       : DDM Level 1E\n"); GSS_LOG();break;
  case DDML1AC: sprintf(log_buf,"Solver type       : DDM Level 1E AC Sweep\n"); GSS_LOG();return 0;
  case QDDML1E: sprintf(log_buf,"Solver type       : DDM Level 1E with Quantum Correction\n"); GSS_LOG();break;
  case DDML1MIX:sprintf(log_buf,"Solver type       : DDM Level 1E Device/Circuit Mixed\n"); GSS_LOG();break;
  case DDML2:   sprintf(log_buf,"Solver type       : DDM Level 2\n");  GSS_LOG();break;
  case DDML2E:  sprintf(log_buf,"Solver type       : DDM Level 2E\n"); GSS_LOG();break;
  case DDML2MIX:sprintf(log_buf,"Solver type       : DDM Level 2E Device/Circuit Mixed\n"); GSS_LOG();break;
  case EBML3E:  sprintf(log_buf,"Solver type       : EBM Level 3E\n"); GSS_LOG();break;
  case EMFEM:   sprintf(log_buf,"Solver type       : EMFEM\n");        GSS_LOG();return 0;
  case POISSON: sprintf(log_buf,"Solver type       : POISSON\n");      GSS_LOG();return 0;
  default: break;
  }


  switch(solve_define.Scheme)
  {
  case DDM_Newton:   sprintf(log_buf,"Numerical scheme  : Newton Iteration\n");GSS_LOG();break;
  default: break;
  }

  switch(solve_define.NS)
  {
  case LineSearch:   sprintf(log_buf,"Nonlinear method  : LineSearch\n");GSS_LOG();break;
  case Basic:        sprintf(log_buf,"Nonlinear method  : Basic\n");GSS_LOG();break;
  case TrustRegion:  sprintf(log_buf,"Nonlinear method  : TrustRegion\n");GSS_LOG();break;
  default: break;
  }
  sprintf(log_buf, "Linear Solver     : %s\n",solve_define.LS.c_str());GSS_LOG();

  if(solve_define.ImpactIonization == true)
  {
    sprintf(log_buf,"Impact Ionization enabled\n");GSS_LOG();
  }
  else
  {
    sprintf(log_buf,"Impact Ionization disabled\n");GSS_LOG();
  }

  if(solve_define.BandBandTunneling == true)
  {
    sprintf(log_buf,"Band to Band Tunneling enabled\n");GSS_LOG();
  }
  else
  {
    sprintf(log_buf,"Band to Band Tunneling disabled\n");GSS_LOG();
  }

  return 0;
}


int  SolveControl::set_solve(list<Cmd>::iterator  pcmd)
{
  //do clear
  solve_define.Type  = 0;
  solve_define.Electrode_Record_Name.clear();
  solve_define.Electrode_Record_Index.clear();
  solve_define.Electrode_Record.clear();
  solve_define.Electrode_VScan_Name.clear();
  solve_define.IVFile.erase();
  solve_define.BDF_Type = BDF2;

  solve_define.VStart    = pcmd->get_number("vstart",0,0.0)*scale_unit.s_volt;
  solve_define.VStep     = pcmd->get_number("vstep",0,0.0)*scale_unit.s_volt;
  solve_define.VStepMax  = pcmd->get_number("vstepmax",0,10*solve_define.VStep)*scale_unit.s_volt;
  solve_define.VStop     = pcmd->get_number("vstop",0,1e10)*scale_unit.s_volt;
  solve_define.IStart    = pcmd->get_number("istart",0,0.0)*scale_unit.s_mA;
  solve_define.IStep     = pcmd->get_number("istep",0,0.0)*scale_unit.s_mA;
  solve_define.IStop     = pcmd->get_number("istop",0,1e10)*scale_unit.s_mA;
  solve_define.TStart    = pcmd->get_number("tstart",0,0.0)*scale_unit.s_second;
  solve_define.TStep     = pcmd->get_number("tstep",0,0.0)*scale_unit.s_second;
  solve_define.TStop     = pcmd->get_number("tstop",0,1e10)*scale_unit.s_second;
  solve_define.VAC       = pcmd->get_number("vac",0,0.0026/scale_unit.s_volt)*scale_unit.s_volt;
  solve_define.FStart    = pcmd->get_number("fstart",0,1e6)/scale_unit.s_second;
  solve_define.FMultiple = pcmd->get_number("fmultiple",0,1.1);
  solve_define.FStop     = pcmd->get_number("fstop",0,1e9)/scale_unit.s_second;
  solve_define.AutoStep  = pcmd->get_bool("autostep",0,true);
  solve_define.Predict   = pcmd->get_bool("predict",0,true);

  solve_define.IVFile = pcmd->get_string("ivfile",0,"");
  solve_define.IVFileAppend = pcmd->get_bool("append",0,false);

  if(pcmd->is_arg_exist("type"))
  {
    if (pcmd->is_arg_value("type",      "EQUILIBRIUM"))       solve_define.Type = EQUILIBRIUM;
    else if (pcmd->is_arg_value("type", "STEADYSTATE"))       solve_define.Type = STEADYSTATE;
    else if (pcmd->is_arg_value("type", "DCSWEEP"))           solve_define.Type = DCSWEEP;
    else if (pcmd->is_arg_value("type", "ACSWEEP"))           solve_define.Type = ACSWEEP;
    else if (pcmd->is_arg_value("type", "TRANSIENT"))         solve_define.Type = TRANSIENT;
    else if (pcmd->is_arg_value("type", "TRACE"))             solve_define.Type = TRACE;
  else {sprintf(log_buf,"line %d SOLVE: Type syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("ode.formula"))
  {
    if (pcmd->is_arg_value("ode.formula", "BDF1"))       solve_define.BDF_Type = BDF1;
    else if (pcmd->is_arg_value("ode.formula", "BDF2"))  solve_define.BDF_Type = BDF2;
  else {sprintf(log_buf,"line %d SOLVE: ODE.Formula syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  for(pcmd->arg_search_begin();!pcmd->arg_search_end();pcmd->goto_next_arg())
  {
    if(!strcmp(pcmd->get_current_arg_label(),"vscan"))
    {
      string vscan_name=pcmd->get_current_arg_value().sval;
      solve_define.Electrode_VScan_Name.push_back(vscan_name);
      solve_define.Electrode_VScan = 1;
      solve_define.Electrode_IScan =-1;
      if(!bc.is_electrode(vscan_name.c_str()))
      {sprintf(log_buf,"line %d SOLVE: Can't find VScan Electrode!!\n",pcmd->lineno);GSS_LOG();return 1;}
    }
  }  

  if(pcmd->is_arg_exist("iscan"))
  {
    solve_define.Electrode_IScan_Name = pcmd->get_string("iscan",0,"");
    solve_define.Electrode_VScan = -1;
    solve_define.Electrode_IScan =  1;
    if(!bc.is_electrode(pcmd->get_string("iscan",0,"")))
    {sprintf(log_buf,"line %d SOLVE: Can't find IScan Electrode!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("acscan"))
  {
    solve_define.Electrode_ACScan_Name = pcmd->get_string("acscan",0,"");
    if(!bc.is_electrode(pcmd->get_string("acscan",0,"")))
    {sprintf(log_buf,"line %d SOLVE: Can't find ACScan Electrode!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  for(pcmd->arg_search_begin();!pcmd->arg_search_end();pcmd->goto_next_arg())
  {
    if(!strcmp(pcmd->get_current_arg_label(),"ivrecord"))
    {
      if(!bc.is_electrode(pcmd->get_current_arg_value().sval))
      {sprintf(log_buf,"line %d SOLVE: Can't find IVRecord Electrode!!\n",pcmd->lineno);GSS_LOG();return 1;}
      solve_define.Electrode_Record_Name.push_back(pcmd->get_current_arg_value().sval);
    }
  }
  return 0;
}


int  SolveControl::set_export(list<Cmd>::iterator  pcmd)
{
  //do clear
  export_define.core  = false;
  export_define.ascii = false;
  export_define.vtk   = false;
  export_define.acceptor = false;
  export_define.donor    = false;

  if(!pcmd->allowed_args(5,"corefile","ascfile","vtkfile","acceptorfile","donorfile"))
  {
    sprintf(log_buf,"line %d EXPORT: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  if(pcmd->is_arg_exist("corefile"))
  {
    export_define.core = true;
    strcpy(export_define.CoreFile,pcmd->get_string("corefile",0,""));
  }

  if(pcmd->is_arg_exist("ascfile"))
  {
    export_define.ascii = true;
    strcpy(export_define.AscFile,pcmd->get_string("ascfile",0,""));
  }

  if(pcmd->is_arg_exist("vtkfile"))
  {
    export_define.vtk = true;
    strcpy(export_define.VTKFile,pcmd->get_string("vtkfile",0,""));
  }
  
  if(pcmd->is_arg_exist("acceptorfile"))
  {
    export_define.acceptor = true;
    strcpy(export_define.AcceptorFile,pcmd->get_string("acceptorfile",0,""));
  }
  
  if(pcmd->is_arg_exist("donorfile"))
  {
    export_define.donor = true;
    strcpy(export_define.DonorFile,pcmd->get_string("donorfile",0,""));
  }
  return 0;
}


int  SolveControl::set_refine(list<Cmd>::iterator  pcmd)
{
  refine_define.Measure = Linear;
  refine_define.Dispersion = pcmd->get_number("dispersion",0,3.0);
  refine_define.DivisionRatio = pcmd->get_number("divisionratio",0,0.25);
  strcpy(refine_define.tri_cmd, pcmd->get_string("triangle", 0, "praq30DzQ"));

  if(!pcmd->allowed_args(5,"variable","measure","dispersion","divisionratio","triangle"))
  {
    sprintf(log_buf,"line %d REFINE: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  if(pcmd->is_arg_exist("variable"))
  {
    if (pcmd->is_arg_value("variable", "Potential"))         refine_define.Variable = Potential;
    else if (pcmd->is_arg_value("variable", "Doping"))       refine_define.Variable = Doping;
  else {sprintf(log_buf,"line %d REFINE: Variable syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("measure"))
  {
    if (pcmd->is_arg_value("measure", "Linear"))              refine_define.Measure = Linear;
    else if (pcmd->is_arg_value("measure", "SignedLog"))      refine_define.Measure = SignedLog;
  else {sprintf(log_buf,"line %d REFINE: Measure syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  return 0;
}


int  SolveControl::set_import(list<Cmd>::iterator  pcmd)
{
  //do clear
  import_define.file_type = 0;

  if(!pcmd->allowed_args(2,"corefile","modelfile"))
  {
    sprintf(log_buf,"line %d IMPORT: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  if(pcmd->is_arg_exist("corefile"))
  {
    import_define.file_type = CORE;
    strcpy(import_define.CoreFile,pcmd->get_string("corefile",0,""));
  }

  if(pcmd->is_arg_exist("modelfile"))
  {
    import_define.file_type = MODEL;
    strcpy(import_define.ModelFile,pcmd->get_string("modelfile",0,""));
  }

  return 0;
}


int  SolveControl::set_attach(list<Cmd>::iterator  pcmd)
{
  //do clear
  attach_define.vsrc_name.clear();
  attach_define.vsrc_index.clear();
  attach_define.isrc_name.clear();
  attach_define.isrc_index.clear();
  attach_define.electrode_type = VoltageBC;

  if(!pcmd->allowed_args(4,"electrode","vapp","iapp","type"))
  {
    sprintf(log_buf,"line %d ATTACH: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  if(pcmd->is_arg_exist("electrode"))
  {
    strcpy(attach_define.Electrode,pcmd->get_string("electrode",0,""));
    if(!bc.is_electrode(attach_define.Electrode))
    {
      sprintf(log_buf,"line %d ATTACH: Can't find Electrode %s!\n",pcmd->lineno,attach_define.Electrode);
      GSS_LOG();
      return 1;
    }
  }

  if(pcmd->is_arg_exist("type"))
  {
    if (pcmd->is_arg_value("type", "Current"))       attach_define.electrode_type = CurrentBC;
    else if(pcmd->is_arg_value("type", "Voltage"))   attach_define.electrode_type = VoltageBC;
  else {sprintf(log_buf,"line %d ATTACH: No such electrode type !\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  for(pcmd->arg_search_begin();!pcmd->arg_search_end();pcmd->goto_next_arg())
  {
    if(!strcmp(pcmd->get_current_arg_label(),"vapp"))
    {
      string vsrc_name=pcmd->get_current_arg_value().sval;
      attach_define.vsrc_name.push_back(vsrc_name);
      int flag=0;
      for(int j=0;j<vsrc.size();j++)
        if(!strcmp(vsrc_name.c_str(),vsrc[j]->label))
        {
          attach_define.vsrc_index.push_back(j);
          flag=1;
        }
      if(flag==0)
      {
        sprintf(log_buf,"line %d ATTACH:Can't find VApp %s!\n",pcmd->lineno,pcmd->get_current_arg_value().sval);
        GSS_LOG();
        return 1;
      }

    }
    else if(!strcmp(pcmd->get_current_arg_label(),"iapp"))
    {
      string isrc_name=pcmd->get_current_arg_value().sval;
      attach_define.isrc_name.push_back(isrc_name);
      int flag=0;
      for(int j=0;j<isrc.size();j++)
        if(!strcmp(isrc_name.c_str(),isrc[j]->label))
        {
          attach_define.isrc_index.push_back(j);
          flag=1;
        }
      if(flag==0)
      {
        sprintf(log_buf,"line %d ATTACH:Can't find IApp %s!\n",pcmd->lineno,pcmd->get_current_arg_value().sval);
        GSS_LOG();
        return 1;
      }
    }
  }

  return 0;

}


int  SolveControl::set_plot(list<Cmd>::iterator  pcmd)
{

  //set default value
  plot_define.MeshOnly = 0;
  plot_define.Resolution = RES_800x600;
  plot_define.Measure = Linear;
  plot_define.GeneratePS = 0;
  plot_define.GenerateTIFF = 0;
  plot_define.Style = F_COLOR;
  plot_define.az = pcmd->get_number("azangle",0,240.0);
  plot_define.el = pcmd->get_number("elangle",0,60.0);
  plot_define.Persp = 0;
  plot_define.ztick = -1;

  if(!pcmd->allowed_args(8,"variable","resolution","ps.out","tiff.out","measure","azangle","elangle","style"))
  {
    sprintf(log_buf,"line %d PLOT: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  if(pcmd->is_arg_exist("variable"))
  {
    plot_define.VariableName=pcmd->get_string("variable",0,"");
    if (pcmd->is_arg_value("variable", "Potential"))          plot_define.Variable = Potential;
    else if (pcmd->is_arg_value("variable", "Doping"))        plot_define.Variable = Doping;
    else if (pcmd->is_arg_value("variable", "NetDoping"))     plot_define.Variable = NetDoping;
    else if (pcmd->is_arg_value("variable", "Na"))            plot_define.Variable = DopingNa;
    else if (pcmd->is_arg_value("variable", "Nd"))            plot_define.Variable = DopingNd;
    else if (pcmd->is_arg_value("variable", "Phosphorus"))    plot_define.Variable = Phosphorus;
    else if (pcmd->is_arg_value("variable", "Arsenic"))       plot_define.Variable = Arsenic;
    else if (pcmd->is_arg_value("variable", "Antimony"))      plot_define.Variable = Antimony;
    else if (pcmd->is_arg_value("variable", "Boron"))         plot_define.Variable = Boron;
    else if (pcmd->is_arg_value("variable", "ElecDensity"))   plot_define.Variable = ElecDensity;
    else if (pcmd->is_arg_value("variable", "HoleDensity"))   plot_define.Variable = HoleDensity;
    else if (pcmd->is_arg_value("variable", "ElecTemp"))      plot_define.Variable = ElecTemp;
    else if (pcmd->is_arg_value("variable", "HoleTemp"))      plot_define.Variable = HoleTemp;
    else if (pcmd->is_arg_value("variable", "Phi.Intrinsic")) plot_define.Variable = Phi_Intrinsic;
    else if (pcmd->is_arg_value("variable", "PhiN"))          plot_define.Variable = PhiN;
    else if (pcmd->is_arg_value("variable", "PhiP"))          plot_define.Variable = PhiP;
    else if (pcmd->is_arg_value("variable", "QuantumEc"))     plot_define.Variable = QuantumEc;
    else if (pcmd->is_arg_value("variable", "QuantumEv"))     plot_define.Variable = QuantumEv;
    else if (pcmd->is_arg_value("variable", "Temperature"))   plot_define.Variable = Temperature;
    else if (pcmd->is_arg_value("variable", "EFieldX"))       plot_define.Variable = EFieldX;
    else if (pcmd->is_arg_value("variable", "EFieldY"))       plot_define.Variable = EFieldY;
    else if (pcmd->is_arg_value("variable", "E.Field"))       plot_define.Variable = EField;
    else if (pcmd->is_arg_value("variable", "OpticalEx"))     plot_define.Variable = OpticalEx;//```
    else if (pcmd->is_arg_value("variable", "OpticalEy"))     plot_define.Variable = OpticalEy;
    else if (pcmd->is_arg_value("variable", "OpticalEz"))     plot_define.Variable = OpticalEz;
    else if (pcmd->is_arg_value("variable", "OpticalHx"))     plot_define.Variable = OpticalHx;
    else if (pcmd->is_arg_value("variable", "OpticalHy"))     plot_define.Variable = OpticalHy;
    else if (pcmd->is_arg_value("variable", "OpticalHz"))     plot_define.Variable = OpticalHz;
    else if (pcmd->is_arg_value("variable", "OpticalG"))      plot_define.Variable = OpticalG;
    else if (pcmd->is_arg_value("variable", "Mesh"))          plot_define.MeshOnly = MeshObject;
    else if (pcmd->is_arg_value("variable", "DeviceMesh"))    plot_define.MeshOnly = MeshObject;
  else {sprintf(log_buf,"line %d PLOT: Variable syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("resolution"))
  {
    if (pcmd->is_arg_value("resolution", "RES.Low"))         plot_define.Resolution = RES_640x480;
    else if (pcmd->is_arg_value("resolution", "RES.Middle")) plot_define.Resolution = RES_800x600;
    else if (pcmd->is_arg_value("resolution", "RES.High"))   plot_define.Resolution = RES_1024x768;
  else {sprintf(log_buf,"line %d PLOT: Resolution syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("ps.out"))
  {
    plot_define.GeneratePS = 1;
    plot_define.PSFileName=pcmd->get_string("ps.out",0,"");
  }

  if(pcmd->is_arg_exist("tiff.out"))
  {
    plot_define.GenerateTIFF = 1;
    plot_define.TIFFFileName=pcmd->get_string("tiff.out",0,"");
  }

  if(pcmd->is_arg_exist("measure"))
  {
    if (pcmd->is_arg_value("measure", "Linear"))         plot_define.Measure = Linear;
    else if (pcmd->is_arg_value("measure", "SignedLog")) plot_define.Measure = SignedLog;
  else {sprintf(log_buf,"line %d PLOT: Measure syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("style"))
  {
    if (pcmd->is_arg_value("style", "Scale"))           plot_define.Style = F_SCALE;
    else if (pcmd->is_arg_value("style", "Color"))      plot_define.Style = F_COLOR;
    else if (pcmd->is_arg_value("style", "GrayLevel"))  plot_define.Style = F_GRAYLEVEL;
  else {sprintf(log_buf,"line %d PLOT: Style syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  return 0;
}


int  SolveControl::set_plotmesh(list<Cmd>::iterator  pcmd)
{
#ifdef HAVE_X11
  plot_define.GenerateTIFF = 0;
  if(!pcmd->allowed_args(1,"tiff.out"))
  {
    sprintf(log_buf,"line %d PLOTMESH: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }
  if(pcmd->is_arg_exist("tiff.out"))
  {
    plot_define.GenerateTIFF = 1;
    plot_define.TIFFFileName=pcmd->get_string("tiff.out",0,"");
  }
#endif  
  return 0;
}

int  SolveControl::set_vtk_plot(list<Cmd>::iterator  pcmd)
{

  //set default value
  plot_define.MeshOnly = 0;
  plot_define.Resolution = RES_800x600;
  plot_define.Measure = Linear;

  if(!pcmd->allowed_args(3,"variable","resolution","measure"))
  {
    sprintf(log_buf,"line %d PLOTVTK: unrecognized parameter(s)!\n",pcmd->get_current_lineno());
    GSS_LOG();
    return 1;
  }

  if(pcmd->is_arg_exist("variable"))
  {
    plot_define.VariableName=pcmd->get_string("variable",0,"");
    if (pcmd->is_arg_value("variable", "Potential"))          plot_define.Variable = Potential;
    else if (pcmd->is_arg_value("variable", "Doping"))        plot_define.Variable = Doping;
    else if (pcmd->is_arg_value("variable", "NetDoping"))     plot_define.Variable = NetDoping;
    else if (pcmd->is_arg_value("variable", "Na"))            plot_define.Variable = DopingNa;
    else if (pcmd->is_arg_value("variable", "Nd"))            plot_define.Variable = DopingNd;
    else if (pcmd->is_arg_value("variable", "Phosphorus"))    plot_define.Variable = Phosphorus;
    else if (pcmd->is_arg_value("variable", "Arsenic"))       plot_define.Variable = Arsenic;
    else if (pcmd->is_arg_value("variable", "Antimony"))      plot_define.Variable = Antimony;
    else if (pcmd->is_arg_value("variable", "Boron"))         plot_define.Variable = Boron;
    else if (pcmd->is_arg_value("variable", "ElecDensity"))   plot_define.Variable = ElecDensity;
    else if (pcmd->is_arg_value("variable", "HoleDensity"))   plot_define.Variable = HoleDensity;
    else if (pcmd->is_arg_value("variable", "ElecTemp"))      plot_define.Variable = ElecTemp;
    else if (pcmd->is_arg_value("variable", "HoleTemp"))      plot_define.Variable = HoleTemp;
    else if (pcmd->is_arg_value("variable", "Phi.Intrinsic")) plot_define.Variable = Phi_Intrinsic;
    else if (pcmd->is_arg_value("variable", "PhiN"))          plot_define.Variable = PhiN;
    else if (pcmd->is_arg_value("variable", "PhiP"))          plot_define.Variable = PhiP;
    else if (pcmd->is_arg_value("variable", "QuantumEc"))     plot_define.Variable = QuantumEc;
    else if (pcmd->is_arg_value("variable", "QuantumEv"))     plot_define.Variable = QuantumEv;
    else if (pcmd->is_arg_value("variable", "Temperature"))   plot_define.Variable = Temperature;
    else if (pcmd->is_arg_value("variable", "EFieldX"))       plot_define.Variable = EFieldX;
    else if (pcmd->is_arg_value("variable", "EFieldY"))       plot_define.Variable = EFieldY;
    else if (pcmd->is_arg_value("variable", "E.Field"))       plot_define.Variable = EField;
    else if (pcmd->is_arg_value("variable", "OpticalEx"))     plot_define.Variable = OpticalEx;//```
    else if (pcmd->is_arg_value("variable", "OpticalEy"))     plot_define.Variable = OpticalEy;
    else if (pcmd->is_arg_value("variable", "OpticalEz"))     plot_define.Variable = OpticalEz;
    else if (pcmd->is_arg_value("variable", "OpticalHx"))     plot_define.Variable = OpticalHx;
    else if (pcmd->is_arg_value("variable", "OpticalHy"))     plot_define.Variable = OpticalHy;
    else if (pcmd->is_arg_value("variable", "OpticalHz"))     plot_define.Variable = OpticalHz;
    else if (pcmd->is_arg_value("variable", "OpticalG"))      plot_define.Variable = OpticalG;
  else {sprintf(log_buf,"line %d PLOTVTK: Variable syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("measure"))
  {
    if (pcmd->is_arg_value("measure", "Linear"))         plot_define.Measure = Linear;
    else if (pcmd->is_arg_value("measure", "SignedLog")) plot_define.Measure = SignedLog;
  else {sprintf(log_buf,"line %d PLOTVTK: Measure syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }

  if(pcmd->is_arg_exist("resolution"))
  {
    if (pcmd->is_arg_value("resolution", "RES.Low"))         plot_define.Resolution = RES_640x480;
    else if (pcmd->is_arg_value("resolution", "RES.Middle")) plot_define.Resolution = RES_800x600;
    else if (pcmd->is_arg_value("resolution", "RES.High"))   plot_define.Resolution = RES_1024x768;
  else {sprintf(log_buf,"line %d PLOTVTK: Resolution syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  }
  return 0;
}


int  SolveControl::set_probe(list<Cmd>::iterator pcmd)
{
  ProbeDefine probe;
  probe.Region = pcmd->get_string("region", 0, "");
  probe.Segment = pcmd->get_string("segment", 0, "");
  probe.zone_index = Get_zone_index(probe.Region.c_str());
  if(probe.zone_index == -1)
    {sprintf(log_buf,"line %d PROBE: I can't find region for probing!\n",pcmd->lineno);GSS_LOG();return 1;}
  probe.bc_index = bc.Get_bc_index(probe.zone_index,probe.Segment.c_str());
  if(probe.bc_index == -1)
    {sprintf(log_buf,"line %d PROBE: I can't find segment for probing!\n",pcmd->lineno);GSS_LOG();return 1;}

  probe.ProbeFile = pcmd->get_string("probefile", 0, "");
  probe.Append = pcmd->get_bool("append", 0, false);

  probe.VariableName = pcmd->get_string("variable",0,"");
  if (pcmd->is_arg_value("variable", "Potential"))          probe.Variable = Potential;
  else if (pcmd->is_arg_value("variable", "Doping"))        probe.Variable = Doping;
  else if (pcmd->is_arg_value("variable", "NetDoping"))     probe.Variable = NetDoping;
  else if (pcmd->is_arg_value("variable", "Na"))            probe.Variable = DopingNa;
  else if (pcmd->is_arg_value("variable", "Nd"))            probe.Variable = DopingNd;
  else if (pcmd->is_arg_value("variable", "Phosphorus"))    probe.Variable = Phosphorus;
  else if (pcmd->is_arg_value("variable", "Arsenic"))       probe.Variable = Arsenic;
  else if (pcmd->is_arg_value("variable", "Antimony"))      probe.Variable = Antimony;
  else if (pcmd->is_arg_value("variable", "Boron"))         probe.Variable = Boron;
  else if (pcmd->is_arg_value("variable", "ElecDensity"))   probe.Variable = ElecDensity;
  else if (pcmd->is_arg_value("variable", "HoleDensity"))   probe.Variable = HoleDensity;
  else if (pcmd->is_arg_value("variable", "ElecTemp"))      probe.Variable = ElecTemp;
  else if (pcmd->is_arg_value("variable", "HoleTemp"))      probe.Variable = HoleTemp;
  else if (pcmd->is_arg_value("variable", "Phi.Intrinsic")) probe.Variable = Phi_Intrinsic;
  else if (pcmd->is_arg_value("variable", "PhiN"))          probe.Variable = PhiN;
  else if (pcmd->is_arg_value("variable", "PhiP"))          probe.Variable = PhiP;
  else if (pcmd->is_arg_value("variable", "QuantumEc"))     probe.Variable = QuantumEc;
  else if (pcmd->is_arg_value("variable", "QuantumEv"))     probe.Variable = QuantumEv;
  else if (pcmd->is_arg_value("variable", "Temperature"))   probe.Variable = Temperature;
  else if (pcmd->is_arg_value("variable", "EFieldX"))       probe.Variable = EFieldX;
  else if (pcmd->is_arg_value("variable", "EFieldY"))       probe.Variable = EFieldY;
  else if (pcmd->is_arg_value("variable", "E.Field"))       probe.Variable = EField;
  else if (pcmd->is_arg_value("variable", "OpticalEx"))     probe.Variable = OpticalEx;//```
  else if (pcmd->is_arg_value("variable", "OpticalEy"))     probe.Variable = OpticalEy;
  else if (pcmd->is_arg_value("variable", "OpticalEz"))     probe.Variable = OpticalEz;
  else if (pcmd->is_arg_value("variable", "OpticalHx"))     probe.Variable = OpticalHx;
  else if (pcmd->is_arg_value("variable", "OpticalHy"))     probe.Variable = OpticalHy;
  else if (pcmd->is_arg_value("variable", "OpticalHz"))     probe.Variable = OpticalHz;
  else if (pcmd->is_arg_value("variable", "OpticalG"))      probe.Variable = OpticalG;
  else {sprintf(log_buf,"line %d PROBE: Variable syntax error!\n",pcmd->lineno);GSS_LOG();return 1;}
  probe_define.push_back(probe);
  return 0;
}

