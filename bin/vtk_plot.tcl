package require vtk
package require vtkinteraction

#----- get startup directory and name -----------------------------------------
set cmd_name [file tail $argv0]
set cmd_dir  [file dirname $argv0]
if {![file exists $argv0] || [file isdirectory $argv0]} {
  if {$tcl_platform(platform) == "windows"} {
    set sep ";"
  } else {
    set sep ":"
  }
  foreach i [split $env(PATH) $sep] {
    if {$sep == ";"} {
      set i [join [split $i \\] /]
    }
    if {[file exists $i/$cmd_name] && ![file isdirectory $i/$cmd_name]} {
      set cmd_dir $i
      break;
    }
  }
}
set curdir [pwd]
if ![catch {cd $cmd_dir}] {
  set cmd_dir [pwd]
  cd $curdir
}
if {$tcl_platform(platform) == "windows"} {
  set cmd_dir [file attributes $cmd_dir -shortname]
}


#----- set path to tcl scripts ------------------------------------------------
if {[info exists env(GSS_DIR)]} {
  set auto_path "$env(GSS_DIR)/lib/gui_script $cmd_dir $cmd_dir/../lib/gui_script . $auto_path"
}

proc error_exit {msg} {
  wm withdraw .
  tk_dialog .error Error $msg error 0 Exit
  exit 1
}

if [catch {do_config -win} msg] {error_exit $msg}

#----- global variable array --------------------------------------------------
array set VTKPlot {
  vtk_data_file ""
  scalarname    ""
  view          "surface"
  min           0
  max           0
  view3D        0
}


vtkInteractorStyleTrackballCamera   style1
vtkUnstructuredGridReader  reader
vtkWarpScalar              rubberPlane
vtkConnectivityFilter      connect
vtkGeometryFilter          parison
vtkOutlineCornerFilter     bbox 
vtkPolyDataNormals         normals
vtkPolyDataMapper          parisonMapper
vtkDataSetMapper           moldMapper

vtkActor                   modelActor
vtkActor                   modelActor3D
vtkActor                   meshActor
vtkActor                   bboxActor
vtkCubeAxesActor2D         modelAxes
vtkAxesActor               axes
vtkLookupTable             lut
vtkOrientationMarkerWidget marker1
vtkTextProperty            tprop1
vtkTextProperty            tprop2
vtkTextProperty            tprop3
vtkTextProperty            tprop4
vtkScalarBarActor          scalarBar

#----- create main window -----------------------------------------------------

wm title . "GSS 3D plot Window"
wm minsize    .  640 480
wm protocol   . WM_DELETE_WINDOW exit

#----- create menu -------------------------------------------------------------
menubar_create {File Display Option Help}

#----- file menu ---------------------------------------------------------------
set m [menubar_get File]
$m add command -label "Save Screen..." -command save_image_file
$m add separator
$m add command -label "Quit" -command exit

proc save_image_file { } {
    # Save window to image
    global renWin
    set file_name [FileSave "Save Image File" "" . \
         {{{jpeg Files} {.png}} {{All Files} {*}}} png]
    if {$file_name == ""} return
  
    vtkWindowToImageFilter renderImage
    renderImage SetInput renWin
    
    vtkPNGWriter writer

    #vtkJPEGWriter writer
    #vtkMetaImageWriter writer
    #vtkPostScriptWriter writer
    #writer ProgressiveOn
    #writer SetQuality 100
    writer SetFileName $file_name
    writer SetInput [renderImage GetOutput]
    writer Write
    renderImage Delete
    writer Delete
}

#----- Option menu ---------------------------------------------------------------
set m [menubar_get Display]
$m add checkbutton -label "3D View" -variable VTKPlot(view3D) \
  -onvalue 1 -offvalue 0 -command {Plot_Update}
$m add separator
$m add command -label "Surface" -command surface_view
$m add command -label "Wire frame" -command wireframe_view
$m add command -label "Points" -command points_view

proc surface_view {} {
  global VTKPlot
  global ren modelActor
  set VTKPlot(view) "surface"
  ren RemoveActor modelActor
  Plot_Update
}

proc wireframe_view {} {
  global VTKPlot
  global ren modelActor
  set VTKPlot(view) "wireframe"
  ren RemoveActor modelActor
  Plot_Update
}

proc points_view {} {
  global VTKPlot
  global ren modelActor
  set VTKPlot(view) "points"
  ren RemoveActor modelActor
  Plot_Update
}

proc Plot_Update {} {
  global VTKPlot
  global renWin
  if {$VTKPlot(view3D)} {
    wrap_z
  }  else {
    plan_z
  }
  model_obj 
  draw_window 
  renWin Render  
}

#---------- images ------------------------------------------------------------
image create photo img_about -file $env(GSS_DIR)/lib/gss.logos.gif

#----- help menu --------------------------------------------------------------
set m [menubar_get Help]
$m add command -label "About..." -underline 0 -command do_about

proc do_about {} {
  global ProgData
  dialog .about -1 -1 "About GSS VTK Plot Window" \
 {    GSS VTK Window 0.1
 A GUI for post-process.
 Date   : 11/10/2006
 Author : Gong Ding
 Email  : gdiso@ustc.edu} \
 img_about 0 Close
}

#----- create VTK display window------------------------------------------------
vtkRenderer                ren
vtkRenderWindow            renWin
   renWin AddRenderer ren

frame .f1
pack .f1 -fill both -expand 1
set vtkw [vtkTkRenderWidget .f1.r1 -width 640 -height 480 -rw renWin]
pack $vtkw -side left -padx 3 -pady 3 -fill both -expand 1
::vtk::bind_tk_render_widget $vtkw
#vtkRenderWindowInteractor  iren
#   iren SetRenderWindow renWin
#   iren Initialize
#   iren AddObserver UserEvent {wm deiconify .vtkInteract}
#   iren SetInteractorStyle style1



proc read_data {} {
  global VTKPlot reader
  reader SetFileName $VTKPlot(vtk_data_file)
  reader SetScalarsName $VTKPlot(scalarname)
	reader Update
	set range [[reader GetOutput] GetScalarRange]
  set VTKPlot(min) [lindex $range 0]
  set VTKPlot(max) [lindex $range 1]        
}

proc wrap_z {} {
  # set the scale value to z. 
  global VTKPlot
  global reader rubberPlane connect parison normals
  rubberPlane SetInput [reader GetOutput]
  rubberPlane SetScaleFactor [expr 2.0/($VTKPlot(max)-$VTKPlot(min))]
  # The threshold filter has been used to extract the parison.
  connect SetInputConnection [rubberPlane GetOutputPort]
  parison SetInputConnection [connect GetOutputPort]
  normals SetInputConnection [parison GetOutputPort]
  normals SetFeatureAngle 0
}

proc plan_z {} {
  global VTKPlot
  global reader rubberPlane connect parison normals
  rubberPlane SetInput [reader GetOutput]
  rubberPlane SetScaleFactor 0.00
  # The threshold filter has been used to extract the parison.
  connect SetInputConnection [rubberPlane GetOutputPort]
  parison SetInputConnection [connect GetOutputPort]
  normals SetInputConnection [parison GetOutputPort]
  normals SetFeatureAngle 0
}

proc build_lut {} {
      global VTKPlot lut
      lut SetHueRange 0.6667 0.0
      lut Build 
}      


proc model_obj {} {
    global VTKPlot
    global parisonMapper normals lut modelActor ren
    parisonMapper SetInputConnection [normals GetOutputPort]
    parisonMapper SetLookupTable lut
    parisonMapper SetScalarRange $VTKPlot(min) $VTKPlot(max)
    modelActor SetMapper parisonMapper
    if {$VTKPlot(view)=="surface"} {
      [modelActor GetProperty] SetRepresentationToSurface
    }
    if {$VTKPlot(view)=="wireframe"} {
      [modelActor GetProperty] SetRepresentationToWireframe
    } 
    if {$VTKPlot(view)=="points"} { 
      [modelActor GetProperty] SetRepresentationToPoints
    }  
    ren AddActor modelActor
}
 
  
proc axes_obj {} {
  global VTKPlot
  global axes tprop1 tprop2 tprop3 marker1
  axes SetShaftTypeToCylinder
  axes SetXAxisLabelText "x"
  axes SetYAxisLabelText "y"
  axes SetZAxisLabelText "z"
  axes SetTotalLength 1.5 1.5 1.5

  tprop1 ItalicOn
  tprop1 ShadowOn
  tprop1 SetFontFamilyToTimes
  [ axes GetXAxisCaptionActor2D ] SetCaptionTextProperty tprop1
  tprop2 ShallowCopy tprop1
  [ axes GetYAxisCaptionActor2D ] SetCaptionTextProperty tprop2
  tprop3 ShallowCopy tprop1
  [ axes GetZAxisCaptionActor2D ] SetCaptionTextProperty tprop3
 
  marker1 SetOutlineColor 0.93 0.57 0.13
  marker1 SetOrientationMarker axes
  marker1 SetViewport 0.0 0.0 0.2 0.3
  marker1 SetInteractor [renWin GetInteractor]
  marker1 SetEnabled 1
  marker1 InteractiveOff
}



proc scalar_bar_obj {} {
  global VTKPlot
  global scalarBar parisonMapper ren
  scalarBar SetLookupTable [parisonMapper GetLookupTable]
  scalarBar SetTitle $VTKPlot(scalarname)
  [scalarBar GetPositionCoordinate] SetCoordinateSystemToNormalizedViewport
  [scalarBar GetPositionCoordinate] SetValue 0.86 0.1
  scalarBar SetOrientationToVertical
  scalarBar SetWidth 0.12
  scalarBar SetHeight 0.9
  scalarBar SetLabelFormat "%-#6.3e"
  ren AddActor2D scalarBar
}

proc model_axes_obj {} {
    global VTKPlot
    global modelAxes rubberPlane ren tprop1
    modelAxes SetInput [rubberPlane GetOutput]
    modelAxes SetCamera [ren GetActiveCamera]
    modelAxes SetScaling 1
    modelAxes SetZLabel "scaled value"
    modelAxes SetLabelFormat "%6.4g"
    modelAxes SetFlyModeToOuterEdges
    modelAxes SetFontFactor 0.8
    modelAxes SetAxisTitleTextProperty tprop1
    modelAxes SetAxisLabelTextProperty tprop1

    ren AddViewProp modelAxes
}    

proc draw_window {} {
  global VTKPlot
  global ren renWin
  ren ResetCamera
  #adjust angle of z and xy-plane 
  [ren GetActiveCamera] Azimuth 0
  #adjust angle roll from x
  [ren GetActiveCamera] Roll 10
  #resize rate
  [ren GetActiveCamera] Dolly 1
  #adjust Camera
  [ren GetActiveCamera] ComputeViewPlaneNormal
  [ren GetActiveCamera] SetViewUp 1 1 1
  [ren GetActiveCamera] OrthogonalizeViewUp
  ren ResetCameraClippingRange
  ren SetBackground 0.33 0.35 0.43
  renWin Render 
}

#---------- produce cmd line arg ----------------------------------------------
if {$argc} {
  global VTKPlot
  set VTKPlot(vtk_data_file) [lindex $argv [expr $argc - 2]]
  set VTKPlot(scalarname)    [lindex $argv [expr $argc - 1]]
  if {[string index $VTKPlot(vtk_data_file) 0] != "-" && [file exists $VTKPlot(vtk_data_file)]} {
    cd [file dirname  $VTKPlot(vtk_data_file)]
    read_data 
    plan_z
    build_lut
    model_obj
    model_axes_obj
    axes_obj
    scalar_bar_obj
    draw_window
  }
}

tkwait window .
