
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

package require Tk
package require ctext


proc error_exit {msg} {
  wm withdraw .
  tk_dialog .error Error $msg error 0 Exit
  exit 1
}

if [catch {do_config -win} msg] {error_exit $msg}


#---- set global data structure -----------------------------------------------
set exe_pid 0

array set ProgData {
  cmd_file_name ""
  change "0"
  pos    "1.0"
}

#----- create main window -----------------------------------------------------
wm title . "GSS UI 0.1"
wm minsize    .  800 600
wm protocol   . WM_DELETE_WINDOW save_quit


proc do_quit {} {
  exit 0
}

#----- create menu -------------------------------------------------------------
menubar_create {File Edit Run Help}

#----- file menu ---------------------------------------------------------------
set m [menubar_get File]
$m add command -label "New"      -command new_cmd_file
$m add command -label "Open ..." -command open_cmd_file
$m add separator
$m add command -label "Save ..." -command save_cmd_file
$m add command -label "Save As..." -command saveas_cmd_file
$m add separator
$m add command -label "Quit" -command save_quit


proc new_cmd_file {} {
  global ProgData
  global editor
  global output
  $editor delete 1.0 end
  $output delete 1.0 end
  set ProgData(cmd_file_name) ""
  $editor edit modified 0
  set ProgData(change) 0
}


proc open_cmd_file {} {
  global ProgData
  global editor
  $editor delete 1.0 end
  set FileList {
  {{GSS Input Files} {.inp .txt}}
  {{All Files} {*}}
  }

  set ProgData(cmd_file_name) [FileOpen "Open GSS Input File" "" . $FileList]

  if {![catch {open $ProgData(cmd_file_name) r} fd]} {
    $editor fastinsert end [read $fd]
    close $fd
  }
  $editor highlight 1.0 end
  $editor edit modified 0
  set ProgData(change) 0
  cd [file dirname  $ProgData(cmd_file_name)]
}

proc save_cmd_file {} {
  global ProgData
  global editor
  if {$ProgData(cmd_file_name) == ""} {
     set ProgData(cmd_file_name) [FileSave "Save GSS Input File" "" . \
         {{{GSS Input Files} {.inp}} {{All Files} {*}}} inp]
     if {$ProgData(cmd_file_name) == ""} return
  }
  set fd [open $ProgData(cmd_file_name) w+]
  puts $fd [string trimright [$editor get 1.0 end]]
  close $fd
  $editor edit modified 0
  set ProgData(change) 0
  cd [file dirname  $ProgData(cmd_file_name)]
}

proc saveas_cmd_file {} {
  global ProgData
  global editor
  set ProgData(cmd_file_name) [FileSave "Save As GSS Input File" "" . \
       {{{GSS Input Files} {.inp}} {{All Files} {*}}} inp]
  if {$ProgData(cmd_file_name) == ""} return
  set fd [open $ProgData(cmd_file_name) w+]
  puts $fd [string trimright [$editor get 1.0 end]]
  close $fd
  $editor edit modified 0
  set ProgData(change) 0
  cd [file dirname  $ProgData(cmd_file_name)]
}

proc save_quit {} {
 global ProgData
 global editor
 if {[$editor edit modified]} then {
     set save [dialog .quitsave -1 -1 Quit \
      "Save the command file?" question 0 Yes No Cancel]
     if {$save == 2} {return 0}
     if {$save == 0} {save_cmd_file}
 }
 do_quit
}


#----- edit menu ---------------------------------------------------------------
set m [menubar_get Edit]
$m add command -label "Undo" -underline 0 -accelerator {Ctrl-z}  -command edit_undo
$m add command -label "Redo" -underline 0 -accelerator {Ctrl-y}  -command edit_redo
$m add separator
$m add command -label "Cut" -underline 0 -accelerator {Ctrl-x}  -command edit_cut
$m add command -label "Copy" -underline 0 -accelerator {Ctrl-c} -command edit_copy
$m add command -label "Paste" -underline 0 -accelerator {Ctrl-v} -command edit_paste
$m add command -label "Select All" -underline 0  -command edit_select_all
#$m add separator
#$m add command -label "Find/Replace" -underline 0 -command edit_find

proc edit_undo {} {
  global editor
  global ProgData
  catch {$editor edit undo} msg
  #puts $msg
  set ProgData(change) [$editor edit modified]
  $editor highlight 1.0 end
}

proc edit_redo {} {
  global editor
  global ProgData
  catch {$editor edit redo} msg
  #puts $msg
  set ProgData(change) [$editor edit modified]
  $editor highlight 1.0 end
}

proc edit_cut {} {
  global editor
  global ProgData
  $editor cut
  set ProgData(change) [$editor edit modified]
  $editor highlight 1.0 end
}

proc edit_copy {} {
  global editor
  $editor copy
}

proc edit_paste {} {
  global editor
  global ProgData
  $editor paste
  set ProgData(change) [$editor edit modified]
  $editor highlight 1.0 end
}

proc edit_select_all {} {
  global editor
  $editor tag add sel 1.0 end
}


#----- run menu ---------------------------------------------------------------
set m [menubar_get Run]
$m add command -label "Start" -underline 0 -command exe_cmd_file
$m add command -label "Stop" -underline 0 -state disabled -command exe_kill 

proc tk_exec_fileevent {id} {
    global tk_exec_data
    global tk_exec_cond
    global tk_exec_pipe
    global output

    if {[eof $tk_exec_pipe($id)]} {
        fileevent $tk_exec_pipe($id) readable ""
        set tk_exec_cond($id) 1
        return
    }
    #append tk_exec_data($id) [gets $tk_exec_pipe($id)]
    #append tk_exec_data($id) "\n"
    set data [gets $tk_exec_pipe($id)]
    append data "\n"
    $output insert end $data
    $output see end
}

proc tk_exec {prog} {
    global tk_exec_id
    global tk_exec_data
    global tk_exec_cond
    global tk_exec_pipe
    global tcl_platform
    global env

    if {![info exists tk_exec_id]} {
        set tk_exec_id 0
    } else {
        incr tk_exec_id
    }

    set keepnewline 0

    for {set i 0} {$i < [llength $prog]} {incr i} {
        set arg [lindex $prog $i]
        switch -glob -- $arg {
            -keepnewline {
                set keepnewline 1
            }
            -- {
                incr i
                break
            }
            -* {
                error "unknown option: $arg"
            }
            ?* {
                # the glob should be on *, but the wiki reformats
                # that as a bullet
                break
            }
        }
    }

    if {$i > 0} {
        set prog [lrange $prog $i end]
    }

    if {$tcl_platform(platform) == "windows" && \
        [info exists env(COMSPEC)]} {
        set prog [linsert $prog 0 $env(COMSPEC) "/c"]
    }


    #set pipe [open "| $prog" r]
    if { [catch {open "| $prog 2>@stdout"} FILEHANDLE] } {
      return "Can't open pipe for '$args'"
    }
    set pipe $FILEHANDLE

    set tk_exec_pipe($tk_exec_id) $pipe
    set tk_exec_data($tk_exec_id) ""
    set tk_exec_cond($tk_exec_id) 0

    fconfigure $pipe -blocking 0 -buffering none

    fileevent $pipe readable "tk_exec_fileevent $tk_exec_id"

    vwait tk_exec_cond($tk_exec_id)

    if {$keepnewline} {
        set data $tk_exec_data($tk_exec_id)
    } else {
        set data [string trimright $tk_exec_data($tk_exec_id) \n]
    }

    unset tk_exec_pipe($tk_exec_id)
    unset tk_exec_data($tk_exec_id)
    unset tk_exec_cond($tk_exec_id)

    if {[catch {close $pipe} err]} {
        error "pipe error: $err"
    }
    return $data
}

proc bgExec {prog readHandler pCount {timeout 0} {toExit ""}} {
      global exe_pid
      upvar #0 $pCount myCount
      if {![string length [auto_execok [lindex $prog 0]]]} {
         # perhaps additional checking with [file executable]
         return -code error "error: could not locate '$prog'"
      }
      set myCount [expr {[info exists myCount]?[incr myCount]:1}]
      set redir [expr {[info patchlevel] >= "8.4.7"?{2>@1}:{2>@stdout}}]
      if [catch {open "| $prog $redir" r} pH] {
         return -code error "error: could not start '$prog' ($pH)"
      }
      set exe_pid [pid $pH]
      fconfigure $pH -blocking 0; # -buffering line (does it really matter?!)
      set tID [expr {$timeout?[after $timeout [list bgExecTimeout $pH $pCount $toExit]]:{}}]
      fileevent $pH readable [list bgExecGenericHandler $pH $pCount $readHandler $tID]
      return $pH
 }

 proc bgExecGenericHandler {chan pCount readHandler tID} {
      upvar #0 $pCount myCount
      if {[eof $chan]} {
         after cancel $tID;   # empty tID is ignored
         catch {close $chan}; # automatically deregisters the fileevent handler
                              # (see Practical Programming in Tcl an Tk, page 229)
         incr myCount -1
      } elseif {[gets $chan line] != -1} {
         # we are not blocked (manpage gets, Practical... page.233)
         lappend readHandler $line
         if {[catch {uplevel $readHandler} rc]} {
            # user-readHandler ended with error -> terminate the processing
            after cancel $tID
            catch {close $chan}
            incr myCount -1
         }
      }
 }

 proc bgExecTimeout {chan pCount toExit} {
      upvar #0 $pCount myCount
      if {[string length $toExit]} {
         catch {uplevel [list $toExit [pid $chan]]}
      }
      catch {close $chan}
      incr myCount -1
 }

proc exe_cmd_record {data}  {
  global output
  append data "\n"
  $output insert end $data
  $output see end
}

proc exe_cmd_file {} {
  global ProgData
  global output
  global editor
  global tcl_platform
  global .toolbar.but.kill
  
  .toolbar.but.run.kill configure -state normal

  foreach i {2 } {
    menubar_state Run normal $i
  }
  if {[$editor edit modified]==1} then {
     set save [dialog .runsave -1 -1 Quit \
      "Save the command file?" question 0 Yes No Cancel]
     if {$save == 2} {return 0}
     if {$save == 0} {save_cmd_file}
  }
  if {$tcl_platform(platform) == "windows" } {
     set cmd [file nativename "$::env(GSS_DIR)\\bin\\gss"]
     set arg [file tail $ProgData(cmd_file_name)]
  } else {
     set cmd "$::env(GSS_DIR)/bin/gss"
     set arg [file tail $ProgData(cmd_file_name)]
  }
  $output delete 1.0 end
  bgExec "gss $arg " exe_cmd_record pCount

# if {$tcl_platform(platform) == "windows" } {
#       tk_exec "gss $ProgData(cmd_file_name)"
#  } else {
#       bgExec "gss $ProgData(cmd_file_name) " exe_cmd_record pCount
#  }

}

proc exe_kill {} {
  global tcl_platform
  global exe_pid
  global output
  exec kill -s SIGINT $exe_pid
  $output insert end "Terminated by signal SIGINT \n"
  # disable kill buttom
  .toolbar.but.run.kill configure -state disabled
  foreach i {2 } {
    menubar_state Run disabled $i
  }
}



#----- help menu --------------------------------------------------------------
set m [menubar_get Help]
$m add command -label "About..." -underline 0 -command do_about

proc do_about {} {
  global ProgData
  dialog .about -1 -1 "About GSS GUI" \
 {            GSS Editor 0.1
 A GUI shell for GSS code.
 Date   : 10/27/2006
 Author : Gong Ding
 Email  : gdiso@ustc.edu} \
 img_about 0 Close

}


#---------- toolbar ------------------------------------------------------------
frame .toolbar
pack .toolbar -side top -pady 2 -fill x

set tbbut [frame .toolbar.but]
pack $tbbut -side left

set b [frame $tbbut.file]
pack $b -side left -padx 5

image create photo img_new -data {\
R0lGODlhEAAQAIMAAPwCBFRShPz6/PT2/Oz2/OTu/MzK/OTy/Nzu/NTq/NTm/Mzm/Mzi/MTi/Lze/Lza/C\
H5BAEAAAAALAAAAAAQABAAAARcEMhJKw04YyuDGCBBbFYQikVgBKV4HEUqrFVwEDCCaHVxIAVEQqGgXYK6\
hHJhnAR0iuFiwWh2hokptWEFBBRaRqPh6Aa044ejzLqIG4+x42E22O/3tlPD13P+fhEAIf5oQ3JlYXRlZC\
BieSBCTVBUb0dJRiBQcm8gdmVyc2lvbiAyLjUNCqkgRGV2ZWxDb3IgMTk5NywxOTk4LiBBbGwgcmlnaHRz\
IHJlc2VydmVkLg0KaHR0cDovL3d3dy5kZXZlbGNvci5jb20AOw==}

image create photo img_open -data {\
R0lGODlhEAAQAIMAAPwCBASCBMyaBPzynPz6nJxmBPzunPz2nPz+nPzSBPzqnPzmnPzinPzenAAAAAAAAC\
H5BAEAAAAALAAAAAAQABAAAARTEMhJq724hp1n8MDXeaJgYtsnDANhvkJRCcZxEEiOJDIlKLWDbtebCBaG\
GmwZEzCQKxxCSgQ4Gb/BbciTCBpOoFbX9X6fChYhUZYU3vB4cXTxRwAAIf5oQ3JlYXRlZCBieSBCTVBUb0\
dJRiBQcm8gdmVyc2lvbiAyLjUNCqkgRGV2ZWxDb3IgMTk5NywxOTk4LiBBbGwgcmlnaHRzIHJlc2VydmVk\
Lg0KaHR0cDovL3d3dy5kZXZlbGNvci5jb20AOw==}

image create photo img_save -data {\
R0lGODlhEAAQAIQAAPwCBFRShGRmzPT6/Oz2/OTy/Nzu/NTm/AQCBMTi/Lze/Mzm/KzW/NTq/LTa/AQCxO\
Ti5Nze3Nza3MzOzLy+xNTS1MTGxMzKzAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQ\
ABAAAAVkICCOZCkGaKqqpxAMRFEYRx0gQfvGRZ0oKZ2MdlgogC6doVc8MgJJADRQaxidU13t8HMwnlGoa4\
UShM2PtFp9fkAikomcQnm0IebKxGKp3/MTF3R2OVICboB7foVSZGQmkCR+IQAh/mhDcmVhdGVkIGJ5IEJN\
UFRvR0lGIFBybyB2ZXJzaW9uIDIuNQ0KqSBEZXZlbENvciAxOTk3LDE5OTguIEFsbCByaWdodHMgcmVzZX\
J2ZWQuDQpodHRwOi8vd3d3LmRldmVsY29yLmNvbQA7}



button $b.new -relief flat -image img_new -takefocus 0 -command new_cmd_file
set_balloon $b.new "New File..."

button $b.open -relief flat -image img_open -takefocus 0 -command open_cmd_file
set_balloon $b.open "Open File..."

button $b.save -relief flat -image img_save -takefocus 0 -command save_cmd_file
set_balloon $b.save "Save File ..."

pack $b.new $b.open $b.save -side left


set b [frame $tbbut.run]
pack $b -side left -padx 5

image create photo img_run -data {\
R0lGODlhEAAQAIEAAAT+BPz+/ASCBAAAACH5BAEAAAAALAAAAAAQABAAAAIrRI6ZFuIPHYr0HXpjhpM3vAmh94m\
fokgn2o3W5pqlCc60Zt8VnoMMqwD4CwAh/mhDcmVhdGVkIGJ5IEJNUFRvR0lGIFBybyB2ZXJzaW9uIDIuNQ0KqS\
BEZXZlbENvciAxOTk3LDE5OTguIEFsbCByaWdodHMgcmVzZXJ2ZWQuDQpodHRwOi8vd3d3LmRldmVsY29yLmNvbQA7}

image create photo img_stop -data {\
R0lGODlhEAAQAIUAAASC/PR6fPyChPyGjPyOjPyKjPxydPy6vPzGzPzCxPy2vPymrLxOVPyurPyanPyOlPx6\
fNRmZPwSFPyepPx2fPyipPx+hPxaXPxSVLQyNPw2NPySlPz+/Pz2/Pze3PxCROQiJEQCBPweJPyytPzm7Pz6/\
Pw2PPwuNMQGDEwCBPwaHPwmJPwOFPwmLPwCBPze5PwKDPzy9PyenPwGBKwCBPQCBLQCBFQGBJQGBLwKDLwGDA\
AAAAAAAAAAAAAAAAAAACH5BAEAAAAALAAAAAAQABAAAAabQIBwSCwWA4IBYVAwGIWCwwGBSCgWDKOi4SgUHgV\
IwBAhSiYCSkUwLlgumMxQU7hsOAiMoePBfEAhQiJ+GAkcIyQlFCYnKClCKissLSsTHB0YKhIuKIEALJsuKi8d\
JQMuqJ1CMDAuGiQxGgIcMjMuNJ41qCsJH6gyBbaqQjSouqjIjkU2yM0oN084OTo6KDQ20E8AIdzd2t9EfkEAI\
f5oQ3JlYXRlZCBieSBCTVBUb0dJRiBQcm8gdmVyc2lvbiAyLjUNCqkgRGV2ZWxDb3IgMTk5NywxOTk4LiBBbGw\
gcmlnaHRzIHJlc2VydmVkLg0KaHR0cDovL3d3dy5kZXZlbGNvci5jb20AOw==}

button $b.run -relief flat -image img_run -takefocus 0 -command exe_cmd_file
set_balloon $b.run "Run GSS"

button $b.kill -relief flat -image img_stop -takefocus 0 -command exe_kill -state disabled
set_balloon $b.kill "Kill process"

pack $b.run $b.kill -side left


#---------- color text editor -------------------------------------------------
pack [frame .f_ctex] -fill both -expand 1
pack [scrollbar .f_ctex.s -command {.f_ctex.t yview}] -side right -fill y
set editor [ctext .f_ctex.t -bg white -fg black -insertbackground red  -undo 1 \
      -yscrollcommand {.f_ctex.s set}]
pack  $editor -pady 2 -fill both -expand 1


#ctext::addHighlightClassForRegexp .f_ctex.t string black {=[ \t]*([a-zA-Z0-9\_][a-zA-Z0-9\.\_]*|\"[a-zA-Z0-9\.\_\ \-]*\")}
ctext::addHighlightClassForRegexp .f_ctex.t number red {=[ \t]*([+-]?\d+\.\d*([EeDd][+-]?\d+)?|[+-]?\d*\.\d+([EeDd][+-]?\d+)?|[+-]?\d+([EeDd][+-]?\d+)|[+-]?\d+)}
ctext::addHighlightClassForRegexp .f_ctex.t bool DeepPink {=[ \t]*(TRUE|FALSE|On|Off)}
ctext::addHighlightClassForRegexp .f_ctex.t equ black {\=}
ctext::addHighlightClassForRegexp .f_ctex.t keyword blue \
      {^[ \t]*(SET|MESH|XMESH|YMESH|ELIMINATE|SPREAD|REGION|SEGMENT|REFINE|PROFILE|\
               ISOURCE|VSOURCE|BOUNDARY|CONTACT|PMIS|ATTACH|METHOD|SOLVE|IMPORT|EXPORT|\
               PLOTVTK|PLOT|PLOTMESH|PROBE|PHOTOGEN|END)}
ctext::addHighlightClassForRegexp .f_ctex.t comment  #007f7f {#[^\n\r]*}

bindtags $editor "Ctext $editor all ."

bind $editor <KeyRelease>        key_release
bind $editor <ButtonRelease>     button_release

bind Text <<Cut>>   {}
bind Text <<Copy>>  {}
bind Text <<Paste>> {}
bind Text <<Undo>>   {}
bind Text <<Redo>>   {}

bind $editor <Control-x>   {edit_cut}
bind $editor <Control-c>   {edit_copy}
bind $editor <Control-v>   {edit_paste}
bind $editor <Control-z>   {edit_undo}
bind $editor <Control-y>   {edit_redo}

proc button_release {} {
  global ProgData
  global editor
  set ProgData(pos) [$editor index insert]
}

proc key_release {} {
  global ProgData
  global editor
  set ProgData(change) [$editor edit modified]
  set ProgData(pos) [$editor index insert]
}

bind $editor <Button-3> {tk_popup [menubar_get Edit] [winfo pointerx . ] [winfo pointery .]}

#---------- text window for log -----------------------------------------------
FrameCreate .f_out -text "GSS Log Message" -font $Font(bold) -pady 0
pack .f_out  -pady 2 -fill both -expand 1
set outwin [FrameGet .f_out]

pack [frame $outwin.f] -fill both -expand 1
pack [scrollbar $outwin.f.s -command {$outwin.f.t yview}] -side right -fill y
set output [text $outwin.f.t -bg white -fg black -height 4  \
      -yscrollcommand {$outwin.f.s set}]
pack  $output -pady 2 -fill both -expand 1



#---------- status bar ----------------------------------------------
frame .status -borderwidth 1 -relief flat
pack  .status  -padx 1 -side bottom -fill x
label .status.l1 -text ""  -relief sunken -borderwidth 1 -width 10
label .status.l2 -text "" -relief sunken -borderwidth 1
label .status.filler -text "" -relief sunken -borderwidth 1
label .status.l3 -text "" -relief sunken -borderwidth 1 -width 13

canvas .status.c -relief sunken -borderwidth 0 -width 20 -height 16 -cursor bottom_right_corner

pack .status.l1 .status.l3 .status.l2 -side left -padx 1
pack .status.filler -side left -fill x -expand yes -padx 1
pack .status.c  -padx 1 -side left


.status.c create line 19 0   1 18 -fill white
.status.c create line 19 1   2 18 -fill gray70
.status.c create line 19 2   3 18 -fill gray70
.status.c create line 19 6   7 18 -fill white
.status.c create line 19 7   8 18 -fill gray70
.status.c create line 19 8   9 18 -fill gray70
.status.c create line 19 12   13 18 -fill white
.status.c create line 19 13   14 18 -fill gray70

proc status_change {args} {
   global ProgData
   global output
   global editor

  .status.l2 configure -text " Filename: $ProgData(cmd_file_name) "

   if {$ProgData(change) == 0} {
     .status.l1 configure -foreground blue -text  "Unchanged"
   } else {
     .status.l1 configure -foreground red -text "Changed"
   }

   .status.l3 configure -text "Line: $ProgData(pos)"
}

trace variable ProgData w status_change

#---------- produce cmd line arg ----------------------------------------------
if {$argc} {
  set file [lindex $argv [expr $argc - 1]]
  if {[string index $file 0] != "-" && [file exists $file]} {
    set ProgData(cmd_file_name) $file
    if {![catch {open $ProgData(cmd_file_name) r} fd]} {
      $editor fastinsert end [read $fd]
      close $fd
    }
    $editor highlight 1.0 end
    $editor edit modified 0
    set ProgData(change) 0
    cd [file dirname  $ProgData(cmd_file_name)]
  }
}

grab .
focus $editor
status_change

#---------- images ------------------------------------------------------------
image create photo img_about -file $env(GSS_DIR)/lib/gss.logo.gif
