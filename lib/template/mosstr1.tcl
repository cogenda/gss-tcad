#-------------------
#  mosstr1
#-------------------
# Template File for MOS Structures

puts "Generate GSS input file..."

#COMMENT   Derived Quantities and Definitions
#    WD       = depletion width under the drain (microns)
#    SUBDEP   = substrate depth (microns)
#    LMID     = distance to middle of gate 
#               from source edge (microns)
#    LMAX     = width of device (microns)
#    EPS      = roundoff error
set   WD        [expr 3.6E7*sqrt((.8+abs($VDBMAX))/$NSUB)]
set   SUBDEP    [expr int(2*($SDJUNC+$WD+1.5))/2]
set   LMID      [expr $LSOURCE+$LGATE/2]
set   LMAX      [expr $LSOURCE+$LGATE+$LDRAIN]
set   EPS       [expr 0.0002]

#COMMENT   Only consider LDD implants if LDDPEAK >= NSUB.
#         LDT ==> LDD implant
#         LDF ==> No LDD implant
if {$LDDPEAK > $NSUB} {
  set  LDT  1
  set  LDF  0
} else {
  set  LDT  0
  set  LDF  1
}

#COMMENT   Transistor type.
if {$TRANTYPE == "PMOS"} {
  set NMOS  0
  set PMOS  1
} else {
  set NMOS  1
  set PMOS  0
}

#COMMENT   Threshold implant impurity type.
if {$VTTYPE == "N"} {
  set VTN  1
  set VTP  0
} else {
  set VTN  0
  set VTP  1
}


#COMMENT   More Derived Quantities
#    SDDIF    = source/drain out-diffusion distance (microns)
#    LDDIF    = lightly doped drain 
#               out-diffusion distance (microns)

set NSRF  [expr $NSUB+$VTPEAK]
if { $NMOS && $VTN } {
  set NSRF  [expr abs($NSUB-$VTPEAK)]
}  
if { $PMOS && $VTP } {
 set NSRF   [expr abs($NSUB-$VTPEAK)]
} 

set NSDJ  [expr $NSUB]
set TMP   [expr $VTPEAK*exp(- pow($SDJUNC/$VTCHAR,2))]
if { $NMOS && $VTP } {
  set NSDJ  [expr $NSUB+$TMP]
}  
if { $PMOS && $VTN } {
 set NSDJ   [expr $NSUB+$TMP]
} 
set SDYC  [expr $SDJUNC/sqrt(log($SDPEAK/$NSDJ))]
set SDXC  [expr $SDYC*$LATD]
set SDDIF [expr $SDXC*sqrt(log($SDPEAK/$NSRF))]

if { $LDT } {
  set NLDJ  [expr $NSUB]
  set TMP   [expr $VTPEAK*exp(- pow($LDDJUNC/$VTCHAR,2))]
  if { $NMOS && $VTP } {
   set NLDJ  [expr $NSUB+$TMP]
  } 
  if { $PMOS && $VTN } {
   set NLDJ  [expr $NSUB+$TMP]
  } 
  set LDYC  [expr $LDDJUNC/sqrt(log($LDDPEAK/$NLDJ))]
  set LDXC  [expr $LDYC*$LATD]
  set LDDIF [expr $LDXC*sqrt(log($LDDPEAK/$NSRF))]
}

if { $LDF } {
 set LDDIF [expr 0.0 ]
}

#COMMENT   More definitions
#    LSSDIF   = distance from source edge 
#               to edge of source diffusion
#    LDDDIF   = distance from drain edge 
#		to edge of drain diffusion
#    SCE2JD   = source contact edge 
#		to junction distance (microns)
#    DCE2JD   = drain contact edge 
#		to junction distance (microns)
#    SCESP    = grid spacing at source contact edge (microns)
#    DCESP    = grid spacing at drain contact edge (microns)
#    LSLDIF   = distance from source edge 
#		to edge of LDD diffusion
#    LDLDIF   = distance from drain edge 
#		to edge of LDD diffusion

set LSSDIF    [expr $LSOURCE-$LSPACER+$SDDIF]
set LDDDIF    [expr $LDRAIN-$LSPACER+$SDDIF]
set SCE2JD    [expr $LSSDIF-$LSCONT]
set DCE2JD    [expr $LDDDIF-$LDCONT]
set SCESP     [expr $SCE2JD-($SCE2JD-$JUNCSP)/$RATIO]
set DCESP     [expr $DCE2JD-($DCE2JD-$JUNCSP)/$RATIO]
set LSLDIF    [expr $LSOURCE+$LDDIF]
set LDLDIF    [expr $LDRAIN+$LDDIF]
 


#COMMENT   Evaluate Statement Mask Logicals
#    A  ==> (SDDIF+JUNCSP) <  LSPACER
#    A1 ==> (SDDIF+JUNCSP) >= LSPACER
#    B  ==> A1 & (SDDIF < (LSPACER+JUNCSP))
#    B1 ==> B & LDT & (2*JUNCSP <  LDDIF)
#    B2 ==> B & LDT & (2*JUNCSP >= LDDIF)
#    C  ==> A1 & (SDDIF >= (LSPACER+JUNCSP))
#    C1 ==> C & LDT & ((SDDIF+2*JUNCSP) <  (LDDIF+LSPACER))
#    C2 ==> C & LDT & ((SDDIF+2*JUNCSP) >= (LDDIF+LSPACER))
#    D  ==> LDT & (SDDIF <  (LSPACER+LDDIF))
#    D1 ==> LDT & (SDDIF >= (LSPACER+LDDIF))
#    E  ==> D & (LDDJUNC <  SDJUNC)
#    E1 ==> D & (LDDJUNC >= SDJUNC)

set A   0
set A1  0
set B   0
set B1  0
set B2  0
set C   0
set C1  0
set C2  0
set D   0
set D1  0
set E   0
set E1  0
set F   0
set F1  0
set G   0
set G1  0
set G2  0
set G3  0

set TMP [expr int(($SDDIF+$JUNCSP+$EPS)/($LSPACER+$EPS))]

if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
  set A   0 
  set A1  1 
} else {
  set A   1
  set A1  0 
}

if { $A &&  $LDT} {
set WD   [expr ($LSPACER-$SDDIF+$LDDIF)/2]
set SPC  [expr $WD-($WD-$JUNCSP)/$RATIO]
set SPT  [expr $LSSDIF+$WD]
set DPT  [expr $LDDDIF+$WD]
}

#COMMENT   More masks
#    G  ==> A & LDT & ((SDDIF+WD+SPC) <  LSPACER)
#    G1 ==> A & LDT & ((SDDIF+WD+SPC) >= LSPACER)
#    G2 ==> G1 & ((SDDIF+WD) <  (LSPACER+SPC))
#    G3 ==> G1 & ((SDDIF+WD) >= (LSPACER+SPC))

if { $A1 } { 
 set TMP  [expr int($SDDIF/($LSPACER+$JUNCSP))]
 if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
   set B   0 
   set C   1 
 } else {
   set B   1
   set C   0 
 }
}

if { $B && $LDT} {
    set TMP  [expr int((2*$JUNCSP)/$LDDIF)]
    if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set B1   0 
     set B2   1 
    } else {
     set B1   1
     set B2   0 
 }
}

if { $C && $LDT} {
    set TMP  [expr int(($SDDIF+2*$JUNCSP)/($LDDIF+$LSPACER))
    if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set C1   0 
     set C2   1 
    } else {
     set C1   1
     set C2   0 
    }
}

if { $LDT } {
  set TMP  [expr int($SDDIF/($LSPACER+$LDDIF))]
  if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set D   0 
     set D1  1 
  } else {
     set D   1
     set D1  0 
  }
}

if { $D }  {
  set TMP  [expr int($LDDJUNC/$SDJUNC) ]
  if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set E   0 
     set E1  1 
  } else {
     set E   1
     set E1  0 
  }
}

if { $A && $LDT} {
  set   TMP      [expr int(($SDDIF+$WD+$SPC)/$LSPACER) ]
  if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set G   0 
     set G1  1 
  } else {
     set G   1
     set G1  0 
  }
}

if { $G1 } {
  set TMP  [expr int(($SDDIF+$WD)/($LSPACER+$SPC))]
  if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set G2   0 
     set G3  1 
  } else {
     set G2   1
     set G3  0 
  }
}

set FileID [open $FileName w]
puts "Write GSS input statements to $FileName"

#-----------------------
#  MESH statement
#-----------------------
puts $FileID "MESH     Type=GSS   ModelFile=$TRANTYPE.cgns Triangle=pzAY"

#-----------------------
#  XMESH statements
#-----------------------
puts $FileID  "#---------------------------------------------------"
set  WD  [expr $SCESP]
set  TMP [expr int($SCESP/$LSCONT) ]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X1   1 
} else {
     set X1   0 
}

if { $X1 } {
 set WD  $LSCONT
}

puts $FileID  "XMESH    X.MAX=$LSCONT   H2=$WD      RATIO=$RATIO"

if {$A} {
 puts $FileID "XMESH    X.MAX=$LSSDIF   H2=$JUNCSP  RATIO=$RATIO"
} 

if {$G} {
 puts $FileID  "XMESH    X.MAX=$SPT      H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID  "XMESH    X.MAX=$LSOURCE  H1=$SPC     RATIO=1/$RATIO"
 puts $FileID  "XMESH    X.MAX=$LSLDIF   H2=$JUNCSP  RATIO=$RATIO"
}

if {$G2} {
 puts $FileID  "XMESH    X.MAX=$LSOURCE  H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID  "XMESH    X.MAX=$LSLDIF   H2=$JUNCSP  RATIO=$RATIO"
}

if {$G3} {
 puts $FileID  "XMESH    X.MAX=$LSOURCE  H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID  "XMESH    X.MAX=$SPT      H2=$SPC     RATIO=1/$RATIO"
 puts $FileID  "XMESH    X.MAX=$LSLDIF   H2=$JUNCSP  RATIO=$RATIO"
}

if {$A && $LDT} {
 puts $FileID  "XMESH    X.MAX=$LMID     H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID  "XMESH    X.MAX=$LMAX-$LDLDIF  H2=$JUNCSP  RATIO=$RATIO"
}

if {$G3} {
 puts $FileID  "XMESH    X.MAX=$LMAX-$DPT     H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID  "XMESH    X.MAX=$LMAX-$LDRAIN  H1=$SPC     RATIO=1/$RATIO"
 puts $FileID  "XMESH    X.MAX=$LMAX-$LDDDIF  H2=$JUNCSP  RATIO=$RATIO"
}

if {$G2} {
 puts $FileID "XMESH    X.MAX=$LMAX-$LDRAIN  H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID "XMESH    X.MAX=$LMAX-$LDDDIF  H2=$JUNCSP  RATIO=$RATIO"
}

if {$G} {
 puts $FileID "XMESH    X.MAX=$LMAX-$LDRAIN  H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID "XMESH    X.MAX=$LMAX-$DPT     H2=$SPC     RATIO=1/$RATIO"
 puts $FileID "XMESH    X.MAX=$LMAX-$LDDDIF  H2=$JUNCSP  RATIO=$RATIO"
}

if {$A && $LDF } {
 set TMP    [expr $LSPACER-$SDDIF]
 set GCESP  [expr $TMP-($TMP-$JUNCSP)/$RATIO]
 puts $FileID "XMESH    X.MAX=$LSOURCE  H1=$JUNCSP  RATIO=$RATIO"
 puts $FileID "XMESH    X.MAX=$LMID     H1=$GCESP   RATIO=$RATIO"
 puts $FileID "XMESH    X.MAX=$LMAX-$LDRAIN  H2=$GCESP   RATIO=$RATIO"
 puts $FileID "XMESH    X.MAX=$LMAX-$LDDDIF  H2=$JUNCSP  RATIO=$RATIO"
}

if {$A}  {  puts $FileID  "XMESH    X.MAX=$LMAX-$LDCONT  H1=$JUNCSP  RATIO=$RATIO" }

if {$B}  {  puts $FileID  "XMESH    X.MAX=$LSOURCE  H2=$JUNCSP  RATIO=$RATIO"}
if {$B1} {  puts $FileID  "XMESH    WIDTH=$LDDIF/2  H1=$JUNCSP  RATIO=$RATIO"}
if {$B1} {  puts $FileID  "XMESH    WIDTH=$LDDIF/2  H2=$JUNCSP  RATIO=$RATIO"}
if {$B}  {  puts $FileID  "XMESH    X.MAX=$LMID     H1=$JUNCSP  RATIO=$RATIO"}
if {$B1} {  puts $FileID  "XMESH    X.MAX=$LMAX-$LDLDIF  H2=$JUNCSP  RATIO=$RATIO"}
if {$B1} {  puts $FileID  "XMESH    WIDTH=$LDDIF/2  H1=$JUNCSP  RATIO=$RATIO"}
if {$B1} {  puts $FileID  "XMESH    WIDTH=$LDDIF/2  H2=$JUNCSP  RATIO=$RATIO"}
if {$B2} {  puts $FileID  "XMESH    X.MAX=$LMAX-$LDRAIN  H2=$JUNCSP  RATIO=$RATIO"}
if {$B && $LDF} {  puts $FileID  "XMESH    X.MAX=$LMAX-$LDRAIN  H2=$JUNCSP  RATIO=$RATIO"}
if {$B}  {  puts $FileID  "XMESH    X.MAX=$LMAX-$LDCONT  H1=$JUNCSP  RATIO=$RATIO"}

if {$C} {
   set TMP    [expr $SDDIF-$LSPACER]
   set GCESP  [expr $TMP-($TMP-$JUNCSP)/$RATIO]
   puts $FileID  "XMESH    X.MAX=$LSOURCE  H1=$SCESP   H2=$GCESP"
   puts $FileID  "XMESH    X.MAX=$LSSDIF   H2=$JUNCSP  RATIO=$RATIO"
}
if {$C1} {puts $FileID  "XMESH    WIDTH=[expr ($LDDIF-$SDDIF+$LSPACER)/2]   H1=$JUNCSP  RATIO=$RATIO"}
if {$C1} {puts $FileID  "XMESH    WIDTH=[expr ($LDDIF-$SDDIF+$LSPACER)/2]   H2=$JUNCSP  RATIO=$RATIO"}
if {$C}  {puts $FileID  "XMESH    X.MAX=$LMID     H1=$JUNCSP  RATIO=$RATIO"}
if {$C1} {puts $FileID  "XMESH    X.MAX=$LMAX-$LDLDIF  H2=$JUNCSP  RATIO=$RATIO"}
if {$C1} {puts $FileID  "XMESH    WIDTH=[expr ($LDDIF-$SDDIF+$LSPACER)/2]  H1=$JUNCSP   RATIO=$RATIO"}
if {$C1} {puts $FileID  "XMESH    WIDTH=[expr ($LDDIF-$SDDIF+$LSPACER)/2]  H2=$JUNCSP   RATIO=$RATIO"}
if {$C2}        {puts $FileID  "XMESH    X.MAX=$LMAX-$LDDDIF  H2=$JUNCSP  RATIO=$RATIO"}
if {$C && $LDF} {puts $FileID  "XMESH    X.MAX=$LMAX-$LDDDIF  H2=$JUNCSP  RATIO=$RATIO"}
if {$C}         {puts $FileID  "XMESH    X.MAX=$LMAX-$LDRAIN  H1=$JUNCSP  RATIO=$RATIO"}
if {$C}         {puts $FileID  "XMESH    X.MAX=$LMAX-$LDCONT  H1=$GCESP   H2=$DCESP"}

set WD  [expr $DCESP]
set TMP [expr int($DCESP/$LDCONT)]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X1   1 
} else {
     set X1   0 
}
if {$X1} {set WD  [expr $LDCONT]}

puts $FileID  "XMESH    X.MAX=$LMAX  H1=$WD  RATIO=$RATIO"

#-----------------------
#  Y.MESH statements
#-----------------------
#COMMENT  JUNC1 = deepest vertical junction location
#         JUNC2 = shallowest vertical junction location 
#                (if it exists)
puts $FileID  "#---------------------------------------------------"
set  JUNC1  [expr $SDJUNC]
set  JUNC2  [expr $SDJUNC]
if {$E}   {set JUNC2  [expr $LDDJUNC]}
if {$E1}  {set JUNC1  [expr $LDDJUNC]}
if {$E1}  {set JUNC2  [expr $SDJUNC]}

#COMMENT   F  ==> D & ((JUNC2+2*JUNCSP) <  JUNC1)
#+         F1 ==> D & ((JUNC2+2*JUNCSP) >= JUNC1)

if {$D} {
  set TMP  [expr int(($JUNC2+2*$JUNCSP)/$JUNC1)]
  if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set F   0
     set F1  1 
  } else {
     set F   1
     set F1  0
  }
}

puts $FileID  "YMESH    Y.TOP=0.1+$TOX  Y.BOTTOM=$TOX  N.SPACES=1"
puts $FileID  "YMESH    Y.BOTTOM=0  N.SPACES=1"

if {$F}   {puts $FileID  "YMESH    Y.BOTTOM=-$JUNC2  H1=$CHANSP  H2=$JUNCSP  #RATIO=$RATIO"}
if {$F}   {puts $FileID  "YMESH    Y.BOTTOM=-$JUNC1  H1=$JUNCSP  H2=$JUNCSP  #RATIO=$RATIO"}
if {$F1}  {puts $FileID  "YMESH    Y.BOTTOM=-$JUNC1  H1=$CHANSP  H2=$JUNCSP  #RATIO=$RATIO"}
if {$D1}  {puts $FileID  "YMESH    Y.BOTTOM=-$JUNC1  H1=$CHANSP  H2=$JUNCSP  #RATIO=$RATIO"}
if {$LDF} {puts $FileID  "YMESH    Y.BOTTOM=-$JUNC1  H1=$CHANSP  H2=$JUNCSP  #RATIO=$RATIO"}

puts $FileID  "YMESH    Y.BOTTOM=-$SUBDEP  H1=$JUNCSP  RATIO=$RATIO"
puts $FileID  "YMESH    Y.BOTTOM=-$SUBDEP-0.1  N.SPACES=1"

#-----------------------
#  ELIMINATE ROWS
#-----------------------
puts $FileID  "#---------------------------------------------------"

set WD   [expr (($RATIO+$EPS)*$JUNCSP-$CHANSP)/($RATIO+$EPS-1)]
set TMP  [expr int(2*$WD/$JUNC2)]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X1  1
  } else {
     set X1  0
}
  

if {$X1} {set WD  [expr $JUNC2/2]}

puts $FileID  "ELIMINATE  Direction=ROWS Y.TOP=-$EPS Y.BOTTOM=-$WD X.MIN=0  X.MAX=$LSOURCE-$EPS"
puts $FileID  "ELIMINATE  Direction=ROWS Y.TOP=-$EPS Y.BOTTOM=-$WD X.MIN=[expr $LMAX-$LDRAIN+$EPS]  X.MAX=$LMAX"

set Z1  [expr 1.5*$JUNCSP]
set Z2  [expr $JUNCSP*(2*($RATIO+$EPS)-1)/($RATIO+$EPS-1)]

set XMN  [expr $LSSDIF+$Z2]
set XMX  [expr $LMAX-$LDDDIF-$Z2]
set YMN  [expr $SDJUNC-$Z1]
set YMX  [expr $SDJUNC+$Z1]

if {$E1} {
 set XMN  [expr $LSLDIF+$Z2]
 set XMX  [expr $LMAX-$LDLDIF-$Z2]
 set YMN  [expr $LDDJUNC-$Z1]
 set YMX  [expr $LDDJUNC+$Z1]
}

#COMMENT   Make sure we can really do this eliminate.
set TMP  [expr int($XMN/$LMID)]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X1  0
  } else {
     set X1  1
}

if {$X1}  { puts $FileID  "ELIMINATE  Direction=ROWS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN  Y.BOTTOM=-$YMX" }

#-----------------------
#  ELIMINATE COLUMNS
#-----------------------
#COMMENT   The program treats the B2 and C2 cases as if there were
#	  only one lateral junction, whether there is an LDD 
#         implant or not. Furthermore, the B2 case puts the 
#         finest lateral spacing at LSOURCE and LMAX-LDRAIN.

set DSSDIF  [expr $LSSDIF]
set DDDDIF  [expr $LDDDIF]
set DSLDIF  [expr $LSLDIF]
set DDLDIF  [expr $LDLDIF]

if {$B2} {
 set DSSDIF  [expr $LSOURCE]
 set DDDDIF  [expr $LDRAIN]
 set DSLDIF  [expr $DSSDIF]
 set DDLDIF  [expr $DDDDIF]
}

if {$C2} {
 set DSLDIF  [expr $DSSDIF]
 set DDLDIF  [expr $DDDDIF]
}

set Z3  [expr 2*$JUNCSP*(1+$RATIO)]
set Z4  [expr $Z2+$JUNCSP*(1+$RATIO+$RATIO*$RATIO)]

set H   1
set H1  0
if {$D} {
 set TMP  [expr int(($DSSDIF+$Z3)/($DSLDIF))]
 if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set H   0
     set H1  1
 } else {
     set H   1
     set H1  0
 }
}

if {$H} {
 set XMN  [expr $DSSDIF-$Z1]
 set XMX  [expr $DSSDIF+$Z1]
 set YMN  [expr $SDJUNC+$Z2]
 if {$E1} {
   set XMN  [expr $DSLDIF-$Z1]
   set XMX  [expr $DSLDIF+$Z1]
   set YMN  [expr $LDDJUNC+$Z2]
 }
 puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"

 set XMN  [expr $LMAX-$DDDDIF-$Z1]
 set XMX  [expr $LMAX-$DDDDIF+$Z1]
 if {$E1} {
  set XMN  [expr $LMAX-$DDLDIF-$Z1]
  set XMX  [expr $LMAX-$DDLDIF+$Z1]
 }
 puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"

 if {$D} { set YMN  [expr $LDDJUNC+$Z2]}

 set X1  0
 if {$E} {
  set TMP  [expr int(($DSLDIF)/($DSSDIF+$Z4))]
  if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X1   0
  } else {
     set X1   1
  }
 }
 
 if {$X1} { set YMN [expr $SDJUNC+$JUNCSP*(1.5+$RATIO)] }
 if {$E} {
   set XMN  [expr $DSLDIF-$Z1]
   set XMX  [expr $DSLDIF+$Z1]
 }
 if {$E1} {
  set XMN  [expr $DSSDIF-$Z1]
  set XMX  [expr $DSSDIF+$Z1]
 }
 if {$D}  { puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"}
 if {$E} {
   set XMN  [expr $LMAX-$DDLDIF-$Z1]
   set XMX  [expr $LMAX-$DDLDIF+$Z1]
 }
 if {$E1} {
  set XMN  [expr $LMAX-$DDDDIF-$Z1]
  set XMX  [expr $LMAX-$DDDDIF+$Z1]
 }
 if {$D} { puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN "}
}

#-----------------------
if {$H1} {
  set XMN  [expr $DSSDIF-$Z1]
  set XMX  [expr $DSLDIF+$Z1]

  if {$E}  {  set YMN  [expr $SDJUNC+$Z2]}
  if {$E1} {  set YMN  [expr $LDDJUNC+$Z2]}

  puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"

  set XMN  [expr $LMAX-$DDLDIF-$Z1]
  set XMX  [expr $LMAX-$DDDDIF+$Z1]
  puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"
}

#-----------------------
#COMMENT   Second pass of eliminates.

if {$LDF} { set DSLDIF  [expr $DSSDIF]}
if {$LDF} { set DDLDIF  [expr $DDDDIF]}
if {$D1}  { set DSLDIF  [expr $DSSDIF]}
if {$D1}  { set DDLDIF  [expr $DDDDIF]}

set WD    [expr ($DSLDIF-$DSSDIF)/2]

set SPC  [expr $JUNCSP]
set TMP  [expr int($WD/$JUNCSP)]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X1   1
} else {
     set X1   0
}

if {$X1}  { set SPC  [expr $WD-($WD-$JUNCSP)/$RATIO] }

set WD  [expr (($RATIO+$EPS)*4*$SPC-$JUNCSP)/($RATIO+$EPS-1)]
set YMN [expr $JUNC1+$WD]

#COMMENT   Only do the eliminate if YMN < SUBDEP.
set TMP  [expr int($YMN/$SUBDEP)]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X1   0
} else {
     set X1   1
}

if {$X1} { set XMX  [expr $DSLDIF+$WD] }

#COMMENT   See if a single eliminate can be done or if we need to
#+         do separate eliminates for the left and right sides.
set X2  0
set X3  0

if {$X1} {
set TMP  [expr int($XMX/$LMID)]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
     set X2   0
     set X3   1
} else {
     set X2   1
     set X3   0
}
}

#COMMENT   Separate eliminates for the left and right.
if {$X2} {
 set XMN  [expr $DSSDIF-$WD]
 puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"
 set XMN  [expr $LMAX-$DDLDIF-$WD]
 set XMX  [expr $LMAX-$DDDDIF+$WD]
 puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"
}

#COMMENT   Single eliminate.
if {$X3} {
 set XMN  [expr $DSSDIF-$WD-$EPS]
 set XMX  [expr $LMAX-$DDDDIF+$WD+$EPS]
 puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"
}


#-----------------------
#  REGIONS
#-----------------------
puts $FileID  "#---------------------------------------------------"
puts $FileID  "REGION  Label=Silicon  Material=Si"
puts $FileID  "REGION  Label=Oxide    Material=Ox  IY.MAX=2"

#-----------------------
#  ELECTRODES
#-----------------------
puts $FileID  "REGION  Label=Drain   IY.MAX=2  X.MIN=$LMAX-$LDCONT Material=Elec"
puts $FileID  "REGION  Label=Gate    IY.MAX=1  X.MIN=$LSOURCE  X.MAX=$LMAX-$LDRAIN Material=Elec" 
puts $FileID  "REGION  Label=Source  IY.MAX=2  X.MAX=$LSCONT Material=Elec"
puts $FileID  "REGION  Label=Substrate  Y.MAX=-$SUBDEP+0.0001 Material=Elec"

#-----------------------
#  PROFILES
#-----------------------
if {$VTTYPE == "P"}  {
  set VTTYP   "Acceptor"
}  else {
  set VTTYP   "Donor"
}

if {$NMOS} {
 set SDTYP   "Donor"
 set SUBTYP  "Acceptor"
} 
if {$PMOS} {
 set SDTYP   "Acceptor"
 set SUBTYP  "Donor"
}

#-----------------------
#  ANALYTIC PROFILES
#-----------------------
puts $FileID  "#---------------------------------------------------"
puts $FileID  "PROFILE  Type=Uniform Ion=$SUBTYP  N.PEAK=$NSUB  X.MIN=0.0  X.MAX=$LMAX Y.TOP=0.0 Y.BOTTOM=-$SUBDEP"

if {$VTPEAK > 1.0} {
    puts $FileID  "PROFILE  Type=Gauss Ion=$VTTYP  N.PEAK=$VTPEAK  Y.CHAR=$VTCHAR X.MIN=0.0  X.MAX=$LMAX"
}

if {$LDT} {
  puts $FileID  "PROFILE  Type=Gauss Ion=$SDTYP  N.PEAK=$LDDPEAK  Y.CHAR=$LDYC X.CHAR=[expr $LATD*$LDYC] \
                          X.MIN=0.0  X.MAX=$LSOURCE"
  puts $FileID  "PROFILE  Type=Gauss Ion=$SDTYP  N.PEAK=$LDDPEAK  Y.CHAR=$LDYC X.CHAR=[expr $LATD*$LDYC] \
                          X.MIN=[expr $LMAX-$LDRAIN]  X.MAX=[expr $LMAX-$LDRAIN+$LDRAIN]"
}

puts $FileID  "PROFILE  Type=Gauss Ion=$SDTYP  N.PEAK=$SDPEAK  Y.CHAR=$SDYC  X.CHAR=[expr $LATD*$SDYC] \
                        X.MIN=0.0  X.MAX=[expr $LSOURCE-$LSPACER]"
puts $FileID  "PROFILE  Type=Gauss Ion=$SDTYP  N.PEAK=$SDPEAK  Y.CHAR=$SDYC  X.CHAR=[expr $LATD*$SDYC] \
                        X.MIN=[expr $LMAX-$LDRAIN+$LSPACER]  X.MAX=[expr $LMAX-$LDRAIN+$LSPACER+$LDRAIN-$LSPACER]"


#-----------------------
#  Boundary Conditions
#-----------------------
#  WF   =  Work Function (V)
if {$GATE == "PPoly"} {
 set WF 5.25
} elseif {$GATE == "NPoly"} {
 set WF 4.17
} elseif {$GATE == "Al"} {
 set WF 4.15
} 

puts $FileID  "#---------------------------------------------------"
puts $FileID  "BOUNDARY Type = InsulatorInterface ID = IF_Oxide_to_Silicon QF=0"
puts $FileID  "CONTACT  Type = GateContact        ID = Gate        WorkFunction=$WF"
puts $FileID  "CONTACT  Type = OhmicContact       ID = Substrate"
puts $FileID  "CONTACT  Type = OhmicContact       ID = Source"
puts $FileID  "CONTACT  Type = OhmicContact       ID = Drain"

#-----------------------
#  PLOTS
#-----------------------
puts $FileID  "#---------------------------------------------------"
puts $FileID  "PLOT Variable=DeviceMesh"

#-----------------------
#  Solve Equilibrium State
#-----------------------
puts $FileID  "METHOD   Type = DDML1   Scheme = Newton  NS=Basic LS=GMRES"
puts $FileID  "SOLVE    Type = TRANSIENT TStart = 0 TStep=1e-11  TStop = 3e-10"
puts $FileID  "SOLVE    Type = EQUILIBRIUM"
puts $FileID  "EXPORT   CoreFile = $TRANTYPE.init.cgns"
puts $FileID  "PLOT     Variable=Na        Resolution=RES.High    AzAngle=120  ElAngle=60"
puts $FileID  "PLOT     Variable=Nd        Resolution=RES.High    AzAngle=120  ElAngle=60"
puts $FileID  "PLOT     Variable=Potential Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "PLOT     Variable=ElecDensity      Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "PLOT     Variable=HoleDensity      Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "END"

close $FileID

