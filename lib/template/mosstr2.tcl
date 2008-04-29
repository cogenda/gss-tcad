#-------------------
#  mosstr1
#-------------------
#  Template File for MOS Structures

puts "Generate GSS input file..."

#COMMENT   Derived Quantities and Definitions
#    WD       = depletion width under the drain (microns)
#    SUBDEP   = substrate depth (microns)
#    LMID     = distance to middle 
#		of gate from source edge (microns)
#    LMAX     = width of device (microns)
#    EPS      = roundoff error
set WD        [expr 3.6E7*sqrt((.8+abs($VDBMAX))/$NSUB)]
set SUBDEP    [expr int(2*($SDJUNC+$WD+1.5))/2]
set LMID      [expr $LSOURCE+$LGATE/2]
set LMAX      [expr $LSOURCE+$LGATE+$LDRAIN]
set EPS       [expr 0.0002]

#COMMENT   Only consider LDD implants if LDDPEAK >= NSUB.
#         LDT ==> LDD implant
#         LDF ==> No LDD implant
if {$LDDPEAK>$NSUB} {
  set LDT  1
  set LDF  0
} else {
  set LDT  0
  set LDF  1
}

#COMMENT   Transistor type.
if  {$TRANTYPE == "PMOS"} {
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
#    LDDIF    = lightly doped drain out-diffusion
#               distance (microns)

set NSRF  [expr $NSUB+$VTPEAK]
if {$NMOS && $VTN} { set NSRF  [expr abs($NSUB-$VTPEAK)]}
if {$PMOS && $VTP} { set NSRF  [expr abs($NSUB-$VTPEAK)]}

set NSDJ  [expr $NSUB]
set TMP   [expr $VTPEAK*exp(- pow($SDJUNC/$VTCHAR,2))]
if {$NMOS && $VTP} { set NSDJ  [expr $NSUB+$TMP]}
if {$PMOS && $VTN} { set NSDJ  [expr $NSUB+$TMP]}
set SDYC  [expr $SDJUNC/sqrt(log($SDPEAK/$NSDJ))]
if {$LDT} {set NLDJ  [expr $NSUB]}
if {$LDT} {set TMP   [expr $VTPEAK*exp(- pow($LDDJUNC/$VTCHAR,2))]}
if {$LDT && $NMOS && $VTP} {set NLDJ  [expr $NSUB+$TMP]}
if {$LDT && $PMOS && $VTN} {set NLDJ  [expr $NSUB+$TMP]}
if {$LDT} {set LDYC  [expr $LDDJUNC/sqrt(log($LDDPEAK/$NLDJ))]}

set FileID [open $FileName w]
puts "Write GSS input statements to $FileName"

#-----------------------
#  MESH GENERATION
#-----------------------
puts $FileID "MESH   Type=GSS   ModelFile=$TRANTYPE.cgns Triangle=pzAY"
puts $FileID "XMESH  WIDTH=$LSCONT           N.SPACES=[expr int($LSCONT/0.1)]"
puts $FileID "XMESH  WIDTH=$LSOURCE-$LSCONT  N.SPACES=[expr int(($LSOURCE-$LSCONT)/0.1)]"
puts $FileID "XMESH  WIDTH=$LGATE/2          N.SPACES=[expr int($LGATE/0.2)]"
puts $FileID "XMESH  WIDTH=$LGATE/2          N.SPACES=[expr int($LGATE/0.2)]"
puts $FileID "XMESH  WIDTH=$LDRAIN-$LDCONT   N.SPACES=[expr int(($LDRAIN-$LDCONT)/0.1)]"
puts $FileID "XMESH  WIDTH=$LDCONT           N.SPACES=[expr int($LDCONT/0.1)]"

puts $FileID "YMESH  Y.TOP=0.1+$TOX  Y.BOTTOM=$TOX  N.SPACES=1"
puts $FileID "YMESH  Y.BOTTOM=0.0    N.SPACES=1"
puts $FileID "YMESH  DEPTH=$SDJUNC   H1=$SDJUNC/8  H2=$SDJUNC/4"
puts $FileID "YMESH  DEPTH=1-$SDJUNC H1=$SDJUNC/4  H2=0.2"
puts $FileID "YMESH  Y.BOTTOM=-$SUBDEP  H1=0.2  RATIO=1.2"
puts $FileID "YMESH  Y.BOTTOM=-$SUBDEP-0.1  N.SPACES=1"

puts $FileID "ELIMINATE  Direction=COLUMNS  Y.TOP=-1.5"

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
puts $FileID  "REGION  Label=Gate    IY.MAX=1  X.MIN=$LSOURCE  X.MAX=$LSOURCE+$LGATE Material=Elec" 
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
puts $FileID  "PROFILE  Ion=$SUBTYP  N.PEAK=$NSUB  Type=Uniform  X.MIN=0.0  X.MAX=$LMAX Y.TOP=0.0 Y.BOTTOM=-$SUBDEP"
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
#  REGRIDS
#-----------------------
# Doping regrids
puts $FileID  "REFINE   Variable=Doping Measure=SignedLog Dispersion=3  Triangle=praz"   

# Potential regrids


#-----------------------
#  PLOTS
#-----------------------			  
puts $FileID  "#---------------------------------------------------"
puts $FileID  "PLOT Variable=DeviceMesh"

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
puts $FileID  "PLOT     Variable=Na        Resolution=RES.High    AzAngle=120  ElAngle=60"
puts $FileID  "PLOT     Variable=Nd        Resolution=RES.High    AzAngle=120  ElAngle=60"
puts $FileID  "PLOT     Variable=Potential Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "PLOT     Variable=ElecDensity      Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "PLOT     Variable=HoleDensity      Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "EXPORT   CoreFile = $TRANTYPE.init.cgns"
puts $FileID  "END"


close $FileID


