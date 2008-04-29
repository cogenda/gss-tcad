#-------------------
#  bipstr
#-------------------
#  Template File for BJT Structures

puts "Generate GSS input file..."

#COMMENT   Only consider X-base implants if XBPEAK >= 1.
if {$XBPEAK>1.0} {
  set  XBT  1
} else {
  set  XBT  0
}

if {$XBT} {
set  XYCH   [expr ($XBJUNC)/sqrt(log($XBPEAK/$NEPI))]
set  XXCH   [expr $LATD*$XYCH]
set  XDIF   [expr $XXCH*sqrt(log($XBPEAK/$NEPI))]
}
#COMMENT   Derived Quantities and Definitions
#    EPS     = roundoff error
#    WD      = collector-base depletion width (microns)
#    BYCH    = base diffusion vertical 
#              characteristic length (microns)
#    BXCH    = base diffusion horizontal 
#	       characteristic length (microns)
#    BDIF    = base lateral diffusion distance (microns)
#    XDIF    = lateral diffusion distance beyond BXOVER (microns)
#    DLBDE   = distance from left device edge 
#              to base diffusion edge
#    DRBDE   = distance from right device edge 
#              to base diffusion edge
#    WBASE   = width of base
#    WMAX    = width of device (microns)
#    LEDE    = location of emitter diffusion edge
#    LECONT  = location of emitter contact
#    LBCONT  = location of base contact
set  EPS    [expr .0002]
set  YBPPE  [expr $YBPEAK+$EPS]
set  WD     [expr 3.6E7*sqrt((.8+abs($VCBMAX))/$NEPI)]
set  WBASE  [expr $BEOVER+$WEMIT+$WEXB+$WXBASE+$BXOVER]
set  BYCH   [expr ($BCJUNC-$YBPEAK)/sqrt(log($BPEAK/$NEPI))]
set  BXCH   [expr $LATD*$BYCH]
set  TMP    [expr $BPEAK*exp(- pow($YBPPE/$BYCH,2))]
set  BDIF   [expr $BXCH*sqrt(log($TMP/$NEPI))]
set  NSRF   [expr $TMP-$NEPI]
set  DLBDE  [expr int(2*($BDIF+$WD+1.0))/2]
set  LEDE   [expr $DLBDE+$BEOVER]
set  LECONT [expr $LEDE+($WEMIT-$WECONT)/2]
set  LBCONT [expr $LEDE+$WEMIT+$WEXB+($WXBASE-$WBCONT)/2]
set  TMP    [expr $BPEAK*exp(- pow(($YBPPE-$EBJUNC)/$BYCH,2))]
set  NE     [expr $TMP-$NEPI]
set  EYCH   [expr ($EBJUNC)/sqrt(log($EPEAK/$NE))]
set  EXCH   [expr $LATD*$EYCH]
set  EDIF   [expr $EXCH*sqrt(log($EPEAK/$NSRF))]

#COMMENT See if XDIF < BDIF+BXOVER+BCSP
if {$XBT} {
 set  TMP  [expr int($XDIF/($BDIF+$BXOVER+$BCSP))]
 if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
  set X1   0 
  set X2   1 
 } else {
  set X1   1
  set X2   0 
 }
}

if {$X1} {  set  XDIF  [expr $BDIF]}
if {$X2} {  set  XDIF  [expr $XDIF-$BXOVER]}
set  DRBDE  [expr int(2*($XDIF+$WD+1.0))/2]
set  WMAX   [expr $DLBDE+$WBASE+$DRBDE]

#COMMENT   Transistor type.
if {$TRANTYPE == "PNP"} {
  set  NPN  0
  set  PNP  1
} else {
  set  NPN  1
  set  PNP  0
}

set FileID [open $FileName w]
puts "Write GSS input statements to $FileName"

#-----------------------
#  MESH statement
#-----------------------
puts $FileID "MESH     Type=GSS   ModelFile=$TRANTYPE.cgns Triangle=pzAY"
puts $FileID  "#---------------------------------------------------"

#-----------------------
#  X.MESH statements
#-----------------------
puts $FileID "XMESH    WIDTH=$DLBDE-$BDIF  H2=$BCSP  RATIO=$RATIO"
set  TMP   [expr 0.5*($BCSP-$EBSP)/($RATIO-1)]
set  LX1   [expr ($BDIF+$BEOVER-$EDIF)/2+$TMP]
set  LX2   [expr ($BDIF+$BEOVER-$EDIF)/2-$TMP]
puts $FileID "XMESH    WIDTH=$LX2  H1=$BCSP  RATIO=$RATIO"
puts $FileID "XMESH    WIDTH=$LX1  H2=$EBSP  RATIO=$RATIO"

set    W1    [expr $EDIF+($WEMIT-$WECONT)/2]
set    SPC1  [expr $W1-($W1-$EBSP)/$RATIO]
puts $FileID "XMESH    WIDTH=$W1  H1=$EBSP  RATIO=$RATIO"
puts $FileID "XMESH    WIDTH=$WECONT/2  H1=$SPC1  RATIO=$RATIO"
puts $FileID "XMESH    WIDTH=$WECONT/2  H2=$SPC1  RATIO=$RATIO"
puts $FileID "XMESH    WIDTH=$W1  H2=$EBSP  RATIO=$RATIO"

set    W1    [expr $WEXB-$EDIF+($WXBASE-$WBCONT)/2]
set    SPC1  [expr $W1-($W1-$EBSP)/$RATIO]
set    W2    [expr $XDIF+$BXOVER+($WXBASE-$WBCONT)/2]
set    SPC2  [expr $W2-($W2-$BCSP)/$RATIO]

#COMMENT   See if SPC1+SPC2 < WBCONT
set    TMP  [expr int(($SPC1+$SPC2)/$WBCONT)]
if { [expr int((2+$TMP)/(1+$TMP))%2] == 1} {
  set X1   0 
  set X2   1 
} else {
  set X1   1
  set X2   0 
}
 

puts $FileID "XMESH    WIDTH=$W1  H1=$EBSP  RATIO=$RATIO"
if {$X1} { puts $FileID "XMESH    WIDTH=$WBCONT  H1=$SPC1  H2=$SPC2  RATIO=$RATIO"}
if {$X2} { puts $FileID "XMESH    WIDTH=$WBCONT  N.SPACES=2"}
puts $FileID "XMESH    WIDTH=$W2  H2=$BCSP  RATIO=$RATIO"
puts $FileID "XMESH    X.MAX=$WMAX  H1=$BCSP  RATIO=$RATIO"
puts $FileID  "#---------------------------------------------------"

#-----------------------
#  Y.MESH statements
#-----------------------
puts $FileID "YMESH    Y.BOTTOM=0  Y.TOP=0.1  N.SPACES=1"
puts $FileID "YMESH    Y.BOTTOM=-$EBJUNC  H2=$EBSP  RATIO=$RATIO"
set  TMP  [expr 0.5*($BCSP-$EBSP)/($RATIO-1)]
set  LY1  [expr ($BCJUNC-$EBJUNC)/2+$TMP]
set  LY2  [expr ($BCJUNC-$EBJUNC)/2-$TMP]
puts $FileID "YMESH    DEPTH=$LY1  H1=$EBSP  RATIO=$RATIO"
puts $FileID "YMESH    Y.BOTTOM=-$BCJUNC  H2=$BCSP  RATIO=$RATIO"
puts $FileID "YMESH    Y.BOTTOM=-$TEPI  H1=$BCSP  RATIO=$RATIO"
puts $FileID "YMESH    Y.BOTTOM=-$TEPI-0.1  N.SPACES=1"
puts $FileID  "#---------------------------------------------------"


#-----------------------
#  ELIMINATE ROWS
#-----------------------
set  W12    [expr $EBSP*(2*$RATIO-1)/($RATIO-1)]
set  W14    [expr $EBSP*(4*$RATIO-1)/($RATIO-1)]
set  W22    [expr $BCSP*(2*$RATIO-1)/($RATIO-1)]
set  W24    [expr $BCSP*(4*$RATIO-1)/($RATIO-1)]
set  WEJBC  [expr $WEXB-$EDIF+($WXBASE-$WBCONT)/2]
set  WXJBC  [expr $XDIF+$BXOVER+($WXBASE-$WBCONT)/2]
set  SFAC   [expr 1.6]
#-----------------------------------------------------------------------


for {set STEPS 0} {$STEPS < 2 } { incr STEPS } {
  if {$STEPS == 0} {
    set  W1   $W12
  } else {
    set  W1   $W14
  }   
  set  TMP  [expr int($W1/$EBJUNC)]
  if {$STEPS == 0} {
    set    X1  1
    set    X2  0
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
    set    X2  [expr !$X2]
  }
  if {$X1} { set  YMN  [expr $EBJUNC-$W1]}
  if {$X2} { set  YMN  [expr $EPS]}
  
  set  TMP  [expr int($W1/$LY1)]
  if {$STEPS == 0} {
    set    X1  1
    set    X2  0
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
    set    X2  [expr !$X2]
  }
  if {$X1} { set  YMX  [expr $EBJUNC+$W1]}
  if {$X2} { set  YMX  [expr $EBJUNC+$LY1-$EPS]}

  set  TMP  [expr int($W1/$LX1)]
  if {$STEPS == 0} {
    set    X1  1
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
  }
  
  set  TMP  [expr int($W1/$WEJBC)]
  if {$STEPS == 0} {
    set    X2  1
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X2  [expr !$X1]
  }
  
  if {$X1} { set  XMX   [expr $LEDE-$EDIF-$SFAC*$W1] }
  if {$X1} { puts $FileID  "ELIMINATE  Direction=ROWS  X.MAX=$XMX  Y.TOP=-$YMN  Y.BOTTOM=-$YMX" }
  if {$X2} { set  XMN  [expr $LEDE+$WEMIT+$EDIF+$SFAC*$W1] }
  if {$X2} { puts $FileID  "ELIMINATE  Direction=ROWS  X.MIN=$XMN  Y.TOP=-$YMN  Y.BOTTOM=-$YMX" }
}

#-----------------------------------------------------------------------
for {set STEPS 0} {$STEPS < 2 } { incr STEPS } {
  if {$STEPS == 0} {
    set  W1   $W22
  } else {
    set  W1   $W24
  }   
  set  TMP  [expr int($W1/$LY2)]
  if {$STEPS == 0} {
    set    X1  1
    set    X2  0
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
    set    X2  [expr !$X2]
  }
  if {$X1} { set  YMN  [expr $BCJUNC-$W1]}
  if {$X2} { set  YMN  [expr $EPS]}
  
  set  YMX  [expr $BCJUNC+$W1]
  
  set  XMX   [expr $DLBDE-$BDIF-$SFAC*$W1]
  puts $FileID  "ELIMINATE  Direction=ROWS  X.MAX=$XMX  Y.TOP=-$YMN  Y.BOTTOM=-$YMX"
  set  XMN  [expr $WMAX-$DRBDE+$XDIF+$SFAC*$W1]
  puts $FileID  "ELIMINATE  Direction=ROWS  X.MIN=$XMN  Y.TOP=-$YMN  Y.BOTTOM=-$YMX"
}

#-----------------------
#  ELIMINATE COLUMNS
#-----------------------
for {set STEPS 0} {$STEPS < 2 } { incr STEPS } {
  if {$STEPS == 0} {
    set  W1   $W12
  } else {
    set  W1   $W14
  }   
  set  YMN  [expr $EBJUNC+$SFAC*$W1]
 
  set  TMP  [expr int($W1/$LX1)]
  if {$STEPS == 0} {
    set    X1  1
    set    X2  0
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
    set    X2  [expr !$X2]
  }
  if {$X1} { set  XMN  [expr $LEDE-$EDIF-$W1]}
  if {$X2} { set  XMN  [expr $DLBDE+$LX2+$EPS]}

  set  TMP  [expr int($W1/($EDIF+$WEMIT/2)) ] 
  if {$STEPS == 0} {
    set    X1  1
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
  }
  
  if {$X1} { set  XMX  [expr $LEDE-$EDIF+$W1]}
  if {$X1} { puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"}
  if {$X1} { set  XMN  [expr $LEDE+$WEMIT+$EDIF-$W1]}

  set  TMP  [expr int($W1/$WEJBC)]
  if {$STEPS == 0} {
    set    X1  1
    set    X2  0
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
    set    X2  [expr !$X2]
  }
  if {$X1} { set  XMX  [expr $LEDE+$WEMIT+$EDIF+$W1]}
  if {$X2} { set  XMX  [expr $LBCONT-$EPS]}
  puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"
}

#-----------------------------------------------------------------------
for {set STEPS 0} {$STEPS < 2 } { incr STEPS } {
  if {$STEPS == 0} {
    set  W1   $W22
  } else {
    set  W1   $W24
  }   
  set  YMN  [expr $BCJUNC+$SFAC*$W1]
  set  XMN  [expr $DLBDE-$BDIF-$W1]
  set  TMP  [expr int($W1/$LX2)]
  if {$STEPS == 0} {
    set    X1  1
    set    X2  0
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
    set    X2  [expr !$X2]
  }
 
  if {$X1} { set  XMX  [expr $DLBDE-$BDIF+$W1]}
  if {$X2} { set  XMX  [expr $DLBDE-$BDIF+$LX2-$EPS]}
  puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"

  set  TMP  [expr int($W1/$WXJBC)]
  if {$STEPS == 0} {
    set    X1  1
    set    X2  0
  }  
  for {set loop 0} { $loop < [expr int((2+$TMP)/(1+$TMP))]} {incr loop} {
    set    X1  [expr !$X1]
    set    X2  [expr !$X2]
  }
  if {$X1} { set  XMN  [expr $WMAX-$DRBDE+$XDIF-$W1]}
  if {$X2} { set  XMN  [expr $WMAX-$DRBDE+$XDIF-$WXJBC+$EPS]}
  set  XMX  [expr $WMAX-$DRBDE+$XDIF+$W1]
  puts $FileID  "ELIMINATE  Direction=COLUMNS  X.MIN=$XMN  X.MAX=$XMX  Y.TOP=-$YMN"
}

#-----------------------
#  REGIONS
#-----------------------
puts $FileID  "#---------------------------------------------------"
puts $FileID  "REGION    Label=Silicon    Material=Si"
puts $FileID  "REGION    Label=Oxide1     Material=Ox    IY.MIN=0 IY.MAX=1 X.MIN=0 X.MAX=$LBCONT"

#-----------------------
#  ELECTRODES
#-----------------------
puts $FileID  "REGION    Label=Collector  Material=Elec Y.TOP=-$TEPI+0.0001" 
puts $FileID  "REGION    Label=Emitter    Material=Elec IY.MAX=1  X.MIN=$LECONT  X.MAX=$LECONT+$WECONT"
puts $FileID  "REGION    Label=Oxide2     Material=Ox   IY.MAX=1  X.MIN=$LECONT+$WECONT  X.MAX=$LBCONT"
puts $FileID  "REGION    Label=Base       Material=Elec IY.MAX=1  X.MIN=$LBCONT  X.MAX=$LBCONT+$WBCONT"
puts $FileID  "REGION    Label=Oxide3     Material=Ox   IY.MAX=1  X.MIN=$LBCONT+$WBCONT  X.MAX=$WMAX"
#-----------------------
#  PROFILES
#-----------------------
set  BLCHAR  [expr $BLDEPTH/sqrt(log($BLPEAK/$NEPI))]
if {$NPN} {
 set  ECTYP  "Donor"
 set  BTYP   "Acceptor"
}

if {$PNP} {
 set  ECTYP  "Acceptor"
 set  BTYP   "Donor"
}

#-----------------------
#  ANALYTIC PROFILES
#-----------------------
puts $FileID  "#---------------------------------------------------"
puts $FileID  "PROFILE  Ion=$ECTYP  N.PEAK=$NEPI  Type=Uniform  X.MIN=0 X.MAX=$WMAX Y.TOP=0.1 Y.BOTTOM=-$TEPI-0.1"

if {$BLPEAK>1.0} {
puts $FileID  "PROFILE  Type=Gauss Ion=$ECTYP  N.PEAK=$BLPEAK  Y.TOP=-$TEPI  Y.CHAR=$BLCHAR X.MIN=0 X.MAX=$WMAX"
}

puts $FileID  "PROFILE  Type=Gauss Ion=$BTYP  N.PEAK=$BPEAK  Y.TOP=-$YBPEAK  Y.Junction=-$BCJUNC XY.RATIO=$LATD   \
                        X.MIN=$DLBDE  X.MAX=[expr $DLBDE+$WBASE]"
  
if {$XBT} {
puts $FileID  "PROFILE  Type=Gauss Ion=$BTYP  N.PEAK=$XBPEAK  Y.CHAR=$XYCH  X.CHAR=[expr $LATD*$XYCH] \
                        X.MIN=[expr $LEDE+$WEMIT+$WEXB]   X.MAX=[expr $LEDE+$WEMIT+$WEXB+$WXBASE]"
}  

puts $FileID  "PROFILE  Type=Gauss Ion=$ECTYP  N.PEAK=$EPEAK  Y.Junction=-$EBJUNC  XY.RATIO=$LATD  \
                        X.MIN=$LEDE  X.MAX=[expr $LEDE+$WEMIT]"

#-----------------------
#  Boundary Conditions
#-----------------------
puts $FileID  "#---------------------------------------------------"
puts $FileID  "CONTACT  Type = OhmicContact       ID = Emitter"
puts $FileID  "CONTACT  Type = OhmicContact       ID = Base"
puts $FileID  "CONTACT  Type = OhmicContact       ID = Collector"

#-----------------------
#  PLOTS
#-----------------------
puts $FileID  "#---------------------------------------------------"
puts $FileID  "PLOT Variable=DeviceMesh"

#-----------------------
#  Solve Equilibrium State
#-----------------------
puts $FileID  "METHOD   Type = DDML1   Scheme = Newton  NS=LineSearch LS=GMRES"
puts $FileID  "SOLVE    Type = TRANSIENT TStart = 0 TStep=1e-12  TStop = 3e-11"
puts $FileID  "SOLVE    Type = EQUILIBRIUM"
puts $FileID  "PLOT     Variable=Na        Resolution=RES.High    AzAngle=120  ElAngle=60"
puts $FileID  "PLOT     Variable=Nd        Resolution=RES.High    AzAngle=120  ElAngle=60"
puts $FileID  "PLOT     Variable=Potential Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "PLOT     Variable=ElecDensity      Resolution=RES.High    AzAngle=240  ElAngle=20"
puts $FileID  "PLOT     Variable=HoleDensity      Resolution=RES.High    AzAngle=240  ElAngle=20"

puts $FileID  "EXPORT   CoreFile = $TRANTYPE.init.cgns"
puts $FileID  "END"

close $FileID

