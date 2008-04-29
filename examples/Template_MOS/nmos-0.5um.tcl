# We use tcl script to generate a 0.5um NMOS with LDD structure
# run "tclsh nmos-0.5um.tcl" to generate "nmos-0.5u.inp"

# set the cmd file name for GSS 
set FileName "nmos-0.5u.inp"

# load default MOS structure data file (defined for 1um MOS)
source  "$env(GSS_DIR)/lib/template/mosdef.tcl"

# Modify the mosdef0 parameters to define a 0.5 micron n-channel LDD
set    TRANTYPE  NMOS
set    LGATE     0.5
set    LSOURCE   0.5
set    LSCONT    0.25
set    LDRAIN    0.5
set    LDCONT    0.25
set    LSPACER   0.15
set    TOX       0.01
set    NSUB      5E15
set    GATE      PPoly
set    VTTYPE    P
set    VTPEAK    5E15
set    VTCHAR    0.10
set    SDPEAK    1E20
set    SDJUNC    0.20
set    LDDPEAK   2E18
set    LDDJUNC   0.1
set    CHANSP    .005
set    JUNCSP    .03
set    RATIO     1.4
set    VDBMAX    3

# load template tcl script to generate input file
# one can try mosstr1.tcl or mosstr2.tcl
source  "$env(GSS_DIR)/lib/template/mosstr1.tcl"
#source  "$env(GSS_DIR)/lib/template/mosstr2.tcl"

# after that, call gss
#exec $env(GSS_DIR)/bin/gss $FileName
