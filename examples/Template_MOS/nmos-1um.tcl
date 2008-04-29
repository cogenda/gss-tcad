# We use tcl script to generate a 1um NMOS with LDD structure
# run "tclsh nmos-1um.tcl" to generate "nmos-1u.inp"

# set the cmd file name for GSS 
set FileName "nmos-1u.inp"

# load default MOS structure data file (defined for 1um MOS)
source  "$env(GSS_DIR)/lib/template/mosdef.tcl"

#set    VTTYPE    P
#set    VTPEAK    [expr 1E14]
#set    VTCHAR    [expr 0.25]
#set    GATE      PPoly
# load template tcl script to generate input file
# one can try mosstr1.tcl or mosstr2.tcl
source  "$env(GSS_DIR)/lib/template/mosstr1.tcl"
#source  "$env(GSS_DIR)/lib/template/mosstr2.tcl"

# after that, call gss
#exec $env(GSS_DIR)/bin/gss $FileName
