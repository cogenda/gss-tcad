# We use tcl script to generate a BJT structure
# run "tclsh bjt.tcl" to generate "bjt.inp"

#set the cmd file name for GSS 
set FileName "bjt.inp"

#load default BJT structure data file
source  "$env(GSS_DIR)/lib/template/bipdef.tcl"

# load template tcl script to generate input file
source  "$env(GSS_DIR)/lib/template/bipstr.tcl"

# after that, call gss
#exec $env(GSS_DIR)/bin/gss $FileName
