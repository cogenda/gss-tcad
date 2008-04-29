#-------------------
#  mosdef
#-------------------
#  Default Value File for MOS Structures

puts "Loading default MOS parameters..."
puts "These settings can be overwrite by user"



#COMMENT   Structure Definitions
#    LGATE    = gate length (microns)
#    LSOURCE  = distance from left device edge 
#               to gate edge (microns)
#    LSCONT   = length of source contact (microns)
#    LDRAIN   = distance from right device edge 
#               to gate edge (microns)
#    LDCONT   = length of drain contact (microns)
#    LSPACER  = spacer thickness (microns)
#    TOX      = gate oxide thickness (microns)

set LGATE     [ expr 1.0]
set LSOURCE   [ expr 1.0]
set LSCONT    [ expr 0.5]
set LDRAIN    [ expr 1.0]
set LDCONT    [ expr 0.5]
set LSPACER   [ expr 0.2]
set TOX       [ expr 0.0250]

#COMMENT   Doping Information
#    TRANTYPE = transistor type (NMOS or PMOS)
#    PROFTYPE = profile type (ANALYTIC, SUPREM3, TSUPREM4, 
#               or SUPRA)
#    LATD     = source/drain and LDD lateral diffusion factor
#  Analytic Profile Parameters
#    NSUB     = substrate doping (#/cm^3)
#    VTTYPE   = doping type for threshold adjust implant (N or P)
#    VTPEAK   = peak doping for threshold adjust implant (#/cm^3)
#    VTCHAR   = characteristic length 
#               for threshold implant (microns)
#    SDPEAK   = peak doping for source/drain (#/cm^3)
#    SDJUNC   = junction depth for source/drain (microns)
#    LDDPEAK  = peak doping for lightly doped drain (#/cm^3)
#    LDDJUNC  = junction depth for lightly doped drain (microns)


set TRANTYPE  NMOS
set PROFTYPE  ANALYTIC
set LATD      [expr 0.80]

set NSUB      [expr 3E15]
set VTTYPE    P
set VTPEAK    [expr 2E16]
set VTCHAR    [expr 0.25]
set SDPEAK    [expr 1E20]
set SDJUNC    [expr 0.25]
set LDDPEAK   [expr 2E18]
set LDDJUNC   [expr 0.35]


#COMMENT   Grid Spacings, Ratio, Maximum Voltage
#    CHANSP   = vertical grid spacing in the channel (microns)
#    JUNCSP   = grid spacing at junctions (microns)
#    RATIO    = grid spacing ratio
#    VDBMAX   = maximum drain-substrate reverse bias (volts)

set CHANSP    [expr 0.0125]
set JUNCSP    [expr 0.0250]
set RATIO     [expr 1.4]
set VDBMAX    [expr 5.0]

# COMMENT   Material of gate 
#    GATE    = gate material, default: PolySi 
if {$TRANTYPE == "PMOS"} {
 set GATE PPoly
} else {
 set GATE NPoly
}
