#-------------------
#  bipdef0
#-------------------
#  Default Value File for BJT Structure

puts "Loading default BJT parameters..."
puts "These settings can be overwrite by user"

#COMMENT   Structure Definitions
#    WEMIT    = emitter width (microns)
#    WECONT   = emitter contact width (microns)
#    WXBASE   = extrinsic base width (microns)
#    WBCONT   = base contact width (microns)
#    WEXB     = emitter to extrinsic base distance (microns)
#    BEOVER   = base diffusion edge 
#               to emitter diffusion edge (microns)
#    BXOVER   = base diffusion edge 
#               to x-base diffusion edge (microns)
#    TEPI     = structure depth (epi thickness) (microns)

set    WEMIT     1.0
set    WECONT    1.0
set    WXBASE    1.0
set    WBCONT    1.0
set    WEXB      0.5
set    BEOVER    0.5
set    BXOVER    0.5
set    TEPI      2.0

#COMMENT   Doping Information
#    TRANTYPE = transistor type (NPN or PNP)
#    PROFTYPE = profile type (ANALYTIC, SUPREM3, TSUPREM4,
#                             or SUPRA)
#    LATD     = lateral diffusion factor
#  Analytic Profile Parameters
#    NEPI     = epitaxial layer doping (#/cm^3)
#    BLPEAK   = peak doping for buried layer (#/cm^3)
#    BLDEPTH  = depth of buried layer (microns)
#    BPEAK    = peak doping for intrinsic base (#/cm^3)
#    YBPEAK   = distance from surface 
#               to peak base doping (microns)
#    BCJUNC   = base-collector junction depth (microns)
#    XBPEAK   = peak doping for extrinsic base (#/cm^3)
#    XBJUNC   = extrinsic base-collector junction depth (microns)
#    EPEAK    = peak doping for emitter (#/cm^3)
#    EBJUNC   = emitter-base junction depth (microns)


set    TRANTYPE  NPN
set    PROFTYPE  ANALYTIC
set    LATD      0.80

set    NEPI      1E16
set    BLPEAK    1E19
set    BLDEPTH   0.5
set    BPEAK     4E17
set    YBPEAK    0.0
set    BCJUNC    0.40
set    XBPEAK    5E19
set    XBJUNC    0.45
set    EPEAK     1E20
set    EBJUNC    0.10


#COMMENT   Grid Spacings, Ratio, Maximum Voltage
#    EBSP     = grid spacing at emitter-base junction (microns)
#    BCSP     = grid spacing at collector-base junction (microns)
#    RATIO    = grid spacing ratio
#    VCBMAX   = maximum collector-base reverse bias (volts)

set    EBSP      0.0125
set    BCSP      0.0250
set    RATIO     1.5
set    VCBMAX    3

