# Imen: This is the custom defaults file for SST-1M.
# == GrOptics
export telescope_def="$gropticsCfg/single.py"
export gro_default="$gropticsCfg/SST-1M"
export gro_geometry="$gropticsCfg/SST-1M.geo"

# == CARE
# CARE defaults
export care_default="$careCfg/SST-1M.cfg"
export telescopes_default="$careCfg/single_telescope.dat"
#export pulse="$careCfg/pulse_SST-1M_AfterPreamp"
export pulse="$careCfg/all_pulseshapes.dat"
export qe="$careCfg/qe-fresnel_care.dat"
export patches="" #"$careCfg/patches.dat" # Does it deserve the effort?
export camera="$careCfg/SST-1M_camera_descriptor.dat"
