#!/bin/bash

# Load .bashrc if $HOME is not accessible. Must be customize for the local batch system or commented out
#source $HOME/.bashrc
#source $HOME/.bash_profile
# Load the paths for this execution. Must be customize for the local tree
if [ "$thislocalsource" == "" ]
then
    localdir=`dirname $0`
    if [ -e $localdir/local.sh ]
    then
	source $localdir/local.sh
    else
	echo
	echo "ERROR: $localdir/local.sh doesn't exist! local.sh must be in the same folder of ${0}..."
	echo
	exit 1
    fi
else
    source "$thislocalsource"
fi
# Expand aliases so far defined in case the batch system doesn't import by default
shopt -s expand_aliases

if [ "$1" == "-h" ] || [ "$1" == "--help" ]
then
    echo
    echo "Launch a GrOptics simulation"
    echo "  usage: $0 <production> <run> <corsika_input> [ <output_name> <defaults> ] [ <cfg> <PyTelescopes> ] [ <shadowing> <atm_extintion> ]"
    echo
    echo "<production_name>   Name of the output (and base of the tmp) folder"
    echo "<run>               Run ID (used also as seed)"
    echo "<corsika_input>     Input file from CORSIKA files ('gunzip'ped!)"
    echo
    echo "<output_name>       Name of the output files (run name!) [Default: gropticsRun<run>]"
    echo "<defaults>          Default of the (following) settings [Default: $runPath/defaults.sh]"
    echo
    echo "<cfg>               GrOptics configuration name (WITHOUT EXTENSIONS!). The extensions '.pilot', '.array',"
    echo "                    '.cfg' and '.geo' will be added automatically [Default defined in <defaults>]"
    echo "<PyTelescopes>      Python file of the telescope(s) configuration(s) [Default defined in <defaults>]"
    echo
    echo "<shadowing>         Shadowing ON (1) or OFF (0) [Default: 1]"
    echo "<atm_extintion>     Eventual atmosphere extintion file [Default: '']"
    echo
    echo "Let's simulate the the ray tracing of the telescopes!"
    echo
    exit 1
elif [ $# -lt 3 ]
then
    echo
    echo "ERROR: at least 3 arguments (<name> <run> <grorptics_input>) must be given!"
    echo
    exit 1
fi

# Base setting (ARGS 1 and 2: NAME and RUN)
NAME="$1"
RUN="$2"
run=`printf "%06d" $RUN`
seed=$((RUN+1))

# Set the input path. N.B.: DO NOT USE corsika.gz!!!
coriska_input="$3"
# Set the output name
output_name="${4-gropticsRun$run}"
#if [ "$output_name" != "" ]
#then
#    output_name="${output_name}."
#fi
#output_name="${output_name}gropticsRun$run"

# Default settings
if [ "$5" == "" ]
then
    settings="$runPath/defaults.sh"
else
    settings="$5"
fi
source "$settings"

# Telescope configurations
if [ "$6" == "" ]
then
    config="$gro_default"
else
    config="$6"
fi
groPilot="$config".pilot
groCfg="$config".cfg
groArray="$config".array

if [ "$7" == "" ]
then
    telescope_def="$telescope_def"
else
    telescope_def="$7"
fi

# Other GrOptics controls
shadowing="${8-1}"

# Atmosphere config
atm_extintion="${9-}"

# Set the temporary folder/configs
tmp_folder="${tmp_dir}/${NAME}_${PBS_JOBID-`date +%d%m%y%H%M%S%N`}"
tmp_cfg="$tmp_folder/cfg"
mkdir -p "$tmp_cfg"
cfgarray="$tmp_cfg/groptics.array"
cfg="$tmp_cfg/groptics.cfg"
pilot="$tmp_cfg/groptics.pilot"
geometry="$tmp_cfg/geometry.geo"
cp "$gro_geometry" "$geometry"
# - only to keep record of this config in the .cfg.tar file ($atm_extintion will be used in the execution)
if [ "$atm_extintion" != "" ]
then
    cp "$atm_extintion" "$tmp_cfg/atm.extintion"
fi

# build array configurations
python "$telescope_def" "$tmp_cfg"
telescope_array="$tmp_cfg/telescopes.corsikaIO"
mirrors="$tmp_cfg/telescopes.mirrors"
arraytels="$tmp_cfg/telescopes.arraytel"
ntel=`grep '* TLLOC' "$telescope_array" | wc -l`

# Set output folder and tree structure
output="$out_dir/$NAME/GrOptics/$output_name"
badgroptics="$out_dir/$NAME/GrOptics/Bads/$output_name"
badcIOr="$out_dir/$NAME/corsikaIOreader/Bads/$output_name"
mkdir -p "$output"
mkdir -p "$badgroptics"
mkdir -p "$badcIOr"

# grOptics must be run in the $gropticsPath folder!
cd "$gropticsPath"

# Starting info
echo "#-------------------------------------------"
echo "#"
echo "#            START (run ${RUN})"
echo "#"
echo "#  grOptics on ${USER}@${HOSTNAME}"
echo "#  time: `date`"
echo "#"
if [ $PBS_JOBID ]
then
#    echo "#         job host: $PBS_JOBID"
    echo "#           job id: $PBS_JOBID"
    echo "#         job name: $PBS_JOBNAME"
    echo "#      job logname: $PBS_O_LOGNAME"
    echo "#       ...@ queue: $PBS_O_QUEUE"
fi
echo "#    CORSIKA input: ${coriska_input}"
echo "#"
echo "#           ...cwd: '`pwd`'"
echo "#"
echo "#-------------------------------------------"
echo ""
echo ""

# Build the local configs
# - .cfg
SEDcfg="s=@shadowing@=${shadowing}=g;s=@outpath@=${tmp_folder}=g;s=@gro_geometry@=${geometry}=g;s=@ntel@=${ntel}=g;/@telescope_array@/ {
r ${telescope_array}
d
};/@mirrors@/ {
r ${mirrors}
d
}"
cat <<EOF 1>&2
 -- Local setting of the configuration:
  - "$groCfg" to "$cfg":
      sed -e "$SEDcfg" "$groCfg" > "$cfg"

EOF
sed -e "$SEDcfg" "$groCfg" > "$cfg"
# - .array
SEDarray="s=@outpath@=${tmp_folder}=g;s=@cfgarray@=${cfgarray}=g;s=@cfg@=${cfg}=g;/@arraytels@/ {
r ${arraytels}
d
}"
cat <<EOF 1>&2
  - "$groArray" to "$cfgarray":
      sed -e "$SEDarray" "$groArray" > "$cfgarray"

EOF
sed -e "$SEDarray" "$groArray" > "$cfgarray"
# - .pilot
SEDpilot="s=@outpath@=${tmp_folder}=g;s=@cfgarray@=${cfgarray}=g;s=@seed@=${seed}=g"
cat <<EOF 1>&2
  - "$groPilot" to "$pilot":
      sed -e "$SEDpilot" "$groPilot" > "$pilot"

EOF
sed -e "$SEDpilot" "$groPilot" > "$pilot"
rm "$mirrors" "$arraytels"

# EXECUTIONS
# - corikaIOreader options
corikaIOreaderOpt="-queff 1 -cors $coriska_input -seed $seed -grisu stdout -cfg $telescope_array"
if [ "$atm_extintion" != "" ]
then
  corikaIOreaderOpt="$corikaIOreaderOpt -abs corsika -absfile $atm_extintion"
fi

# - Run
echo " -- Running:"
echo "  ../corsikaSimulationTools/corsikaIOreader $corikaIOreaderOpt | ./grOptics -p $pilot &> $tmp_folder/log"
echo ""
../corsikaSimulationTools/corsikaIOreader $corikaIOreaderOpt | ./grOptics -p "$pilot" &> "$tmp_folder/log"
corsikaIOreaderStatus=${PIPESTATUS[0]}
grOpticsStatus=$?
echo ""

# Move results in the proper folder
cd "$tmp_folder" # This to create a local path in the tarball
tar -czf cfg.tar.gz cfg
cd "$gropticsPath" # This to go back to the previous folder before deleting the tmp
rm -rf "$tmp_cfg"
gzip "$tmp_folder/log"

if [ $corsikaIOreaderStatus -ne 0 ]
then
    echo " *** BAD EXECUTION: something went wrong with corsikaIOreader ***"
    echo "      (all in \"$badcIOr\"... check the log!)"
    mv "$tmp_folder"/* "$badcIOr"
elif [ $grOpticsStatus -ne 0 ]
then
    echo " *** BAD EXECUTION: something went wrong with grOptics ***"
    echo "      (all in \"$badgroptics\"... check the log!)"
    mv "$tmp_folder"/* "$badgroptics"
else
    echo " ==> GOOD EXECUTION (all in \"$output\") <=="
    mv "$tmp_folder"/* "$output"
fi

# Cleaning
rm -rf "$tmp_folder"

# Ending info
echo ""
echo "#-------------------------------------------"
echo "#"
echo "#            END (run ${RUN})"
echo "#"
echo "#  grOptics on ${USER}@${HOSTNAME}"
echo "#  time: `date`"
echo "#"
echo "#-------------------------------------------"
echo ""
