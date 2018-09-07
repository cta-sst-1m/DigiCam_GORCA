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
    echo "Launch a biascurve simulation from CARE"
    echo "  usage: $0 <production_name> <run> [ <output_name> <defaults> ] [ <cfg> <telescopes> ] [ <extras>... ]"
    echo
    echo "<production_name>   Name of the output (and base of the tmp) folder"
    echo "<run>               Run ID (used also as seed)"
    echo
    echo "<output_name>       Name of the output files (run name!) [Default: careRun<run>]"
    echo "<defaults>          Default of the (following) settings [Default: $runPath/defaults.sh]"
    echo
    echo "<cfg>               CARE configuration file [Default defined in <defaults>]"
    echo "<telescopes>        Telescope array configuration file [Default defined in <defaults>]"
    echo
    echo "<extras>            All possible available arguments in CameraAndReadout executable"
    echo
    echo "Let's try to estimate the bias curve!"
    echo
    exit 1
elif [ $# -lt 2 ]
then
    echo
    echo "ERROR: at least 2 arguments (<name> <run>) must be given!"
    echo
    exit 1
fi

# Base setting (ARGS 1 and 2: NAME and RUN)
NAME="$1"
RUN="$2"
run=`printf "%06d" $RUN`
seed=$((RUN+1))
biascurve=1
loopoverevents=0

# Set the output name
output_name="${3-careRun$run}"
#if [ "$output_name" != "" ]
#then
#    output_name="${output_name}."
#fi
#output_name="${output_name}careRun$run"

# Default settings
if [ "$4" == "" ]
then
    settings="$runPath/defaults.sh"
else
    settings="$4"
fi
source "$settings"

# Other configs
if [ "$5" == "" ]
then
    CFG="$care_default"
else
    CFG="$5"
fi
if [ "$6" == "" ]
then
    telescopes="$telescopes_default"
else
    telescopes="$6"
fi
extras="${@:7}"
if [ "$extras" != "" ]
then
    extras="$extras"" "
fi

# Set the temporary folder/configs
tmp_folder="${tmp_dir}/${NAME}_${PBS_JOBID-`date +%d%m%y%H%M%S%N`}"

cfgFolder="$tmp_folder/cfg"
mkdir -p "$cfgFolder"
tmp_cfg="$cfgFolder/care.cfg"
cp "${pulse}HighGain.dat" "$cfgFolder/pulseHighGain.dat"
cp "${pulse}LowGain.dat" "$cfgFolder/pulseLowGain.dat"
cp "$telescopes" "$cfgFolder/telescopes.dat"
telescopes="$cfgFolder/telescopes.dat"
cp "$qe" "$cfgFolder/qe.dat"
qe="$cfgFolder/qe.dat"
if [ "$patches" == "" ]
then
    patches="# No patches!"
    patches_sed="s=@patches@=${patches}=g"
else
    cp "$patches" "$cfgFolder/patches.dat"
    patches="$cfgFolder/patches.dat"
    patches_sed="/@patches@/ {
r ${patches}
d
}"
fi
cp "$camera" "$cfgFolder/camera.dat"
camera="$cfgFolder/camera.dat"

# Autoset from config files
ntel=`grep '* TLCFG' "$telescopes" | wc -l`
npixels=`grep '* PMPIX' "$camera" | wc -l`
ngroups=`grep '* GROUP' "$camera" | wc -l`

# Set output folder and tree structure
outbase="${out_dir}/$NAME/BiasCurve"
logs="$outbase/Logs"
output="$outbase/Roots"
badcare="$outbase/Bads"
mkdir -p "$logs"
mkdir -p "$output"
mkdir -p "$badcare"

cd "$carePath"

# Starting info
echo "#-------------------------------------------"
echo "#"
echo "#            START (run ${RUN})"
echo "#"
echo "#  CameraAndReadout on ${USER}@${HOSTNAME}"
echo "#    (BiasCurve simulation)"
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
echo "#"
echo "#           ...cwd: '`pwd`'"
echo "#"
echo "#-------------------------------------------"
echo ""
echo ""

# Build the local configs
today=`date +%F`
SED="s=@ntel@=${ntel}=g;s=@npixels@=${npixels}=g;s=@ngroups@=${ngroups}=g;s=@cfgFolder@=${cfgFolder}=g;s=@loopoverevents@=${loopoverevents}=g;s=@biascurve@=${biascurve}=g;s=@today@=${today}=g;/@telescopes@/ {
r ${telescopes}
d
};/@qe@/ {
r ${qe}
d
};${patches_sed};/@camera@/ {
r ${camera}
d
}"
cat <<EOF 1>&2
 -- Local setting of the configuration:
  - "$CFG" to "$tmp_cfg":
      sed -e "$SED" "$CFG" > "$tmp_cfg"

EOF
sed -e "$SED" "$CFG" > "$tmp_cfg"

# Make exe file with arguement parsing
ARGS=(-s "$seed" --configfile "$tmp_cfg" --outputfile "$tmp_folder/BC" "${extras[@]}")
echo "./CameraAndReadout ${ARGS[@]} &> $tmp_folder/log" > "$tmp_folder/exe"
chmod +x "$tmp_folder/exe"

# Run
echo " -- Running:"
echo "  "`cat "$tmp_folder/exe"`
#echo "  ./CameraAndReadout -s $seed --configfile $tmp_cfg --outputfile $tmp_folder/CARE --inputfile $groptics_input ${extras[@]} &> $tmp_folder/log"
echo ""
echo " -- Custom configuration of the configuration file: ${extra_fields[@]}"
echo ""

"$tmp_folder/exe"
#./CameraAndReadout "${extras}" -s "$seed" --configfile "$tmp_cfg" --outputfile "$tmp_folder/BC" &> "$tmp_folder/log"

careStatus=$?

cd "$tmp_folder" # This to create a local path in the tarball
tar -czf cfg.tar.gz cfg
cd "$carePath" # This to go back to the previous folder before deleting the tmp
rm -rf "$tmp_folder/cfg"
gzip "$tmp_folder/log"

if [ $careStatus -ne 0 ]
then
    echo " *** BAD EXECUTION: something went wrong with biascurve ***"
    echo "      (main dir: \"$badcare\"... check the $badcare/$output_name.log.gz!)"
    mv "$tmp_folder/BC.root" "$badcare/$output_name.root"
    mv "$tmp_folder/log.gz" "$badcare/$output_name.log.gz"
    mv "$tmp_folder/cfg.tar.gz" "$badcare/$output_name.cfg.tar.gz"
else
    echo " ==> GOOD EXECUTION (main dir: \"$outbase\") <=="
    mv "$tmp_folder/BC.root" "$output/$output_name.root"
    mv "$tmp_folder/log.gz" "$logs/$output_name.log.gz"
    mv "$tmp_folder/cfg.tar.gz" "$logs/$output_name.cfg.tar.gz"
fi

# Cleaning
rm -rf "$tmp_folder"

# Ending info
echo ""
echo "#-------------------------------------------"
echo "#"
echo "#            END (run ${RUN})"
echo "#"
echo "#  CameraAndReadout on ${USER}@${HOSTNAME}"
echo "#    (CARE simulation)"
echo "#  time: `date`"
echo "#"
echo "#-------------------------------------------"
echo ""
