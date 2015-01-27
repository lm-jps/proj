#!/bin/sh -f
# (line below is for qsub)
#$ -S /bin/sh
#
# track_hmi_harp_gap_stage: HARP gap-filler -- stage and test data for fill
#
# sh driver to stage data and make movies as part of the HARP gap-filler
#
# Usage:
#   track_hmi_harp_gap_stage.sh [-n] [-d] [-g N] run_number dummy_harp dest_dir
#
# where arguments are:
#   run_number -- TKP_RUNC field of hmi.Mharp_log_720s just before the gap
#   dummy_harp -- a data series from the HARP jsd which holds the patched 
#                 HARP information
#   dest_dir   -- the full path of the HARP staging area, no Tracks/jsoc part
# 
# and options are:
#   -n,  flags dry run; the expensive/destructive operations are not done
#   -d,  flag to use developer (~turmon) path, as opposed to ~jsoc
#   -g,  introduces integer to add to region count (positive or negative)
#
# * Script is run by the developer, because its results must be 
# checked by hand before ingestion.  A test run is made to stage the data 
# and check results of the gap-filling.  Then the JSOC operator can
# use the staged files to ingest results into DRMS (hmi.Mharp_720s).
# * This script writes another tailored script to dest_dir/fix. The second
# script ingests the staged data (in dest_dir/Tracks/jsoc) into a series
# named by its argument, then making a summary movie from the series.  
# This should be used once on a dummy series, and again on hmi.Mharp_720s.
#
# turmon 2013-feb

set -e

progname=`basename $0`
# under SGE, $0 is set to a nonsense name, so use our best guess instead
if [ `expr "$progname" : '.*track.*'` == 0 ]; then progname=track_hmi_harp_gap_fill; fi
USAGE="Usage: $progname [-n] [-d] [-g n] run_number dummy_harp dest_dir"

# echo the arguments, for repeatability
echo "${progname}: Invoked as:" "$0" "$@"

# get options
noop=0
developer_path=0
gap_fill=0
while getopts "ndg:" o; do
    case "$o" in
	n)    noop=1;;
	g)    gap_fill="$OPTARG";;
	d)    developer_path=1;;
	[?])  echo "$USAGE" 1>&2
	      exit 2;;
    esac
done
shift `expr $OPTIND - 1`

# get arguments
if [ "$#" -ne 3 ]; then
    echo "$USAGE" 1>&2
    exit 2
fi

###############################################################
## PATH STUFF BEGINS
# (need developer_path set for this)
# SGE/OpenMP setup
#SGE_ROOT=/SGE;     export SGE_ROOT
OMP_NUM_THREADS=1; export OMP_NUM_THREADS
KMP_BLOCKTIME=10;  export KMP_BLOCKTIME

PATH="${PATH}:$SGE_ROOT/bin/$SGE_ARCH"

#uname=`uname -m`
#if [ "X$uname" = Xi686 ]; then
#    PATH="${PATH}:$SGE_ROOT/bin/lx24-x86"
#elif [ "X$uname" = Xx86_64 ]; then
#    PATH="${PATH}:$SGE_ROOT/bin/lx24-amd64"
#elif [ "X$uname" = Xx86_64 ]; then
#    PATH="${PATH}:$SGE_ROOT/bin/lx24-ia64"
#else
#    echo "Could not find system '$uname' to set up SGE" 1>&2; exit 2
#fi
# ensure track_hmi_production_driver_stable.sh is in PATH
if [ "$developer_path" -eq 1 ]; then
  # put script dir of developer into path
  # NOTE: changed to .../proj/mag/harp/scripts
  PATH="${PATH}:$HOME/cvs/JSOC/proj/mag/harp/scripts"
else
  # put script dir of JSOC into path
  PATH="${PATH}:/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts"
fi
# for show_info, describe_series
# there is no sh version of .setJSOCenv, so we use this to get JSOC_MACHINE
jsoc_mach=`csh -f -c 'source /home/jsoc/.setJSOCenv && echo $JSOC_MACHINE'`
PATH="${PATH}:$HOME/cvs/JSOC/bin/$jsoc_mach"
## PATH STUFF ENDS
###############################################################

echo "${progname}: Beginning."
runt0=`date +%s`

# ok to get args
run_no="$1"
dummy_harp="$2"
dest_dir="$3"

#
# Validate arguments
#

# validate first_track is an int
if [[ $run_no =~ ^[0-9]+$ ]]; then 
    : 
else 
    echo "${progname}: bad run number integer (got $run_no)"
    echo "$USAGE" 1>&2
    exit 2
fi

# validate given series
if ! describe_series -j "$dummy_harp" > /dev/null ; then 
    echo "${progname}: dummy harp series \`$dummy_harp' is not in DRMS"
    echo "$USAGE" 1>&2
    exit 2
fi

# try to create it now
mkdir -p $dest_dir

# ========================================================================

#
# Process arguments
#

# handle no-op mode
if [ $noop -eq 1 ]; then
  pre="echo DRYRUN:"
  descrip="dry run"
else
  pre=
  descrip="regular run"
fi

# track-and-ingest options
TI_opts=""
# harp-movie options
HM_opts=""
if [ $developer_path -eq 1 ]; then
  TI_opts="$TI_opts -d"
  HM_opts="$HM_opts -d"
fi

# directory where logs produced by this script go
out_dir=$dest_dir/stage
$pre rm -r $out_dir
mkdir -p $out_dir

# subsequent run number
run_nu=`expr $run_no + 1`

echo "${progname}: Run type is *${descrip}*."
echo "${progname}: Gap is between $run_no and $run_nu."

# ========================================================================

# basic series
mag=hmi.M_720s
mask=hmi.Marmask_720s
harp=hmi.Mharp_720s
harp_log=hmi.Mharp_log_720s
dummy_log=su_turmon.Mharp_logv5_720s

# t00 = last frame seen before gap
t00=`show_info -q ds="${harp_log}[][?TKP_RUNC=${run_no}?]" key=T_LAST`
t00_index=`index_convert ds=$mag T_REC=$t00`
# t0 = t00 + delta = first frame within the gap
t0_index=`expr $t00_index + 1`
t0=`index_convert ds=$mag T_REC_index=$t0_index`

# t11 = first frame seen after gap
t11=`show_info -q ds="${harp_log}[][?TKP_RUNC=${run_nu}?]" key=T_FRST | tail -1`
t11_index=`index_convert ds=$mag T_REC=$t11`
# t1 = t11 - delta = final frame within the gap
t1_index=`expr $t11_index - 1`
t1=`index_convert ds=$mag T_REC_index=$t1_index`

# tmax: final end-time of the gap (last detection)
tmax=`show_info -q ds="${harp}[][${t00}-${t11}]" key=T_LAST1 | sort | tail -1`

# gap length in slots
gap_length=$(( $t1_index - $t0_index + 1 ))

# Find existing checkpoint file
checkpoint=`show_info -q -P ds="$harp_log[][?TKP_RUNC=${run_no}?]" seg=Checkpoint`
if [ ! -r "$checkpoint" ]; then
  echo "${progname}: No checkpoint (offline?)"
  exit 1
fi

# ========================================================================

# Make a parameter file with these times
param_file=$dest_dir/gap-params.dat
cat >$param_file <<EOF
# gap times 
#   by $USER
#   on `date`
#   from $progname
t0=$t0
t1=$t1
t00=$t00
t11=$t11
run_no=$run_no
run_nu=$run_nu
gap_length=$gap_length
EOF

# ========================================================================

# Summarize gap-filling run we're about to do
runfile=$dest_dir/gap-interval.txt
cat <<EOF | tee $runfile
${progname}: Gap summary:
  Run $run_no end: $t00 ($checkpoint)
  Gap from $t0 - $t1 ($gap_length slots)
  Run $run_nu begin: $t11
  Subsequent scan out to $tmax
  Results: $dest_dir and $dummy_harp
EOF

# ========================================================================

# Create a staging area and place the last (before-gap) checkpoint there:
stage=$dest_dir/Tracks/jsoc
echo "${progname}: Putting checkpoint in place"
# note, regenerating is HARD!
$pre rm -rf $stage
mkdir -p $stage
# superior to cp if the file is already there
rsync -t $checkpoint $stage/track-prior.mat

# Run track-and-ingest in gap-filling mode between t0 and t1:
#   (dummy HARP series to hold the ingested data)
echo "${progname}: Tracking across gap"
logfile="$out_dir/gap-track-within.log" && [ $noop -eq 1 ] && logfile=/dev/null
$pre track_and_ingest_mharp.sh $TI_opts -g "$gap_fill" "$dest_dir" "$mask[${t0}-${t1}]" $dummy_harp $dummy_log | tee $logfile 

# preserve the logs and checkpoint of the within-gap run for later ingest
log_save=$stage/Gap-Run
$pre rm -rf $log_save
mkdir -p $log_save
$pre rsync -t $stage/*.* $log_save/
# negate the run_number parameter of the param file that we will ingest
pfile=track-param.txt
if [ $noop -ne 1 ]; then
  awk '/^run_number:/ {print $1, -1*$2; next};{print}' $stage/$pfile > $log_save/$pfile
fi

# Run track-and-ingest from t11 out to tmax, so that all bitmaps are combined as necessary
echo "${progname}: Tracking after gap"
logfile="$out_dir/gap-track-after.log" && [ $noop -eq 1 ] && logfile=/dev/null
$pre track_and_ingest_mharp.sh $TI_opts -a $dest_dir "$mask[${t11}-${tmax}]" $dummy_harp $dummy_log | tee $logfile 

# Make a movie of the new frames between t0-tmax:
echo "${progname}: Making HARP movie"
logfile="$out_dir/gap-movie.log" && [ $noop -eq 1 ] && logfile=/dev/null
$pre track_hmi_harp_movie_driver.sh -m $HM_opts "$mask[${t0}-${tmax}]" $dummy_harp $out_dir | tee $logfile

echo "${progname}: Main processing done"

# ========================================================================

# Produce a standard status report
#   (only if this was a real run)

keys=HARPNUM,T_REC,LATDTMIN,LONDTMIN,LATDTMAX,LONDTMAX,NPIX,DATE

fillfile="$out_dir/gap-fill-status.txt" && [ $noop -eq 1 ] && fillfile=/dev/null
cat >$fillfile <<EOF
Status of gap-fill between $run_no and $run_nu (run on `date`)

Gap is between (inclusive):
  t0 = $t0
  t1 = $t1

Gap was captured in series $dummy_harp

The following HARP lists should align:

== New HARPs created during the gap ${t0}-${t1}:
`show_info ds="$dummy_harp[][${t0}-${t1}][?T_REC=T_FRST1?]" key="$keys"`

== New HARPs created after the gap:
`show_info ds="$dummy_harp[][${t11}/2d][?T_REC=T_FRST1?]" key="$keys"`

== Existing HARPs ($harp) created after the gap:
`show_info ds="$harp[][${t11}/2d][?T_REC=T_FRST1?]" key="$keys"`

EOF

# ========================================================================

runt1=`date +%s`
runtdiff=$(( $runt1 - $runt0 ))
echo "${progname}: Gathering info took $runtdiff seconds."

echo "${progname}: Done."

