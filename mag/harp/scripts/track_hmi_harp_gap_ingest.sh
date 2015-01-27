#!/bin/sh -f
# (line below is for qsub)
#$ -S /bin/sh
#
# track_hmi_harp_gap_ingest: HARP gap-filler -- ingest already-staged data
#
# sh driver to ingest data and make movies as part of the HARP gap-filler
#   ...to be run after the companion script, track_hmi_harp_gap_stage, has run
#
# Usage:
#   track_hmi_harp_gap_ingest.sh [-n] [-d] harp harp_log dest_dir
#
# where arguments are:
#   harp      -- a data series from the HARP jsd which holds the patched 
#                HARP information
#   harp_log  -- harp log data series, for the ingestion log
#   dest_dir  -- the full path of the HARP staging area, no Tracks/jsoc part
# 
# and options are:
#   -d,  flag to use developer (~turmon) path, as opposed to ~jsoc
#
# * Script is test-run by the developer, because its results should be 
# checked by hand before the final ingestion into DRMS (hmi.Mharp_720s) by
# the JSOC operator.
# * So, use a dummy HARP series to test, and hmi.Mharp_720s to apply
# the tested fix.  Check the movie and the dummy series to evaluate the fix.
# * This script contains commands to ingest harps contained in files under
# the staging area (in dest_dir/Tracks/jsoc) into a data and log series
# to fill in the gap between (inclusive)
#   t0 = $t0
#   t1 = $t1
# left between regular runs $run_no and $run_nu.
#
# turmon 2013-mar

set -e

progname=`basename $0`
# under SGE, $0 is set to a nonsense name, so use our best guess instead
if [ `expr "$progname" : '.*track.*'` == 0 ]; then progname=track_hmi_harp_gap_ingest; fi
USAGE="Usage: $progname [-d] harp harp_log dest_dir"

# echo the arguments, for repeatability
echo "${progname}: Invoked as:" "$0" "$@"

# get options
developer_path=0
while getopts "d" o; do
    case "$o" in
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
harp="$1"
harp_log="$2"
dest_dir="$3"

###
### Validate arguments
###

# validate given series
if ! describe_series -j "$harp" > /dev/null ; then 
    echo "${progname}: given harp series \`$harp' is not in DRMS"
    echo "$USAGE" 1>&2
    exit 1
fi
if ! describe_series -j "$harp_log" > /dev/null ; then 
    echo "${progname}: given harp log series \`$harp_log' is not in DRMS"
    echo "$USAGE" 1>&2
    exit 1
fi

# validate dest_dir
if [ ! -d "$dest_dir" ]; then
    echo "${progname}: given dest_dir \`$dest_dir' is not readable"
    echo "$USAGE" 1>&2
    exit 1
fi

# staging area
stage=$dest_dir/Tracks/jsoc
log_save=$stage/Gap-Run
pfile=$log_save/track-param.txt

# more validation
if [ ! -d "$stage" ]; then
    echo "${progname}: staging area \`$stage' is not readable"
    echo "$USAGE" 1>&2
    exit 1
fi
if [ ! -d "$log_save" ]; then
    echo "${progname}: given saved-log area \`$log_save' is not readable"
    echo "$USAGE" 1>&2
    exit 1
fi
if [ ! -r "$pfile" ]; then
    echo "${progname}: saved-log area does not have a param file \`$pfile'"
    echo "$USAGE" 1>&2
    exit 1
fi


###
### Process arguments
###

# track-and-ingest options
TI_opts=""
# harp-movie options
HM_opts=""
if [ $developer_path -eq 1 ]; then
  TI_opts="$TI_opts -d"
  HM_opts="$HM_opts -d"
fi

# basic series
mag=hmi.M_720s
mask=hmi.Marmask_720s
harp_ds=hmi.Mharp_720s

# load t00, t11, etc., etc.
param_file=$dest_dir/gap-params.dat
if [ ! -r "$param_file" ]; then
  echo "${progname}: No parameter file ($param_file). Exiting."
  exit 1
fi
source $param_file

echo "${progname}: Gap is between $run_no and $run_nu."

# additional parameters for the movie interval
t00_index=`index_convert ds=$mag T_REC=$t00`
t11_index=`index_convert ds=$mag T_REC=$t11`

# t00day = t0 - 1 day
t00day_index=$(( $t00_index - 120 ))
t00day=`index_convert ds=$mag T_REC_index=$t00day_index`

# t11day = t11 + 2 days
#   ensure t11day is late enough so a HARP is created between t11 and t11day;
#   this enables us to check that the ROI numbers are synched
t11day_index=$(( $t11_index + 240 ))
t11day=`index_convert ds=$mag T_REC_index=$t11day_index`

cat <<EOF
${progname}: Gap summary:
  Ingesting from $dest_dir 
  Ingesting to ($harp, $harp_log)
  Run $run_no end: $t00
  Gap from $t0 - $t1 ($gap_length slots)
  Run $run_nu begin: $t11
  Diagnostic movie from $t00day - $t11day
EOF

#
# ========================================================================
#
# Produce a standard status report

keys=HARPNUM,T_REC,LATDTMIN,LONDTMIN,LATDTMAX,LONDTMAX,NPIX,DATE

fillfile=$dest_dir/gap-fill-status.txt
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

# fix directory is specific to this destination series
#  (runs with dummy series and harp series should not collide)
fix_dir=$dest_dir/fix/$harp
rm   -rf $fix_dir
mkdir -p $fix_dir

# ingest over T_REC range spanning gap
#   these parameter values are the same as in the standard usage for definitive 
#   HARPs in track_and_ingest_mharp.sh, except that we use the big list in
#   track-final.txt to ingest from
echo "${progname}: Ingesting into $harp to fill gap"
ingest_mharp root=$stage out=${harp} mask=${mask} lists=track-final.txt,track-pending.txt tmin=3 verb=1 href=$harp "t0cut=$t0" "t1cut=$t1" | tee $fix_dir/gap-fix-ingest.log

# ingest log as well
#   note, log ingest uses preserved files *below* $stage, because the 
#   second track_and_ingest run 
echo "${progname}: Ingesting logs into $harp_log"
ingest_mharp_log root=$log_save out=$harp_log trec=$t0 | tee $fix_dir/gap-fix-ingest-log.log

# re-make the harp movie across the gap to double-check the fix
echo "${progname}: Making HARP movie spanning former gap in $harp"
track_hmi_harp_movie_driver.sh -m -d "${mask}[${t00day}-${t11day}]" $harp $fix_dir | tee $fix_dir/gap-fix-movie.log

runt1=`date +%s`
runtdiff=$(( $runt1 - $runt0 ))
echo "${progname}: Ingestion and movie took $runtdiff seconds."

echo "${progname}: Done."

