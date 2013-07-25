#!/bin/sh
# (line below is for qsub)
#$ -S /bin/sh
# 
# sh driver to run the Matlab tracker on HMI data from JSOC
#   Note: This is not the top-level driver.  Ingestion is done in a 
#   different script.  See track_and_ingest_mharp.sh for the
#   top-level driver integrating tracking (this script) and 
#   ingestion (loading tracker results into JSOC).
#
# If this file has _stable in its name, it is the production driver script,
#   and is intended to operate from the cvs/JSOC tree.
# If not, it is the in-development driver, and may reside in a developer's 
#   personal directory.
# The stable one is intended for production.
#
# usage:
#   track_hmi_production_driver_stable.sh ...
#       [ -i N ] [ -r R ] [ -p PARS ] [-m] [-g N] [-d] mask_series mag_series dest_dir
#
# where:
#   mask_series is a mask data series from JSOC (includes T_REC range)
#   mag_series is a magnetogram data series from JSOC (no T_REC range)
#   dest_dir is a destination directory for file and state input and output
#
# and:
#   -i introduces an optional argument indicating an initial run with 
#      new tracks numbered starting from N.  Otherwise, a run that picks
#      up where an earlier one left off, according to a checkpoint file
#      in `dest_dir', is assumed.
#      IF GIVEN, `-i N' CLEARS ALL EXISTING RESULTS.
#   -r introduces an optional argument telling how many chips of history
#      to retain.  If not given, all chips are retained.  Otherwise, only
#      an integer `R' chips of history are retained.  The default is
#      to retain all history, which is equivalent to giving R = inf.
#   -p introduces an argument allowing modification of the default 
#      tracker parameters.  This is intended for testing purposes.
#      PARS is a series of semicolon-terminated assignments using field
#      names from `tracker_params' in track_hmi_production.m, for example:
#      -p "tau=0.19;final_time=0.5;"  
#      For experts only.
#   -m is a flag to generate a movie from the results (default is no movie).
#   -g says to use gap-filling mode for new tracks, and supplies a 
#      required argument N.  In this mode, new tracks will be checked
#      for overlap with future tracks; if so, the track ID will be assigned
#      to match.  If there is no overlap, new tracks will get an ID that 
#      was left unused due to earlier merges.  This means new tracks 
#      during the run will (mostly) not affect numbering of later tracks.  
#      The argument, typically 0, is a number to add to the track ID
#      counter (ROI_rgn.Ntot) at the end of the run, in case more help 
#      is needed to synch the gap-fill count with the desired count.
#   -d is a flag which signals to use the developer MATLABPATH (~turmon) rather
#      than the regular production path (~jsoc).
#
# * Both masks and corresponding magnetograms are needed.  The magnetogram 
# series used is typically hmi.M_720s, but you can use hmi.M_720s_nrt.  
# You DO NOT give a T_REC range to the magnetograms, only the series name.
# The mask_series is typically of the form: 
#     any_series[any_range] 
# for example,
#     hmi.Marmask_720s[2011.01.01/7d]
# * If dest_dir does not exist, it will be created.
# * As mentioned above, by default, results of the run are added to the 
# existing dest_dir.  If you give -i N, the dest_dir is wiped clear and 
# tracking is restarted.
#
# Usage:
#   If tracking was already started, and an earlier run ended on 
# Sept 19 at midnight:
#
#  track_hmi_production_driver_stable.sh 'hmi.Marmask_720s[2010.09.20_TAI/24h]' hmi.M_720s /tmp/sept
#
# Michael Turmon, JPL, December 2010, June 2011, July 2013

set -e

progname=`basename $0`
# under SGE, $0 is set to a nonsense name, so use our best guess instead
if [ `expr "$progname" : '.*track.*'` == 0 ]; then progname=track_hmi_production_driver; fi
USAGE="Usage: $progname [-i N] [-r R] [-p PARS] [-m] [-g N] [-d] mask_series mag_series dest_dir"
# echo the arguments, for repeatability
echo "${progname}: Invoked as:" "$0" "$@"

# Trap error exits through this routine
#   especially important for shell-command exits (trap ... ERR below)
#   usage: die LINENO MESSAGE
#     (arguments optional)
function die() {
    echo "${progname}: Fatal error near line ${1:- \?}: ${2:-(no message)}" 1>&2
    echo "${progname}: Exiting with error." 1>&2
    echo "${progname}: Exiting with error."
    exit 1
}

# Error message selection upon error exit (NB: ordinary exit does not trap)
cmd="shell code"
trap 'die $LINENO "Running $cmd"'     ERR
trap 'die $LINENO "Received SIGHUP"'  SIGHUP
trap 'die $LINENO "Received SIGINT"'  SIGINT
trap 'die $LINENO "Received SIGTERM"' SIGTERM

# get options
first_track=0
retain_history=inf
make_movie=0
gap_fill=nan
developer_path=0
matlab_ver=r2010b
params=""
while getopts "dmg:p:i:r:v:" o; do
    case "$o" in
	i)    first_track="$OPTARG";;
	r)    retain_history="$OPTARG";;
	p)    params="$OPTARG";;
	g)    gap_fill="$OPTARG";;
	m)    make_movie=1;;
	d)    developer_path=1;;
	v)    matlab_ver="$OPTARG";;
	[?])  echo "$USAGE" 1>&2
	      exit 2;;
    esac
done
shift `expr $OPTIND - 1`

# get arguments
if [ "$#" -ne 3 ]; then
    echo "$USAGE" 1>&2
    die $LINENO "Got $# args, expected 3"
fi

echo "${progname}: Beginning."
# ok to get args -- note, we add {mask} on to the mask series
mask_series="$1{mask}"
mag_series="$2"
dest_dir="$3"


# validate first_track is an int
if [[ $first_track =~ ^[0-9]+$ ]]; then 
    : 
else 
    die $LINENO "bad first track integer (got $first_track)"
fi

# validate retain_history is int or inf
if [[ $retain_history =~ ^[0-9]+$ || $retain_history =~ ^inf$ ]]; then 
    : 
else 
    die $LINENO "bad retain history integer (got $retain_history)"
fi

# validate gap_fill is integer (negative OK) or nan
if [[ $gap_fill =~ ^[+-]?[0-9]+$ || $gap_fill =~ ^nan$ ]]; then 
    : 
else 
    die $LINENO "bad gap-fill integer (got $gap_fill)"
fi

# there is no sh version of .setJSOCenv, so we use this to get JSOC_MACHINE
if [ -z "$JSOC_MACHINE" ]; then
  JSOC_MACHINE=`csh -f -c 'source /home/jsoc/.setJSOCenv && echo $JSOC_MACHINE'`
fi
# set up matlab path
if [ "$developer_path" -eq 1 ]; then
    # (this is the code path when -d is given)
    # the tree under ~turmon/matlab/mfile -- development versions
    mtHome=`echo ~turmon`
    MATLABPATH=`find $mtHome/matlab/mfile -maxdepth 1 -type d | paste -d : -s -`
else
    # (this is the default code path)
    # the production path
    rootD=/home/jsoc/cvs/Development/JSOC
    rootDbin=$rootD/_$JSOC_MACHINE/proj/mag/harp/libs/matlab/mfile-mex
    rootDsrc=$rootD/proj/mag/harp/libs/matlab/mfile-plain
    if [ ! -d "$rootDbin" ]; then
	die $LINENO "Bad path to mex shared libs (got $rootDbin)"
    fi
    if [ ! -d "$rootDsrc" ]; then
	die $LINENO "Bad path to matlab scripts (got $rootDsrc)"
    fi
    # all matlab m-files and shared libraries (.m and .mexFOO) are in these dirs
    MATLABPATH="$rootDbin/fits:$rootDbin/assignment:$rootDbin/hmi-mask-patch:$rootDbin/standalone:$rootDsrc"
fi
 
export MATLABPATH

#
# Matlab version
#
# The old JSOC matlab is r14sp3, which is certified to work with glibc 2.3.4, 
# which is the kernel version on n02.
# The new JSOC matlab is R2010b, which wants glibc 2.7, but appears to work
# OK with our older kernel.  R2010b does a check at startup which
# issues a prompt to continue or abort.  The prompt can be disabled with
# the OSCHECK_ENFORCE_LIMITS line.  It still gives a stern warning.
if [ "$matlab_ver" = "r14" ]; then
    ## r14sp3 = version 7.1 from 2005
    matlab_binary=/usr/local/matlabr14sp3/bin/matlab
elif [ "$matlab_ver" = "r2010b" ]; then
    # R2010b = version 7.11 from 2010
    matlab_binary=/usr/local/MATLAB/R2010b/bin/matlab
else
    die $LINENO "Only recognize matlab version r14 or r2010b"
fi
if [ ! -x "$matlab_binary" ]; then
    die $LINENO "Bad path to matlab binary (got $matlab_binary)"
fi
# needed for R2010b (see above); no harm otherwise
OSCHECK_ENFORCE_LIMITS=0; export OSCHECK_ENFORCE_LIMITS

# the basic matlab call
# notes:
#  the try/catch is good matlab to catch trouble 
#    (but, "catch ME" is newer than r14sp3)
#    there is no way to return error status from matlab,
#    so we resort to looking for lines with ERROR
#    Thus, to play nice, the script must start its error messages this way.
tmpnam="/tmp/$progname.$$.err"
rm -f "$tmpnam"
t0=`date +%s`
cmd=$matlab_binary
$matlab_binary -nodesktop -nosplash -nodisplay -logfile "$tmpnam" <<EOF
% use outer try/catch to catch *all* exceptions
try,
  mask_series='$mask_series';
  mag_series='$mag_series';
  dest_dir='$dest_dir';
  first_track=$first_track;
  retain_history=$retain_history;
  gap_fill=$gap_fill;
  make_movie=$make_movie;
  params='$params';
  track_hmi_production;
  fprintf(1, '\\nEXIT_STATUS_OK: clean exit from Matlab.\\n');
catch ME,
  % ensure these errors start on a new line
  fprintf(1, '\\nERROR: Report follows.\\n%s\\n', getReport(ME));
  fprintf(2, '\\nERROR: Report follows.\\n%s\\n', getReport(ME));
  rethrow(ME);
end;
quit;
EOF
cmd="shell code"
t1=`date +%s`
tdiff=$(( $t1 - $t0 ))
echo "${progname}: Matlab call took $tdiff seconds."

# put the log file back
tmpnam2="$dest_dir/Tracks/jsoc/track-latest.log"
mkdir -p `dirname $tmpnam2` # dir present if cmd ran ok, maybe not if error
rm -f "$tmpnam2"
mv -f "$tmpnam" "$tmpnam2"

# generate exit status
if grep -q "^ERROR" "$tmpnam2" ; then
    # grep above matched error sentinel
    die $LINENO "Error caught in Matlab."
elif grep -q "^EXIT_STATUS_OK" "$tmpnam2" ; then
    # grep above matched OK sentinel
    echo "${progname}: Exiting OK."
    echo "${progname}: Done."
    exit 0
else
    # no OK indicator found
    die $LINENO "Error running Matlab, caught by default." 
fi

