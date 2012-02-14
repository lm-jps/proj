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
#       [ -i N ] [ -r R ] [ -p PARS ] [-m] [-d] mask_series mag_series dest_dir
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
#   If tracking was already started, and an earlier run ended on Sept 19 at midnight:
#
#  track_hmi_production_driver_stable.sh 'hmi.Marmask_720s[2010.09.20_TAI/24h]' hmi.M_720s /tmp/sept
#
# Michael Turmon, JPL, December 2010, June 2011


progname=`basename $0`
USAGE="Usage: $progname [ -i N ] [ -r R ] [ -p PARS ] [ -m ] [ -d ] mask_series mag_series dest_dir"

# get options
first_track=0
retain_history=inf
make_movie=0
developer_path=0
matlab_ver=r2010b
params=""
while getopts "dmp:i:r:v:" o; do
    case "$o" in
	i)    first_track="$OPTARG";;
	r)    retain_history="$OPTARG";;
	p)    params="$OPTARG";;
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
    exit 2
fi
# ok to get them -- note, we add {mask} on to the mask series
mask_series="$1{mask}"
mag_series="$2"
dest_dir="$3"

# validate first_track is an int
if [[ $first_track =~ ^[0-9]+$ ]]; then 
    : 
else 
    echo "$progname: bad first track integer (got $first_track)"
    echo "$USAGE" 1>&2
    exit 2
fi

# validate retain_history is int or inf
if [[ $retain_history =~ ^[0-9]+$ || $retain_history =~ ^inf$ ]]; then 
    : 
else 
    echo "$progname: bad retain history integer (got $retain_history)"
    echo "$USAGE" 1>&2
    exit 2
fi

# set up matlab path
if [ "$developer_path" -eq 1 ]; then
    # (this is the code path when -d is given)
    # the tree under ~turmon/matlab/mfile -- development versions
    mtHome=`echo ~turmon`
    MATLABPATH=`find $mtHome/matlab/mfile -maxdepth 1 -type d | paste -d : -s -`
else
    # (this is the default code path)
    # believe this is the right production path
    rootD=/home/jsoc/cvs/Development/JSOC
    # for a time, was using arta's path
    # rootD=/home/arta/jsoctrees/JSOC
    rootDbin=$rootD/_$JSOC_MACHINE/proj/mag/harp/libs/matlab/mfile-mex
    rootDsrc=$rootD/proj/mag/harp/libs/matlab/mfile-plain
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
# R2010b seems not to be on the o.q machines.
if [ "$matlab_ver" = "r14" ]; then
    ## r14sp3 = version 7.1 from 2005
    matlab_binary=/usr/local/matlabr14sp3/bin/matlab
elif [ "$matlab_ver" = "r2010b" ]; then
    # R2010b = version 7.11 from 2010
    matlab_binary=/usr/local/MATLAB/R2010b/bin/matlab
else
    echo "$progname: only recognize matlab version r14 or r2010b"
    echo "$USAGE" 1>&2
    exit 2
fi
# needed for R2010b (see above); no harm otherwise
OSCHECK_ENFORCE_LIMITS=0; export OSCHECK_ENFORCE_LIMITS

# the basic matlab call
# notes:
#  the try/catch is good matlab to catch trouble
#    but, there is no way to return error status from matlab,
#    so we resort to looking for lines with ERROR
#    Thus, to play nice, the script has to start its error messages this
#    way.
#  fprintf(2,...) seems not to print to stderr.
#  sh builtin echo has trouble passing \n on some platforms
tmpnam="/tmp/$progname.$$.err"
rm -f "$tmpnam"
t0=`date +%s`
/bin/echo \
    "try," \
      "mask_series='$mask_series';" \
      "mag_series='$mag_series';" \
      "dest_dir='$dest_dir';" \
      "first_track=$first_track;" \
      "retain_history=$retain_history;" \
      "make_movie=$make_movie;" \
      "params='$params';" \
      "track_hmi_production;" \
    "catch," \
      "e=lasterror;" \
      "fprintf(1, '\\nERROR: %s in %s (Line %d)\\n', e.message, e.stack(1).name, e.stack(1).line);" \
    "end;" \
    "quit" | \
    $matlab_binary -nodesktop -nosplash -nodisplay -logfile "$tmpnam" 1>&2
t1=`date +%s`
tdiff=`expr $t1 - $t0`
echo "${progname}: Matlab call took $tdiff seconds."

# put the log file back
tmpnam2="$dest_dir/Tracks/jsoc/track-latest.log"
mkdir -p `dirname $tmpnam2` # dir present if cmd ran ok, maybe not if error
mv -f "$tmpnam" "$tmpnam2" 

# generate exit status
if grep -q "^ERROR" "$tmpnam2" ; then
    # grep above matched an error
    exit 1
else
    # no error found
    exit 0
fi

