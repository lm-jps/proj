#!/bin/sh -e
# (line below is for qsub)
#$ -S /bin/sh
# 
# sh driver to run the Matlab HARP-frame finder on data from JSOC
#
# usage:
#   track_hmi_harp_movie_driver.sh [-f] [-m] [-d] mask_series harp_series dest_dir
# where:
#   mask_series is a mask data series from JSOC, with a T_REC qualifier
#   harp_series is a HARP series where the times above are looked up
#   dest_dir is a destination directory for file and state input and output
# and:
#   -f is a flag to generate individual frames from the results (default is no frames)
#   -m is a flag to generate a movie from the results (default is no movie)
#   -d is a flag which signals to use the developer MATLABPATH (~turmon) rather
#      than the regular production path.
# One or both of -f or -m should be supplied, otherwise this is a no-op.
#
# * Both masks and corresponding HARPs are needed.  Times are defined using the
# mask series, and the HARPs at those times are looked up in the HARP series.  
# Thus, the HARP series is NOT qualified by HARP number or time.
# * Thus, the mask_series is typically of the form: 
#     any_series[any_range] 
# for example,
#     hmi.Marmask_720s[2011.01.01_TAI/7d]
# The mask series sets the time range.
# * The harp series used is typically either hmi.Mharp_720s or hmi.Mharp_720s_nrt, 
# depending on whether NRT masks were specified or not.  For testing, other series
# may be used.
# * If dest_dir does not exist, it will be created.
# * Sample usage:
#
#  track_hmi_harp_movie_driver.sh -f -m \
#     'hmi.Marmask_720s[2010.09.20_12:00_TAI/1h]' hmi.Mharp_720s /tmp/sept20
#
# Michael Turmon, JPL, October 2012

progname=`basename $0`
USAGE="Usage: $progname [ -f ] [ -m ] [ -d ] mask_series harp_series dest_dir"

# get options
make_movie=0
make_frame=0
developer_path=0
matlab_ver=r2010b
while getopts "dmfv:" o; do
    case "$o" in
	f)    make_frame=1;;
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

echo "${progname}: Beginning."
# ok to get args -- note, we add {mask} on to the mask series
mask_series="$1{mask}"
harp_series="$2"
dest_dir="$3"

# set up mode -- beginning , is OK
mode=""
if [ "$make_frame" -eq 1 ]; then
    mode="$mode,frame"
fi
if [ "$make_movie" -eq 1 ]; then
    mode="$mode,movie"
fi
if [ -z "$mode" ]; then
    echo "${progname}: Need one of -m or -f" 1>&2
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
tmpnam="/tmp/$progname.$$.err"
rm -f "$tmpnam"
t0=`date +%s`
$matlab_binary -nodesktop -nosplash -nodisplay -logfile "$tmpnam" <<EOF
% use outer try/catch to catch *all* exceptions
try,
  mask_series='$mask_series';
  harp_series='$harp_series';
  dest_dir='$dest_dir';
  mode='$mode';
  track_hmi_harp_movie_script;
catch ME,
  % ensure these errors start on a new line
  fprintf(1, '\\nERROR: %s in %s (Line %d)\\n', ME.message, ME.stack(1).name, ME.stack(1).line);
  fprintf(2, '\\nERROR: %s in %s (Line %d)\\n', ME.message, ME.stack(1).name, ME.stack(1).line);
  rethrow(ME);
end;
quit;
EOF
t1=`date +%s`
tdiff=`expr $t1 - $t0`
echo "${progname}: Matlab call took $tdiff seconds."

# put the log file back
tmpnam2="$dest_dir/frame-latest.log"
mkdir -p `dirname $tmpnam2` # dir present if cmd ran ok, maybe not if error
mv -f "$tmpnam" "$tmpnam2" 

# generate exit status
if grep -q "^ERROR" "$tmpnam2" ; then
    # grep above matched an error
    echo "${progname}: Exiting with error." 
    echo "${progname}: Exiting with error." 1>&2
    exit 1
else
    # no error found
    echo "${progname}: Done."
    exit 0
fi

