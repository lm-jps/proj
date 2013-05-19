#!/bin/sh
# (line below is for qsub)
#$ -S /bin/sh
# 
# Top-level sh driver for movie-making using tracker diagnostic files.
# Useful if you forgot to make movies in a normal tracker run.
#
# usage:
#   track_hmi_production_driver_movie.sh [-p] tag dest_dir
# where:
#   tag is a text tag used to name the movie, typically containing a T_REC range
#   dest_dir is the destination directory for state input, and movie output
# and:
#   -d is a flag which signals to use the developer MATLABPATH (~turmon) rather
#      than the regular production path.
#
# * The dest_dir is the same one given to track_hmi_production_driver.
# It should contain subdirectories called "Tracks" and "Tracks/jsoc".
# * tags should contain no filesystem metacharacters like [, (, {, etc.
#
# For example,
#
#  track_hmi_production_driver_movie.sh hmi.Marmask_720s_2012_sept /tmp/sept-run
#
# See also: track_hmi_production_driver.sh
#
# Michael Turmon, JPL, December 2010, June 2011


progname=`basename $0`
USAGE="Usage: $progname [ -d ] tag dest_dir"

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
if [ "$#" -ne 2 ]; then
    echo "$USAGE" 1>&2
    exit 2
fi

echo "${progname}: Beginning."
# ok to get them -- note, we add {mask} on to the mask series
tag="$1"
dest_dir="$2"

# set up matlab path
if [ "$developer_path" -eq 1 ]; then
  # the tree under ~turmon/matlab/mfile -- development versions
  mtHome=`echo ~turmon`
  MATLABPATH=`find $mtHome/matlab/mfile -maxdepth 1 -type d | paste -d : -s -`
else
  # the tree under CVS -- production versions
  rootDir=$HOME/cvs/JSOC/proj/myproj/libs
  MATLABPATH="$rootDir/matlab/mfile-mex/fits:$rootDir/matlab/mfile-mex/assignment:$rootDir/matlab/mfile-mex/hmi-mask-patch:$rootDir/matlab/mfile-mex/standalone:$rootDir/matlab/mfile-plain"
fi
 
export MATLABPATH

#
# Matlab version
#
# The old JSOC matlab is r14sp3, which is certified to work with glibc 2.3.4, 
# which is the kernel version on n02.  The movie-maker works on n02 with the 
# old matlab, but not on n06 etc. (nov 2012)  The NOAA AR location determination
# queries JSOC, requiring the Java-based URL-getter, which throws an error, perhaps
# due to a JAVA path issue (LD_LIBRARY_PATH) on n06 et al.
# The new JSOC matlab is R2010b, which wants glibc 2.7, but appears to work
# OK with our older kernel.  R2010b does a check at startup which
# issues a prompt to continue or abort.  The prompt can be disabled with
# the OSCHECK_ENFORCE_LIMITS line.  It still gives a stern warning.
## r14sp3 = version 7.1 from 2005
#matlab_binary=/usr/local/matlabr14sp3/bin/matlab
# R2010b = version 7.11 from 2010
matlab_binary=/usr/local/MATLAB/R2010b/bin/matlab
# needed for R2010b (see above)
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
  tag='$tag';
  dest_dir='$dest_dir';
  track_hmi_production_movie_only;
catch ME,
  % ensure these errors start on a new line
  fprintf(1, '\\nERROR: %s in %s (Line %d)\\n', ME.message, ME.stack(1).name, ME.stack(1).line);
  fprintf(2, '\\nERROR: %s in %s (Line %d)\\n', ME.message, ME.stack(1).name, ME.stack(1).line);
  rethrow(ME);
end;
quit
EOF
t1=`date +%s`
tdiff=$(( $t1 - $t0 ))
echo "${progname}: Matlab call took $tdiff seconds."

# generate exit status
if grep -q "^ERROR" "$tmpnam" ; then
    # grep above matched an error
    echo "${progname}: Exiting with error." 
    echo "${progname}: Exiting with error." 1>&2
    exit 1
else
    # no error found
    echo "${progname}: Done."
    exit 0
fi

