#!/bin/sh
#
# track_revert_run.sh: Help revert to known tracker state after bad run.
#
# If a tracker run does not terminate cleanly, the file system may contain
# partial results.  After fixing the problem, one would wish to re-run the
# tracker starting at the last known-good state.  This script helps to 
# determine that state.
#
# Usage:
#
#  track_revert_run.sh dest_dir
#
# where:
#   dest_dir is a staging area in the filesystem, 
#     typically with a subdirectory .../Tracks/jsoc
#
#

# turmon nov 2011

# checkpoint series in SUMS
sums_ckpt=hmi.Mharp_log_720s

# usage: die LINENO MESSAGE CODE
# all arguments optional
function die() {
  echo "${progname}: Error near line ${1:- \?}: ${2:-(no message)}" 1>&2
  exit ${3:- 1}
}

progname=`basename $0`
USAGE="Usage: $progname [ -x ] dest_dir"

# get options
execute=0
execute_preface="  "
while getopts "x" o; do
    case "$o" in
	x)    execute=1
              execute_preface="Running: ";;
	[?h]) echo "$USAGE" 1>&2
	      exit 2;;
    esac
done
shift `expr $OPTIND - 1`

# get arguments
if [ "$#" -ne 1 ]; then
    echo "$USAGE" 1>&2
    exit 2
fi
# ok to set up args
dest_dir="$1"

#############################################################
#
# Announce intentions
#
if [ $execute -eq 1 ]; then
    echo ""
    echo "*** Executing shell code to revert tracks!"
    echo ""
fi

#############################################################
#
# Checkpoint file
#

root_dir=$dest_dir/Tracks/jsoc

if [ ! -d $root_dir ]; then
    die $LINENO "Could not find Tracks/jsoc within $dest_dir"
fi

# locate the checkpoint file
cat <<EOF 
Looking for a checkpoint file.  This is the most important step.
EOF
ckpt=$root_dir/track-prior.mat
if [ ! -r $ckpt ]; then
    echo "Could not find prior checkpoint within $root_dir"
    if [ ! -d $root_dir/State ]; then
	echo "No checkpoint directory $root_dir/State either."
	echo "You must check SUMS for the most recent checkpoint file."
	echo "Try $sums_ckpt"
    else 
	echo "Find the newest checkpoint in $root_dir/State"
	echo "and move it to $ckpt, or check SUMS for a checkpoint file."
	echo "Try $sums_ckpt"
	echo ""
	echo "Here are the checkpoints now in $root_dir/State:"
	ls -lt $root_dir/State/*.mat
    fi
    die $LINENO "Could not find checkpoint file within $root_dir"
fi
cat <<EOF
Found the default checkpoint file:

    $ckpt

This is good, you can re-run from it.  Here's its listing:

EOF
ls -l $ckpt
cat <<EOF

Or, check the directory

    $root_dir/State

for more checkpoints.  Here are a few of the recent checkpoints now there:

EOF
ls -lt $root_dir/State/*.mat | head
cat <<EOF

You may find the default checkpoint file

    $ckpt

duplicated with the same contents (but a different name) in this list.  
That is OK.  Also, you can check the checkpoint series in SUMS
for more checkpoint files.  Try $sums_ckpt.

If you can revert to the default checkpoint file listed above, 
you can restore the filesystem completely by following the 
next steps about removing orphaned tracks.  If you wish to use
an earlier checkpoint, you will have to remove all orphaned 
tracks from that checkpoint forward.  The following steps
do not help with that, i.e. recovery from non-default checkpoints.

Proceeding.
EOF

#############################################################
#
# Orphaned tracks (1)
#

cat <<EOF

========================================

EOF

# remove the finalized tracks from the partial run
cat <<EOF
Looking for finalized tracks from the partial run that have been
orphaned by its early termination.

EOF
badtk=$root_dir/track-new.txt
if [ ! -r $badtk ]; then
    echo "Could not find this list in $badtk"
    echo "You can try to resume without figuring out why this isn't there."
    die $LINENO "Could not find finalized-track file within $root_dir"
fi

cat <<EOF
Found the list:

    $badtk

We can use it to remove the orphaned finalized tracks.
Try executing the following remove commands (if any):

EOF
ok=1
empty=1
for tk in `grep -v '^#' $badtk` ; do
    empty=0
    tk0=`printf '%06d' $tk`
    tkfull=$root_dir/track-$tk0
    if [ ! -d $tkfull ]; then
	ok=0
	echo "No directory $tkfull for $tk, surprising."
    else
	echo "$execute_preface" rm -rf $tkfull
	if [ $execute -eq 1 ]; then
	  rm -rf $tkfull
	fi
    fi
done
if [ $empty -eq 1 ]; then
    echo "There were no orphaned tracks.  This is OK."
fi
if [ $ok -eq 0 ]; then
    echo "Found one or more missing track directories."
    echo "You probably want to take a closer look before removing anything."
    echo ""
    die $LINENO "Missing track directories within $root_dir"
fi
cat <<EOF

Existing orphaned tracks seem consistent with info from

    $badtk

It is probably ok to do the rm(s) above, if any.
Proceeding.

EOF

#############################################################
#
# Orphaned tracks (2)
#

cat <<EOF

========================================

EOF

# help remove the finalized tracks from the cumulative file
cat <<EOF
Looking for the cumulative finalized tracks from the partial run.
This list can also have orphaned entries.

EOF
alltk=$root_dir/track-final.txt
if [ ! -r $alltk ]; then
    echo "Could not find this list in $alltk"
    echo "You can try to resume without figuring out why this isn't there."
    die $LINENO "Could not find cumulative finalized-track file within $root_dir"
fi
cat <<EOF
Found it.  This is good.  You should edit the file 

    $alltk

and remove the last, unfinished block of tracks at EOF.
These entries, if any, are expected to duplicate those in

    $badtk

Note: this editing step is not necessary.

EOF

# Exit OK
echo "$progname: Done.  Looks OK."
exit 0

