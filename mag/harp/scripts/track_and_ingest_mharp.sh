#!/bin/sh
# (line below is for qsub)
#$ -S /bin/sh
#
# Top-level driver for deriving HARPs from a series of HMI masks (e.g., hmi.Marmask_720s).
#
# Invokes HARP finder/tracker and ingests resulting HARPs into a harp data series.  
# Uses a filesystem directory (dest_dir) to transfer results.  Essential aspects of the 
# dest_dir results ("checkpoint files") are stored into harp_log_series in case dest_dir 
# is accidentally mangled.
#
# Exit status: 0 for no problems, 1 for an incorrect parameter (e.g., time range) that
# needs to be corrected, 2 for an error part-way through.  If the exit status is 2, you
# must be sure that the next run includes the same time range as the failed run.  You 
# must not skip the error-causing time range and continue with subsequent times.
#
# Usage:
#
#  track_and_ingest_mharp.sh [-nMaefd] [-i N] [-g N] dest_dir ...
#    mask_series harp_series harp_log_series
#
# where:
#   dest_dir is a staging area in the filesystem
#   mask_series is an input data series qualified by T_REC 
#   harp_series is an output data series to put tracks
#
# and:
#   -i indicatea an initial run, having new tracks numbered starting from 
#      the argument N.  Otherwise, the run picks up where an earlier one 
#      left off, according to a checkpoint file in `dest_dir'.
#      IF GIVEN, `-i N' CLEARS ALL EXISTING RESULTS.
#   -g says to use gap-filling mode for new tracks, and supplies a 
#      required argument N.  In this mode, new tracks will be checked
#      for overlap with future tracks; if so, the track ID will be assigned
#      to match.  If there is no overlap, new tracks will get an ID that 
#      was left unused due to earlier merges.  This means new tracks 
#      during the run will (mostly) not affect numbering of later tracks.  
#      The argument, typically 0, is a number to add to the track ID
#      counter (ROI_rgn.Ntot) at the end of the run, in case more help 
#      is needed to synch the gap-fill count with the desired count.
#      -g is not currently used with -n.
#   -n indicates near-real-time (NRT) mode: (1) NRT magnetograms in
#      hmi.M_720s_nrt are used (otherwise, hmi.M_720s is used).  
#      (2) Less track history is retained, allowing faster startup and 
#      more compact checkpoint files.  (3)  Both
#      finalized and currently-pending tracks are ingested into DRMS.
#   -M indicates MDI mode: (1) MDI magnetograms in mdi.fd_M_96m_lev182
#      are used; (2) force-run (-f) and developer paths (-d) implied;
#      (3) Some tracker parameters change.
#   -a says that currently-pending tracks should be ingested also.
#      (-a = ingest all tracks).  -n implies -a, but sometimes -a
#      is useful without -n for testing definitive tracks, to ingest 
#      even incomplete tracks into DRMS for inspection.
#   -e means to skip running the tracker, and only do the ingestion of 
#      tracks in the filesystem (dest_dir) into DRMS.  This can be used
#      if the ingestion phase was interrupted after tracking was complete,
#      and you don't want to run the tracker again.  It can also be used to
#      test ingestion while holding the filesystem constant.
#      Mnemonic: -e = ingest (eat) only.  Both tracks and logs are ingested.
#   -f means to force the run, without checking before-hand for gaps in the
#      relevant Marmask series.  For experts.
#   -d is a flag which signals to use the developer MATLABPATH (~turmon) rather
#      than the regular production path (~jsoc).  For experts.
#
# Note that we require a T_REC to be supplied with mask_series.  
# 
# Typical command lines:
# [definitive]
#   track_and_ingest_mharp.sh -i 1 /tmp/harp/definitive hmi.Marmask_720s[2011.10.01/10d]
#     hmi.Mharp_720s hmi.Mharp_log_720s
# [nrt]
#   track_and_ingest_mharp.sh -i 1 -n /tmp/harp/nrt hmi.Marmask_720s_nrt[2011.10.01_12:36_TAI]
#     hmi.Mharp_720s_nrt hmi.Mharp_log_720s_nrt
# Don't forget to omit the -i 1 argument on subsequent runs, in both cases!
#
# This routine calls: track_hmi_production_driver_stable.sh, ingest_mharp, ingest_mharp_log
#

# turmon jul-sep 2011, oct 2012, feb 2013, july 2013

# echo command for debugging
#set -x
# exit-on-error: cannot put this in #! line because qsub does not read that line
set -e

progname=`basename $0`
# under SGE, $0 is set to a nonsense name, so use our best guess instead
default_progname=track_and_ingest_mharp.sh
if [ `expr "$progname" : '.*track.*'` == 0 ]; then progname=$default_progname; fi
USAGE="Usage: $default_progname [-nMaefd] [-i N] [-g N] dest_dir mask_series harp_series harp_log_series"

## # tempfile stuff: commented out
## TEMPROOT=/tmp/harps-${progname}
## # creates the dir, and returns the dir name
## TEMP_DIR=`mktemp -d $TEMPROOT-XXXXXXXX`
## # filenames, e.g. below, must match the deletion pattern in cleanup()
## TEMPMASK="$TEMP_DIR/mask.list"

# we have not begun tracking yet, so errors are not problematic
#   (after we begin, will set exit status to 2)
HAVE_BEGUN=0 && EXIT_STATUS=1

# echo the arguments, for repeatability
echo "${progname}: Invoked as:" "$0" "$@"
# elapsed time
time_started=`date +'%s'`

# This routine can be entered via implicit exit from script due to a failed command,
# or via die(), which is invoked explicitly at a few places below.
function cleanup() {
    # clean the whole temp dir (commented out)
    # rm -rf "$TEMP_DIR"
    exit $EXIT_STATUS
}
trap cleanup EXIT

# Trap error exits through this routine
#   especially important for shell-command exits (trap ... ERR below)
#
# usage: die LINENO MESSAGE
#   arguments optional
function die() {
    echo "${progname}: Fatal error near line ${1:- \?}: ${2:-(no message)}" 1>&2
    if [ $HAVE_BEGUN = 0 ]; then
	echo "${progname}: Exiting with simple error."
	echo "${progname}: Correct input parameters and re-run."
	echo "${progname}: Correct input parameters and re-run." 1>&2
    else
	echo "${progname}: Exiting with error: Filesystem state is OK (but inconsistent)."
	echo "${progname}: Can run again, but must use same (or longer) T_REC range."
	echo "${progname}: Do *not* proceed to a new T_REC range."
	echo "${progname}: Can run again, but must use same (or longer) T_REC range." 1>&2
	echo "${progname}: Do *not* proceed to a new T_REC range." 1>&2
    fi
    # will trap through cleanup()
    exit
}

# Error message selection upon exit
cmd="shell code"
trap 'die $LINENO "Running $cmd"'     ERR
trap 'die $LINENO "Received SIGHUP"'  SIGHUP
trap 'die $LINENO "Received SIGINT"'  SIGINT
trap 'die $LINENO "Received SIGTERM"' SIGTERM

# log_and_run: write command to stdout stream, and then execute it
#   writes a full command line plus parameters to stdout, and then 
#   executes it.  Doing this allows repeatability of the subsidiary 
#   command.  
#   The log output needs to go to stdout, not stderr, because 
#   outputs to stderr are reserved for error conditions.  Thus, set -x
#   will not work (it always goes to stderr).  We could get around that by
#   forcing the set -x and the subsidiary command into a subshell, and sending
#   stderr from the subshell to stdout, but this would redirect stderr from the
#   subsidiary command as well, which is unacceptable.  
#   Also note: error exit from subshells does not result in exit from 
#   this script even though -e is set in this script.
#   Note that the construct here will not work if the command itself contains
#   redirection.
function log_and_run() {
  echo "${progname}: Executing command:" "$@"
  # "trap ERR" is not inherited in functions, so explicitly call die on failure
  "$@" || die $LINENO "Running $cmd"
  echo "${progname}: Returned OK from command $cmd."
}

# get options -- pass them down to tracker/ingestor
make_movie=1
first_track=0
developer_path=0
track_opts=""
ingest_all=0
nrt_mode=0
mdi_mode=0
skip_roi=0
gap_fill=0
force_run=0
eat_only=0
while getopts "hnMadfeg:i:" o; do
    case "$o" in
	i)  first_track="$OPTARG"
	    track_opts="$track_opts -i $first_track";;
	g)  skip_roi="$OPTARG"
	    gap_fill=1;;
	d)  developer_path=1
	    track_opts="$track_opts -d";;
	a)  ingest_all=1;;
	e)  eat_only=1;;
	f)  force_run=1;;
	n)  mdi_mode=0
	    nrt_mode=1;;
	M)  nrt_mode=0
	    mdi_mode=1;;
     [?h])  echo "$USAGE" 1>&2
	    exit 2;;
    esac
done
shift `expr $OPTIND - 1`

# Disable globbing from here forward, to prevent wildcard chars, esp. within
# $mask_series, from expanding into files
set -f

# handle setup for NRT versus definitive mode
if [ "$nrt_mode" -eq 1 ]; then
    # tracker options:
    # nrt mags + pgrams
    mag_series=hmi.M_720s_nrt
    pgram_series=hmi.Ic_noLimbDark_720s_nrt
    #   retain only r days of history
    track_opts="$track_opts -r 1"
    # ingest both new and pending regions -- always do pending in this case
    listfiles="track-new.txt,track-pending.txt"
    # ingester options:
    #   ingest only most recent trec; no temporal padding
    #   do not filter HARPs based on number of appearances (tmin=1)
    #   do not output match information
    ingest_opts="trec=1 tpad=0 tmin=1 match=0"
    do_match=0
elif [ "$mdi_mode" -eq 1 ]; then
    # MDI tracker options:
    #   definitive mags + pgrams
    mag_series=mdi.fd_M_96m_lev182
    pgram_series=mdi.fd_M_96m_lev182 # unused in presence of force_run=1
    #   (set up mdi-specific tracker params this way)
    # track_opts="$track_opts -p XXX"
    #   mdi-specific tracker parameters
    track_opts="$track_opts -s mdi_track_setup"
    #   don't check for gaps
    force_run=1
    #   developer paths
    developer_path=1
    track_opts="$track_opts -d"
    #   retain all history
    track_opts="$track_opts -r inf"
    if [ "$ingest_all" -eq 0 ]; then
        # ingest only new regions -- usual path
        listfiles="track-new.txt"
    else
        # ingest both new and pending regions -- useful for tests
        listfiles="track-new.txt,track-pending.txt"
    fi
    # ingester options:
    #   accept default trec
    #   set tpad to 1 day (15 images)
    #   set cadence to 96 minutes (5760s)
    #   filter M-TARPs having fewer than tmin appearances
    #   do output match information
    ingest_opts="tmin=3 tpad=15 cadence=5760 match=1"
    do_match=1
else
    # regular HMI tracker options:
    # definitive mags + pgrams
    mag_series=hmi.M_720s
    pgram_series=hmi.Ic_noLimbDark_720s
    #   retain all history
    track_opts="$track_opts -r inf"
    if [ "$ingest_all" -eq 0 ]; then
        # ingest only new regions -- usual path
        listfiles="track-new.txt"
    else
        # ingest both new and pending regions -- useful for tests
        listfiles="track-new.txt,track-pending.txt"
    fi
    # ingester options:
    #   accept default trec, tpad
    #   filter HARPs having fewer than tmin appearances
    #   do output match information
    ingest_opts="tmin=3 match=1"
    do_match=1
fi

# movie vs. not
if [ "$make_movie" -eq 1 ]; then
    # -m makes a movie, its absence does not
    track_opts="$track_opts -m"
fi

# gap-fill vs. not
if [ "$gap_fill" -eq 1 ]; then
    # -g N implies gap-filling, its absence implies not
    #   (note that -g 0 is not the same as no -g at all)
    track_opts="$track_opts -g $skip_roi"
fi

# get arguments
if [ "$#" -ne 4 ]; then
    echo "$USAGE" 1>&2
    exit 2
fi
# ok to set up args
dest_dir="$1"
mask_series="$2"
harp_series="$3"
harp_log_series="$4"

# input data series, without T_REC qualifier
mask_series_only=`echo $mask_series | sed 's/\[.*$//'`

# SGE/OpenMP setup
SGE_ROOT=/SGE;     export SGE_ROOT
OMP_NUM_THREADS=1; export OMP_NUM_THREADS
KMP_BLOCKTIME=10;  export KMP_BLOCKTIME

uname=`uname -m`
if [ "X$uname" = Xi686 ]; then
    PATH="${PATH}:$SGE_ROOT/bin/lx24-x86"
elif [ "X$uname" = Xx86_64 ]; then
    PATH="${PATH}:$SGE_ROOT/bin/lx24-amd64"
elif [ "X$uname" = Xx86_64 ]; then
    PATH="${PATH}:$SGE_ROOT/bin/lx24-ia64"
else
    die $LINENO "Could not find system '$uname' to set up SGE"
fi

#############################################################
#
# Simple error checks, while it's painless
#
#############################################################

# if not -i, check that dest_dir exists and is writable, etc.
if [ $first_track = 0 ]; then
    if [ ! -d "$dest_dir/Tracks/jsoc" ]; then
	die $LINENO "Not initial run, but could not find Tracks/jsoc within dest_dir"
    fi
    if [ ! -w "$dest_dir/Tracks/jsoc" -o ! -x "$dest_dir/Tracks/jsoc" ]; then
	die $LINENO "Not initial run: need write+execute permission on Tracks/jsoc within dest_dir"
    fi
    if [ ! -w "$dest_dir/Tracks/jsoc/track-prior.mat" ]; then
	die $LINENO "Not initial run: could not find/write Tracks/jsoc/track-prior.mat within dest_dir"
    fi
    if [ -e "$dest_dir/Tracks/jsoc/track-post.mat" -a ! -w "$dest_dir/Tracks/jsoc/track-post.mat" ]; then
	die $LINENO "Need write permission on existing file Tracks/jsoc/track-post.mat within dest_dir"
    fi
fi

# ensure track_hmi_production_driver_stable.sh is in PATH
if [ "$developer_path" -eq 1 ]; then
    # put script dir of developer into path
    # (could equally be .../proj/mag/harp/scripts)
    PATH="${PATH}:$HOME/cvs/JSOC/proj/myproj/scripts"
else
    # put script dir of JSOC into path
    PATH="${PATH}:/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts"
fi

# for show_info, ingestor, and log-ingestor
# there is no sh version of .setJSOCenv, so we use this to get JSOC_MACHINE
jsoc_mach=`csh -f -c 'source /home/jsoc/.setJSOCenv && echo $JSOC_MACHINE'`
PATH="${PATH}:$HOME/cvs/JSOC/bin/$jsoc_mach"

echo "${progname}:   show_info is at:" `which show_info`
echo "${progname}:   ingest_mharp is at:" `which ingest_mharp`
echo "${progname}:   tracker driver is at:" `which track_hmi_production_driver_stable.sh`

# check that relevant data series exist
cmd="show_info -s"
$cmd "$mag_series"      > /dev/null
$cmd "$mask_series"     > /dev/null
$cmd "$harp_series"     > /dev/null
$cmd "$harp_log_series" > /dev/null

#############################################################
#
# Main processing starts
#
#############################################################

# stage data (include T_REC)
#   (this does not fail if the series is empty)
echo "${progname}: Staging masks."
cmd="show_info -q -ip"
$cmd "$mask_series" > /dev/null
echo "${progname}: Finished staging masks."
cmd="shell code"

# recover a "first" and "last" T_REC, used to ingest the log.
# use the mag series, not the mask series:
#   reason 1: mask series can be empty while mag is not (bad mag => no mask)
#   reason 2: empty mask series means no harps, but still need to ingest the log
echo "${progname}: Looking up T_REC boundaries in mags."
mag_range=`echo "$mask_series" | sed "s/[^[]*/${mag_series}/"`
# find first T_REC, for later (and last for info)
cmd="show_info -q key=T_REC"
trec_frst=`$cmd "$mag_range" | head -n 1`
trec_last=`$cmd "$mag_range" | tail -n 1`
cmd="shell code"
if [ -z "$trec_frst" ]; then
  die $LINENO "Could not find a valid first T_REC within $mask_series (mags: ${mag_range})"
fi
echo "${progname}: Processing from $trec_frst to $trec_last"
echo "${progname}: Associated T_REC in \`$harp_log_series' will be $trec_frst"

# find the most recent prior ingest using harp_log 
## set -x
DAYSBACK=60 # days back in the harp_log's to look for the prior ingest
CADENCE=120 # 720s images/day
trec_frst_index=`index_convert "ds=$harp_log_series" "T_REC=$trec_frst"`
trec_prior0_index=$(( $trec_frst_index - $DAYSBACK * $CADENCE ))
trec_prior1_index=$(( $trec_frst_index - 1 ))
trec_prior0=`index_convert "ds=$harp_log_series" "T_REC_index=$trec_prior0_index"`
trec_prior1=`index_convert "ds=$harp_log_series" "T_REC_index=$trec_prior1_index"`
trec_prior=`show_info -q "ds=${harp_log_series}[${trec_prior0}-${trec_prior1}]" key=T_REC | tail -n 1`
# if this is an initial ingest, blank the prior ingest time -- check it later
if [ $first_track -ne 0 ]; then
  trec_prior=
fi
# if above show_info was empty, it outputs "No records in query ..." -- convert to empty
if expr "$trec_prior" : ".*No records" > /dev/null ; then 
  trec_prior=
fi
## set +x
# announce the results
if [ -z "$trec_prior" ]; then
  echo "${progname}: Found no prior ingest in \`$harp_log_series' (before ${trec_frst})"
  prior_dt=
else
  echo "${progname}: Most recent prior ingest in \`$harp_log_series' was ${trec_prior}"
  # T_REC of the last frame in that ingest ("last mask we saw")
  trec_last_prior=`show_info -q "ds=${harp_log_series}[${trec_prior}]" key=T_LAST`
  # trec_last_prior_index=`index_convert "ds=$mag_series" "T_REC=$trec_last_prior"`
  # count how many *good* mags between trec_last_prior and trec_frst-1 (inclusive)
  #   1: expected (the 1 frame = trec_last_prior: it had a mask, so it had a good mag)
  #  >1: not OK (intervening valid mags)
  #   0: not good, but OK (probably overlapping intervals)
  prior_dt=`show_info -c "ds=${mag_series}[${trec_last_prior}-${trec_prior1}][? QUALITY>=0 ?]" key=T_REC | awk '{print $1}'`
fi

# ensure there are no un-explained mask gaps
#   (skip on nrt-mode because it will be too late to fix anyway)
if [ $nrt_mode = 0 -a $force_run = 0 ]; then
  echo "${progname}: Ensuring no gaps in masks."
  # 1: ensure no gap between end of prior run and start of current run
  #    (except if initial run, or gap-fill)
  if [ $first_track = 0 -a $gap_fill = 0 ]; then
    echo "${progname}:   Checking for gaps between ingests."
    if [ -z "$prior_dt" ]; then
      die $LINENO "Unable to find a prior ingest in \`$harp_log_series'."
    fi
    echo "${progname}:     Prior ingest final frame at ${trec_last_prior}"
    echo "${progname}:     Current start at ${trec_frst}"
    if [ "$prior_dt" -gt 1 ]; then
      die $LINENO "Prior ingest final frame was $prior_dt frames before current start."
    elif [ "$prior_dt" -lt 1 ]; then
      echo "${progname}:   Prior ingest end seems ahead of current start."
      echo "${progname}:   Perhaps overlapping T_REC ranges, which is OK.  Continuing."
    else
      echo "${progname}:   Detect no gap between ingests.  Continuing."
    fi
  fi
  # 2: ensure no gaps within masks, due to un-made masks
  echo "${progname}:   Checking for gaps within the given T_REC range."
  cmd="track_hmi_check_masks"
  if ! $cmd xm="$mag_series" xp="$pgram_series" y="$mask_series" ; then
    die $LINENO "Found a gap in mask series $mask_series.  Fill and re-run."
  fi
  cmd="shell code"
  echo "${progname}: Finished ensuring no gaps in masks."
else
  echo "${progname}: Skipping check for gaps in masks."
fi


# after this point, the filesystem may get out of sync because tracker alters $dest_dir
HAVE_BEGUN=1 && EXIT_STATUS=2

# driver invokes matlab to perform tracking
#  -- if not first run, reads checkpoint from track-prior.mat and writes to track-post.mat
#  -- note, track_opts is unquoted
if [ $eat_only = 0 ]; then
  echo "${progname}: Beginning tracking."
  cmd=track_hmi_production_driver_stable.sh 
  log_and_run "$cmd" $track_opts "$mask_series" "$mag_series" "$dest_dir"
  # "$cmd" $track_opts "$mask_series" "$mag_series" "$dest_dir" 
  cmd="shell code"
  echo "${progname}: Finished tracking."
else
  echo "${progname}: Skipping tracking, ingest-only selected."
fi

echo "${progname}: Ingesting tracks."

# ingest data
#  -- note, ingest_opts is unquoted
cmd=ingest_mharp 
log_and_run "$cmd" $ingest_opts \
    "mag=$mag_series" "root=$dest_dir/Tracks/jsoc" "lists=$listfiles" \
    "out=$harp_series" "mask=$mask_series_only" verb=1 
cmd="shell code"
# preserve NOAA AR match data
if [ $do_match = 1 ]; then
    # name for saved match info
    ext="$trec_frst"
    # current and new location for NOAA AR match data
    matchdir="$dest_dir/Tracks/jsoc/Match"
    matchold="$dest_dir/Tracks/jsoc/Match-old"
    if [ -d "$matchdir" ]; then
        # -p also suppresses errors if the dir exists already
	mkdir -p "$matchold"
        # ensure nothing is there from an earlier run
        rm -rf "$matchold/Match-$ext"
        mv "$matchdir" "$matchold/Match-$ext"
    else
        echo "${progname}: NOAA-match diagnostics \`$matchdir' unsaved because not present."
        echo "${progname}: (Earlier run was interrupted?)  DRMS unaffected.  Continuing."
    fi
fi

echo "${progname}: Finished ingesting tracks."
echo "${progname}: Ingesting log and checkpoint to ${harp_log_series}[${trec_frst}]"

# combine regular logs
logCbase="track-and-ingest.log"
logCfile="$dest_dir/Tracks/jsoc/$logCbase"
log1file="$dest_dir/Tracks/jsoc/track-latest.log"
log2file="$dest_dir/Tracks/jsoc/ingest-latest.log"
rm -f "$logCfile"
touch "$logCfile"
[ -r "$log1file" ] && cat "$log1file"           >> $logCfile
echo "========================================" >> $logCfile
[ -r "$log2file" ] && cat "$log2file"           >> $logCfile

# combine errror logs -- finicky, want it to not exist if there were no errors
errCbase="track-and-ingest.err"
errCfile="$dest_dir/Tracks/jsoc/$errCbase"
err1file="$dest_dir/Tracks/jsoc/track-latest.err"
err2file="$dest_dir/Tracks/jsoc/ingest-latest.err"
rm -f "$errCfile"
if [ -s "$err1file" -o -s "$err2file" ]; then
    touch "$errCfile"
    [ -s "$err1file" ] && cat "$err1file"           >> $errCfile
    echo "========================================" >> $errCfile
    [ -s "$err2file" ] && cat "$err2file"           >> $errCfile
fi

# ingest log and checkpoint file 
#   -- use modified log/err files created above
cmd=ingest_mharp_log
log_and_run "$cmd" "root=$dest_dir/Tracks/jsoc" "log=$logCbase" "err=$errCbase" \
    "out=$harp_log_series" "trec=$trec_frst" 
cmd="shell code"

echo "${progname}: Finished ingesting log and checkpoint."

# put newly-generated post-run checkpoint file into next-run checkpoint location
# (formerly cp -p, but it was complaining, which went to stderr and gave impression 
# of actual trouble)
cp -f "$dest_dir/Tracks/jsoc/track-post.mat" "$dest_dir/Tracks/jsoc/track-prior.mat"

time_ended=`date +'%s'`
time_taken=$(( $time_ended - $time_started ))
# Exit OK 
echo "${progname}: Clock time elapsed = $time_taken seconds"
echo "${progname}: Status OK."
echo "${progname}: Done."
# establish clean exit
EXIT_STATUS=0
# traps thru cleanup()
exit

