#! /bin/csh -f

# Name     generate_definitive_image.csh
#
# Purpose  This script generates .png files showing all the definitive HARPs on the disk
#          using Mike Turmon's script, track_hmi_harp_movie_driver.sh.
#
# Machine  Must run on n02, n04, or solar3 (only machines where MATLAB runs). 
#
#          This script runs in two modes, forward and backward. 
#          Forward mode -- this is the mode to run every day after the tracker run. 
#          It creates a definitive HARP .png file every hour between T_HORIZ at T_REC = $
#          and T_HORIZ at T_REC = $ - 1. 
#          The forward mode only takes the arguments [TYPE] (see below).
#
#          Backward mode -- This is the mode to run to generate definitive HARP .pngs
#          for any period in HMI's history. The backward mode takes the arguments
#          [TYPE] [TSTART], and [TSTOP].
# 
# Usage    > ruby generate_definitive_harp_pngs.rb [TYPE] [TSTART] [TSTOP]  
#          where [TYPE]   equals 'forward' or 'backward'
#          where [TSTART] equals the start time, in TAI
#          where [TSTOP]  equals the stop time, in TAI
#
# E.g.     > ./generate_definitive_image.csh forward
#          > ./generate_definitive_image.csh backward 2012.01.01_00:00:00_TAI 2012.01.01_06:00:00_TAI 
#

set type   = $argv[1]

if ( $#argv < 1 ) then
  echo ""
  echo "See header: usage is ./generate_definitive_image.csh [TYPE] [TSTART] [TSTOP]"
  echo ""
  exit 1
endif

# Run in backward mode.
if ($type == 'backward') then
set tstart = $argv[2]
set tstop  = $argv[3]

/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts/track_hmi_harp_movie_driver.sh -f 'hmi.Marmask_720s['$tstart'-'$tstop'@1h]' hmi.Mharp_720s /web/jsoc/htdocs/doc/data/hmi/harp/harp_definitive/
endif

# Run in forward mode.
if ($type == 'forward') then

  # Determine the range of times for which to compute images
  set horizon = `show_info -q key=T_HORIZN ds=hmi.Mharp_log_720s"[]" n= -2`

  if ($horizon[1] == $horizon[2]) then
    echo ""
    echo "There have been no updates to hmi.Marmask_720s and hmi.Mharp_720s data since this script was last run."
    echo ""
    exit 1 
  endif

  # TIME A
  # Determine the last T_REC at which hmi.Marmask_720s and hmi.Mharp_720s data are available 
  # and subtract enough minutes such that it is on the hour.
  set tA  = `echo $horizon[1] | cut -c 0-14`'00:00_TAI'

  # TIME B
  # Determine the first T_REC at which hmi.Marmask_720s and hmi.Mharp_720s data are available 
  # This is the first result of show_info -q key=T_REC ds=hmi.Mharp_log_720s"[]" n= -2.
  # However, since this may not be on the hour (i.e., at _HH:00:00_TAI), go back one hour.
  set tB  = `echo $horizon[2] | cut -c 0-14`'00:00_TAI'

  # Generate once-per hour images between TIME B and TIME A:
  /home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts/track_hmi_harp_movie_driver.sh -f 'hmi.Marmask_720s['$tA'-'$tB'@1h]' hmi.Mharp_720s /web/jsoc/htdocs/doc/data/hmi/harp/harp_definitive/

endif