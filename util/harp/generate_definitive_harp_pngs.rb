#! /usr/bin/ruby

# Name     generate_definitive_harp_pngs.rb
#
# Purpose  This script generates .png files showing all the definitive HARPs on the disk
#          using Mike Turmon's script, track_hmi_harp_movie_driver.sh.
#
# Machine  Use only on solar3 (the only machine that supports both ruby and MATLAB). 
#
#          This script runs in two modes, forward and backward. 
#          Forward mode -- this is the mode to run every day after the tracker run. 
#          It creates a definitive HARP .png file at T_REC = (now - 27 days). 
#          The forward mode only takes the arguments [TYPE] and [DIR] (see below).
#
#          Backward mode -- This is the mode to run to generate definitive HARP .pngs
#          for any period in HMI's history. The backward mode takes the arguments
#          [TYPE] [TSTART], [TSTOP], and [DIR].
#          Note: It takes 4 minutes to generate 24 images (a day's worth, at once per hour).
#                Run in backward mode such that it does not swamp solar3. 
# 
# Usage    > ruby generate_definitive_harp_pngs.rb [TYPE] [DIR] [TSTART] [TSTOP]  
#          where [TYPE]   equals 'forward' or 'backward'
#          where [DIR]    equals the path at which the images will be stored
#          where [TSTART] equals the start time, in TAI
#          where [TSTOP]  equals the stop time, in TAI
#
# E.g.     > ruby generate_definitive_harp_pngs.rb forward /web/jsoc/htdocs/doc/data/hmi/harp/harp_definitive/
#          > ruby generate_definitive_harp_pngs.rb backward /web/jsoc/htdocs/doc/data/hmi/harp/harp_definitive/ 2012.01.01_00:00:00_TAI 2012.01.01_06:00:00_TAI 
#
# Written by Monica Bobra 05 June 2013
require 'date'

type   = ARGV[0].to_s
dir    = ARGV[1].to_s
tstart = ARGV[2].to_s
tstop  = ARGV[3].to_s

if ARGV.empty?
 $stderr.puts "See header: usage is ruby generate_definitive_harp_pngs.rb [TYPE] [DIR] [TSTART] [TSTOP]"
 exit
end

# Run in backward mode.
if type == 'backward'
`/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts/track_hmi_harp_movie_driver.sh -f 'hmi.Marmask_720s[#{tstart}-#{tstop}@1h]' hmi.Mharp_720s #{dir}`
end


# Run in forward mode.
if type == 'forward'

# Determine the range of times for which to compute images
horizon = `show_info -q key=T_HORIZN ds=hmi.Mharp_log_720s"[]" n= -2`.split

if horizon[0] == horizon[1] 
 $stderr.puts "There have been no updates to hmi.Marmask_720s and hmi.Mharp_720s data since this script was last run."
 exit
end

# TIME A
# Determine the last T_REC at which hmi.Marmask_720s and hmi.Mharp_720s data are available 
# This is 12 minutes before the result of show_info 'ds=hmi.Mharp_log_720s[$]' key=T_HORIZN
t1a      = DateTime.strptime(horizon[1], '%Y.%m.%d_%H:%M:%S_TAI').strftime('%s').to_i - (60 * 12)
t2a      = Time.at(t1a).utc
t3a      = t2a.strftime('%Y')+'.'+t2a.strftime('%m')+'.'+t2a.strftime('%d')+'_'+t2a.strftime('%H')+':00:00_TAI'

# TIME B
# Determine the first T_REC at which hmi.Marmask_720s and hmi.Mharp_720s data are available 
# This is the first result of show_info -q key=T_REC ds=hmi.Mharp_log_720s"[]" n= -2.
# However, since this may not be on the hour (i.e., at _HH:00:00_TAI), go back one hour.
t1b      = DateTime.strptime(horizon[0], '%Y.%m.%d_%H:%M:%S_TAI').strftime('%s').to_i - (60 * 60)
t2b      = Time.at(t1b).utc
t3b      = t2b.strftime('%Y')+'.'+t2b.strftime('%m')+'.'+t2b.strftime('%d')+'_'+t2b.strftime('%H')+':00:00_TAI'

# Generate once-per hour images between TIME B and TIME A:
`/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts/track_hmi_harp_movie_driver.sh -f 'hmi.Marmask_720s[#{t3b}-#{t3a}@1h]' hmi.Mharp_720s #{dir}`

end
