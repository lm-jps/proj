#! /usr/bin/ruby

# Name     generate_definitive_harp_pngs.rb
#
# Purpose  This script generates .png files showing all the definitive HARPs on the disk
#          using Mike Turmon's script, track_hmi_harp_movie_driver.sh.
#
# Machine  Use only on solar3 (the only machine that supports both ruby and MATLAB). 
#
#          This script runs in two modes, forward and backward. 
#          Forward mode -- this is the mode to run in cron every hour. 
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

t1  = Time.now.utc
t2  = t1 - (60 * 60 * 24 * 27) #27 days in seconds
t3  = t2.strftime('%Y')+'.'+t2.strftime('%m')+'.'+t2.strftime('%d')+'_'+t2.strftime('%H')+':00:00_TAI'

`/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts/track_hmi_harp_movie_driver.sh -f 'hmi.Marmask_720s[#{t3}]' hmi.Mharp_720s #{dir}`

end