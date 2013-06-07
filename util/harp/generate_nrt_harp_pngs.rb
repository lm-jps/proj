#! /usr/bin/ruby

# Name     generate_nrt_harp_pngs.rb
#
# Purpose  This script generates .png files showing all the nrt HARPs on the disk
#          using Mike Turmon's script, track_hmi_harp_movie_driver.sh.
#
#          This script generates three files that are overridden per run of this module:
#          latest.png -- a negative-color latest nrt image used for nrt sharp data visualization
#          thumbnail.png -- a thumbnail latest nrt image used for the webpage
#          latest_nrt.png  -- a latest nrt image
#
#          This script also keeps two months of harp.png files and throws away all
#          harp nrt .png files older than two months.
# 
#          This script should be run in cron every 12 minutes using the path below.
#
# Machine  Use only on solar3 (the only machine that supports both ruby and MATLAB). 
#
# Usage    > ruby generate_definitive_harp_pngs.rb [DIR]  
#          where [DIR]    equals the path at which the images will be stored
#
# E.g.     > ruby generate_nrt_harp_pngs.rb /web/jsoc/htdocs/doc/data/hmi/harp/harp_nrt/
#
# Written by Monica Bobra 05 June 2013
require 'date'

dir   = ARGV[0].to_s
Dir.chdir ARGV[0]

if ARGV.empty?
 $stderr.puts "See header: usage is ruby generate_nrt_harp_pngs.rb [DIR]"
 exit
end

trec = `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=T_REC ds='hmi.Mharp_720s_nrt[][$]' | uniq`.strip

`/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts/track_hmi_harp_movie_driver.sh -f 'hmi.Marmask_720s_nrt[#{trec}]' hmi.Mharp_720s_nrt #{dir}`

filename = "#{dir}harp.#{trec}.png"

# Create latest nrt image
`cp #{filename} latest_nrt.png`
`convert latest_nrt.png -fill white -gravity North -pointsize 36 -font Helvetica -annotate 0 'near real-time (nrt) data' latest_nrt.png`

# Create negative color image for nrt data visualization
`cp #{filename} latest.png`
`convert latest.png -negate latest.png`

# Create thumbnail
`convert -define png:size=1024x1024 #{filename} -thumbnail 256x256 -unsharp 0x.5 thumbnail.png`

#t     = DateTime.strptime(trec, '%Y.%m.%d_%H:%M:%S_TAI').to_time - (60*60*24*60) # would work in ruby 1.9.2
t     = DateTime.strptime(trec, '%Y.%m.%d_%H:%M:%S_TAI').strftime('%s').to_i - (60*60*24*60) # BUT we are dealing with ruby 1.8.7

files = Dir['harp.*.png']

files.each do |file|
  if DateTime.strptime(file[5,19],'%Y.%m.%d_%H:%M:%S').strftime('%s').to_i < t
    `rm -f #{file}`
  end
end
