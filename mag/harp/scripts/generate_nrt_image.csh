#! /bin/csh -f

# Name     generate_nrt_image.csh
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
#          This script should be run in cron every hour using the path below.
#
# Machine  Must run on n02, n04, or solar3 (only machines where MATLAB runs).
#
# Usage    > ./generate_nrt_image.csh
#
# E.g.     > ./generate_nrt_image.csh

# Set directory
cd /web/jsoc/htdocs/doc/data/hmi/harp/harp_nrt/

# Find latest T_REC
set TREC = `/home/jsoc/cvs/JSOC/bin/linux_x86_64/show_info -q key=T_REC ds='hmi.Mharp_720s_nrt[][$]' | uniq`

# Run tracker
/home/jsoc/cvs/Development/JSOC/proj/mag/harp/scripts/track_hmi_harp_movie_driver.sh -f 'hmi.Marmask_720s_nrt['$TREC']' hmi.Mharp_720s_nrt /web/jsoc/htdocs/doc/data/hmi/harp/harp_nrt/

# Identify latest png file
set filename = /web/jsoc/htdocs/doc/data/hmi/harp/harp_nrt/harp.$TREC.png

# Create latest nrt image
cp $filename latest_nrt.png
convert latest_nrt.png -fill white -gravity North -pointsize 36 -font Helvetica -annotate 0 'near real-time (nrt) data' latest_nrt.png

# Create negative color image for nrt data visualization
cp $filename latest.png
convert latest.png -negate latest.png

# Create thumbnail
convert -define png:size=1024x1024 $filename -thumbnail 256x256 -unsharp 0x.5 thumbnail.png

# Delete all .png files older than 60 days
find /web/jsoc/htdocs/doc/data/hmi/harp/harp_nrt/harp.*.png* -type f -atime +60 -exec rm -f {} \;