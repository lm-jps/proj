#! /bin/csh -f

#################################################################################
#  What this script does:                          
#
#  This script: 
#  (1) takes the output of show_info,
#  (2) chunks the output into 5-day intervals, and
#  (3) feeds each chunk into ingest_dsds_test_1.
#
#  How to use this script:
#  (1) run show_info -P on a dsds series and send the output into a file, e.g.:
#      show_info -P dsds.mdi__lev1_5__fd_V_01h"[86904-87647]" >> tmp.txt
#
#  (2) run this script, with two input varables -- the file generated from 
#      step (1) and the output series, e.g.:
#      ingest_script_final.csh textfile.txt mdi.blah
#                                             
##################################################################################

if ( $#argv < 1 ) then
  echo ""
  echo "Usage:  ingest_script_final <file> <output_drms_series>"
  echo "Example:ingest_script_final textfile.txt mdi.blah"
  echo "Remember the script must be modified if you want to run ingest_dsds_to_drms with flags!"
  echo ""
  exit 1
endif

set dir = `pwd`
@ LINES = `wc -l $argv[1] | awk '{print $1}'`
@ n = `expr $LINES / 50`
@ file_no = 0

@ i = 1
@ j = 1

while ( $i <= ($n + 1) )
  @ first= $j
  @ last = ($i * 50)
  echo "$first $last"
  sed -n "$first,$last p" $argv[1] > $argv[1]"_"$i
  @ j = $j + 50
  @ i++
end

@ k = 1
while ( $k <= ($n + 1) )
  /home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/ingest_dsds_to_drms in=@$argv[1]_$k out=$argv[2] map=/home/jsoc/cvs/Development/JSOC/proj/dsdsmigr/apps/fd_test.map SCALE_CORRECTIONS=mdi.scale_corrections -M -p
  @ k++
end
