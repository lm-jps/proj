#

set BATCHSIZE = 500
# set echo 
ingest_test_config_files.csh

set LOG = /tmp/q_multi_log.$$
set SUMMARY = /tmp/q_multi_summary.$$
echo "#! /bin/csh -f" > $LOG

set f = $1
set l = $2
set set = 1
set LEVEL
if ($#argv == 3) then
  set LEVEL = $3
  endif

while ($f <= $l)
  set tl = $l
  @ n = $tl - $f
  @ n = $n + 1
  if ($n > $BATCHSIZE) then
    @ tl = $f + $BATCHSIZE - 1
  endif
  drms_run -L -lifetime 30 q_import_lev0_CIF_from_DSDS.csh $f $tl $LEVEL >& set.$set.$$ &
  echo "(echo -n XXXXXXX set $set ; date ; show_keys ds=hmi_ground.lev0[$f-$tl] -q key=FSN seg=file -p ; echo XXXXXXXX )i" >> $LOG
  @ set = $set + 1
  @ f = $tl + 1
end

wait
 
csh $LOG > $SUMMARY

set SUMMARY_PLACE = /home/jsoc/hmi/ground/q_multi_summary
set NOW = `date +%d_%b_%Y_%H_%M` 
echo $NOW >> $SUMMARY_PLACE/Summary
cat $LOG >> $SUMMARY_PLACE/Summary
cat $SUMMARY >> $SUMMARY_PLACE/Summary
mv $SUMMARY $SUMMARY_PLACE/Summary_$NOW
echo "Please examine  $SUMMARY_PLACE/Summary_$NOW for correctness"

rm -f $LOG
