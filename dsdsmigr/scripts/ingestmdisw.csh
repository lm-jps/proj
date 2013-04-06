#!/bin/tcsh

if ($#argv != 1) then
  echo must supply a single argument giving the number of days to run in each qsub script
  exit 1
endif
set daychunk=$1

setenv JSOCROOT /home/jsoc/cvs/Development/JSOC

set qsubtmp=/tmp27/$USER/qsubtmp
mkdir -p $qsubtmp
set q=j.q

set curdir=$PWD
set sourcedir=/surge40/PAS/D136007/mdi_log/lev0/MDI_log_01d
set jsocdir=/home/jsoc/cvs/Development/JSOC/proj/dsdsmigr/apps
set firstday=1065
set lastday=7013

set day=$firstday
while ($day <= $lastday)
  if (-e $sourcedir/00$day.record.rdb) then
    mkdir $day
    ln -s $sourcedir/00$day.record.rdb $sourcedir/00$day.overview.fits $day/
  endif
  @ day++
end

echo directories created and populated with symlinks

set day1=$firstday
while ($day1 <= $lastday)
  @ day2=$day1 + $daychunk - 1
  if ($day2 > $lastday) set day2=$lastday
  set subfile=subingest.$day1
  echo '#\!/bin/csh' > $subfile
  echo 'setenv PATH' $JSOCROOT'/bin/$JSOC_MACHINE' >> $subfile
  echo 'cd' $PWD >> $subfile
  set day=$day1
  while ($day <= $day2)
    if (-e $sourcedir/00$day.record.rdb) then
      echo ingest_dsds_to_drms in=$curdir/$day map=$jsocdir/nothing.map out=mdi.statusword -sv '>&' log.$day >> $subfile
    endif
    @ day++
  end
  qsub -q $q -e $qsubtmp -o $qsubtmp $subfile
  @ day1 = $day1 + $daychunk
end

echo jobs submitted, start waiting

set njobsrunning = `qstat -r -u $USER | grep "Full jobname:" | grep subingest | wc -l`
while($njobsrunning > 0)
  sleep 60
  set njobsrunning = `qstat -r -u $USER | grep "Full jobname:" | grep subingest | wc -l`
end

echo jobs finished, checking for errors
grep -v "input records found" log.* > errlog
if (-s errlog) then
  echo errors found, check file errlog
  exit 1
else
  echo successful completion
  exit 0
endif
