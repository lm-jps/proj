#!/bin/csh

# notification list
set $notlist = "arta@sun.stanford.edu jennifer@sun.stanford.edu rock@sun.stanford.edu"

# save path
set cp = echo $PWD

set logFile = "/home/jsoc/sdo/fds/fds.log"
set scriptPath = "/home/jsoc/cvs/JSOC/scripts"
set fdsPermDataPath = "/home/jsoc/sdo/fds"
set fdsDataPath = "/surge/sdo/mocprods"
set fdsSeries = "sdo.moc_fds"

# download files from moc server to scratch disk
set cmdStr = "$scriptPath/dlMOCDataFiles.pl -c $scriptPath/mocDlFdsSpec.txt -s $fdsPermDataPath/mocDlFdsStatus.txt -r $fdsDataPath -t 30 $1"

cd $fdsDataPath
$cmdStr >& $logFile

# ingest into fds data series - don't delete unless ingestion was successful
$scriptPath/fdsIngest.pl $fdsDataPath -r >>& $logFile

# process logfile for email
$scriptPath/fdsNotification.pl -l $logFile -n $notlist

# restore path
cd $cp
