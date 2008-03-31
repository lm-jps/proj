#!/bin/csh

# save path
set cp = `pwd`

# actual bits
set scriptPath = "/home/jsoc/cvs/JSOC/scripts"
set fdsPermDataPath = "/home/jsoc/sdo/fds"
set fdsDataPath = "/surge/sdo/mocprods"
set fdsSeries = "sdo.moc_fds"

# download files from moc server to scratch disk
set cmdStr = "$scriptPath/dlMOCDataFiles.pl -c $scriptPath/mocDlFdsSpec.txt -s $fdsPermDataPath/mocDlFdsStatus.txt -r $fdsDataPath -t 30 $1"

cd $fdsDataPath
$cmdStr

# ingest into fds data series - don't delete unless ingestion was successful
$scriptPath/fdsIngest.pl $fdsDataPath -r

# restore path
cd $cp
