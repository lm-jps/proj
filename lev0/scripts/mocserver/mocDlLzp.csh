#!/bin/csh

# save path
set cp = `pwd`

# actual bits
set scriptPath = "/home/jsoc/cvs/JSOC/scripts"
set lzpDataPath = "/home/jsoc/sdo/lzp"
set cmdStr = "$scriptPath/dlMOCDataFiles.pl -c $scriptPath/mocDlLzpSpec.txt -s $lzpDataPath/mocDlLzpStatus.txt -r $lzpDataPath/MOCFiles/ -t 120 $1"

cd $lzpDataPath
perl $cmdStr

# restore path
cd $cp
