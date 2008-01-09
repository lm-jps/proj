#!/bin/csh

# save path
set cp = `pwd`

# actual bits
set scriptPath = "/home/jsoc/cvs/JSOC/scripts"
set lzpDataPath = "/home/jsoc/sdo/lzp"

# Create the spec file on-the-fly.  Specify all files that could possibly exist
# in the last three months.  dlMOCDataFiles.pl will look for those files only.
set MONTHS = (January February March April May June July August September October November December)

# get system times
set yr = `date +%Y`
set mo = `date +%m`
set doy = `date +%j`
@ pyr = $yr - 1

set fullspec = ""
set pspec = ""
set cspec = ""
set filesuffix = "[0-9][0-9][0-9]_[0-9][0-9]\.hkt\S*"

set pbegin = "none"
set pend = 366
set cbegin = "none"
set cend = `date +%j`

# To test last case below
# set mo = "11"
# set cend = `date -dNovember-18-2008 +%j`

if ($mo == "01") then 
  # Add last 3 months of the previous year
  set pbegin = `date -dOctober-01-$pyr +%j`
  set cbegin = `date -dJanuary-01-$yr +%j`
else if ($mo == "02") then
  # Add last 2 months of the previous year
  set pbegin = `date -dNovember-01-$pyr +%j`
  set cbegin = `date -dJanuary-01-$yr +%j`
else if ($mo == "03") then
  # Add last 1 month of the previous year
  set pbegin = `date -dDecember-01-$pyr +%j`
  set cbegin = `date -dJanuary-01-$yr +%j`
else
  # Add last 3 months of the current year
  @ stmo = $mo - 3
  set mobegin = $MONTHS[$stmo]
  set cbegin = `date -d${mobegin}-01-$yr +%j`
endif

if ($pbegin != "none") then 
  set pspec = "${pyr}_|s[$pbegin-$pend]::000[1-9]_${pyr}_$filesuffix\n${pyr}_|s[$pbegin-$pend]::00[1-5][0-9]_${pyr}_$filesuffix\n${pyr}_|s[$pbegin-$pend]::006[0-3]_${pyr}_$filesuffix\n${pyr}_|s[$pbegin-$pend]::0129_${pyr}_$filesuffix\n" 
endif 

# Add current year days
if ($cbegin != "none") then
  set cspec = "${yr}_|s[$cbegin-$cend]::000[1-9]_${yr}_$filesuffix\n${yr}_|s[$cbegin-$cend]::00[1-5][0-9]_${yr}_$filesuffix\n${yr}_|s[$cbegin-$cend]::006[0-3]_${yr}_$filesuffix\n${yr}_|s[$cbegin-$cend]::0129_${yr}_$filesuffix\n"
else
  # error
  echo "Error getting today's date."
endif

set fullspec = "$pspec$cspec"

set noglob
printf "root:moc/lzp\n\n" > $lzpDataPath/mocDlLzpSpec.txt
printf $fullspec >> $lzpDataPath/mocDlLzpSpec.txt

set cmdStr = "$scriptPath/dlMOCDataFiles.pl -c $lzpDataPath/mocDlLzpSpec.txt -s $lzpDataPath/mocDlLzpStatus.txt -r $lzpDataPath/MOCFiles/ -t 120 $1"

cd $lzpDataPath
perl $cmdStr
# printf "%s" $cmdStr

# restore path
cd $cp
