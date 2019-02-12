#!/bin/csh

# This script starts both the internal- and external-DB export manager wrappers. 
# If the wrappers are currently running, this script first attempts to terminate
# them. 

# uncomment to use the avx cluster
source /SGE2/default/common/settings.csh

# run the export managers
cd /home/jsoc/exports

# run the internal export manager wrapper
if (-f keep_running) then
   rm keep_running
endif

set count=8
while ($count > 0)
   if (! -f keep_running) then
      /home/jsoc/cvs/Development/JSOC/proj/util/scripts/exportmanage.pl -jsocdev procser=jsoc.export_procs &
      break
   endif
   echo 'internal export manager wrapper is still running...will try again in 1 second'
   @ count = ($count - 1)
   sleep 1
end

if ($count == 0) then
   echo 'failure starting internal export-manager wrapper'
   exit 1
endif


# run the external export manager wrapper
if (-f keep_running_web) then
   rm keep_running_web
endif

set count=8
while ($count > 0)
   if (! -f keep_running) then
      /home/jsoc/cvs/Development/JSOC/proj/util/scripts/exportmanage.pl -jsocweb procser=jsoc.export_procs &
      break
   endif
   echo 'external export manager wrapper is still running...will try again in 1 second'
   @ count = ($count - 1)
   sleep 1
end

if ($count == 0) then
   echo 'failure starting external export-manager wrapper'
   exit 1
endif

exit 0
