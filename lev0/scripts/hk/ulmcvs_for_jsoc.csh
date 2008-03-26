#!/bin/csh
# Name: ulmcvs_for_jsoc.csh  Update LMSAL CVS for JSOC
# Description: Runs only on HMIFSW1 machine at LMSAL to update LMSAL CVS
# Date: 3-27-2008 moved to JSOC
source $HOME/.cshrc;
cd;
cd SANDBOX/DB_FILES;
/usr/local/bin/make update >> /home/cimilluca/SendFile2JSOCLog;
