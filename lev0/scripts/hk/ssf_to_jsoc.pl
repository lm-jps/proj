#!/bin/perl
##############################################################################
# Name:        ssf_to_jsoc.pl  - Send Stanford File To JSOC                  #
# Description: Send the Stanford and Ground file to solserv.lmsal.com where  #
#              files are pickup by script on JSOC machine.                   #
# Execution:   ssf_to_jsoc.pl    on  LMSAL HMIFSW1 machine only              #
##############################################################################
# set Environment Variables
&set_env_variables();

# log status in SendFile2JSOCLog
open(SLF, ">>$sflog");
print SLF `date`;
print SLF "->checking if need to send update\n";
print SLF "->(1)start processing, do cvs update of files in lmsal cvs\n";

# Update CVS on Lockheed CVS on hmifsw1.atc.lmco.com
$log=`cd; csh ulmcvs_for_jsoc.csh`;
print SLF "->(2)\n$log";
sleep(2);

# Check if Different STHA found before last update 
$log=`/usr/bin/chmod 777 $hm/JSOC-slog2;/usr/bin/rm -f $hm/JSOC-slog2;`; 
print SLF "->(3)\n$log";
$log=`cd $hm/$sfd; /usr/local/bin/cvs status STANFORD_TLM_HMI_AIA.txt | grep Working > $hm/JSOC-slog2`;
$log=`/usr/bin/chmod 777 $hm/JSOC-slog2`; 
print SLF "->(4)\n$log";
$diffy=`/usr/local/bin/diff -qs $hm/JSOC-slog1  $hm/JSOC-slog2`;
print SLF "->difference is $diffy";

# Check if Different GTCIDS found before last update 
$log=`/usr/bin/chmod 777 $hm/JSOC-glog2 ;/usr/bin/rm -f $hm/JSOC-glog2`;
print SLF "->(5)\n$log";
$log=`cd $hm/$gfd; /usr/local/bin/cvs status GROUND_to_CODE_ids.txt | grep Working > $hm/JSOC-glog2`;
$log=`/usr/bin/chmod 777 $hm/JSOC-glog2 `;
print SLF "->(6)\n$log";
$diffy2=`/usr/local/bin/diff -qs $hm/JSOC-glog1  $hm/JSOC-glog2`;
print SLF "->(7)\n->differences is $diffy2";

# Based on difference found, if both are different then send files and send message for queue
if ((( index $diffy, "differ") != -1) and (( index $diffy2, "differ") != -1))
{
  ##log status and get latest file version number for latest stanford file
  print SLF "->(8)\n->STANFORD and GROUND files are new.\n";
  $fvid=`cd $hm/$sfd;/usr/local/bin/cvs status STANFORD_TLM_HMI_AIA.txt | grep Working | cut -f2`;
  print SLF "->(9)\n->Stanford file version:$fvid:";
  $log= `cp $hm/JSOC-slog2 $hm/JSOC-slog1`;
  print SLF "->(10)\n$log";

  ##Copy STANFORD file to data import account on Stanford machine
  $log=`cd $hm/$sfd;/usr/local/bin/cvs status`; #log information
  print SLF "->(11)\n$log";
  $log=`cd $hm/$sfd; cp -f STANFORD_TLM_HMI_AIA.txt STANFORD_TLM_HMI_AIA.txt-$fvid `;
  print SLF "->(12)\n$log";
  $log=`cd $hm/$sfd; scp -o user=jim STANFORD_TLM_HMI_AIA.txt solserv.lmsal.com:/home/jim/stanford-files/STANFORD_TLM_HMI_AIA.txt-$fvid`;
  print SLF "->(13)\n$log";

  ##set the file version for STHA file
  $sfvid=$fvid;
  print SLF "->(14)\n";

  ##log status and get latest file version number for latest GTCIDS file
  $fvid=`cd $hm/$gfd;/usr/local/bin/cvs status GROUND_to_CODE_ids.txt | grep Working | cut -f2`;
  print SLF "->(15)\n->GROUND file version number is $fvid";
  $log=`cp $hm/JSOC-glog2 $hm/JSOC-glog1`;
  print SLF "->(16)\n$log";

  ##Copy GTCIDS file to data import account on Stanford machine
  $log=`cd $hm/$gfd;/usr/local/bin/cvs status GROUND_to_CODE_ids.txt`;#show information in log
  print SLF "->(17)\n$log";
  $log=`cd $hm/$gfd;cp -f GROUND_to_CODE_ids.txt GROUND_to_CODE_ids.txt-$fvid `;
  print SLF "->(18)\n$log";
  $log=`cd $hm/$gfd; scp -o user=jim GROUND_to_CODE_ids.txt solserv.lmsal.com:/home/jim/stanford-files/GROUND_to_CODE_ids.txt-$fvid`;

  ##set file version for GTCIDS file
  print SLF "->(19)$log\n";
  $gfvid=$fvid;
  print SLF "->(20)\n";

  ##Create item on queue file to send to solserv which is used by JSOC carl@yeti or production@d00
  $qline=sprintf("%s%s%s %s%s%s","stha=", "STANFORD_TLM_HMI_AIA.txt-",substr($sfvid,0,-1), "gtcids=","GROUND_to_CODE_ids.txt-",substr($gfvid,0,-1));
  print SLF "->(21)\n";
  $log=`ssh  -l jim solserv.lmsal.com \"echo $qline >> /home/jim/stanford-files/UpdateLMQueue.txt\" `;
  print SLF "->(22)\n$log";
  print SLF "->(23)\n->send file updates required.\n";
}
else
{
  print SLF "->no file updates required.\n";
}
close SLF;
sub set_env_variables()
{
 $hm=$ENV{'HOME_DIR'}= '/home/cimilluca';
 $ENV{'STHA_DIR'}= 'SANDBOX/DB_FILES/STANFORD';
 $ENV{'GTCIDS_DIR'}= 'SANDBOX/DB_FILES/GSE';
 $sfd=$ENV{'STHA_DIR'};
 $gfd=$ENV{'GTCIDS_DIR'};
 $sflog="$hm/SendFile2JSOCLog";

 # important variable to get ssh to work
 $ENV{'CVS_RSH'}='/bin/ssh';

}
