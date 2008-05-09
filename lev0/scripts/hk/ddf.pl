#!/usr/bin/perl
##############################################################################
# Name:        ddf.pl - decode day files                                     #
#              CRON to Get  hmi,sdo,aia dayfiles from DRMS & send to for     #
#              decoding keywords and then write keywords to DRMS hk data     #
#              series by apid. This script sets values to process dayfiles   #
#              for either hmi,aia, or sdo  dayfiles for today. Sets up       #
#              environment variable for hmi,aia, or sdo input dayfile series.#
#              Can turn on debug flag using DF_GDFDRMS_DEBUG and this is     #
#              passed to gdfdrms.pl. Can turn on  report flag which will     #
#              turn on reporting in gdfdrms.pl. Key item to set is           #
#              DF_SERIES_NAME correctly. Set to <proj name>_ground.hk_dayfile#
#              currently. This day file series name is used as               #
#              input  data series to gather dayfiles need to process.        #
#              APIDs to process determined by input arguments to gdfdrms.pl. #
#              Currently input args to gdfdrms.pl are for today, for apids   #
#              set in apidlist argument, and for source set in src argument. #
#              Log file is set using $logfile and is passed to gdfdrms.pl.   #
#              Setup and create apidlist files:ddf_apid_list_hmi_egse,etc    #
#              This script can be run at command line too.                   #
# Execution:   (1)Run option to process dayfiles for today:                  #
#                                                                            #
#                 perl hmi_dfd.pl <project name of day file series> <source> #
#                                                                            #
#              (2)The run options are show in help listing:                  #
#                                                                            #
#                 perl dfd.pl  -h                                            #
#                                                                            #
# Examples Execution(Possible cases as of today):                            #
#              perl ddf.pl hmi egsefm                                        #
#              perl ddf.pl hmi egseem                                        #
#              perl ddf.pl hmi hsb                                           #
#              perl ddf.pl hmi moc                                           #
#              perl ddf.pl aia egsefm                                        #
#              perl ddf.pl aia egseem                                        #
#              perl ddf.pl aia hsb                                           #
#              perl ddf.pl aia moc                                           #
#              perl ddf.pl sdo moc                                           #
#              perl ddf.pl sdo rtmon                                         #
# Limitation:  Setup required environment variables at top of file.          #
#              Must have input dayfile series(i.e.,hmi.hk_dayfile,           #
#              sdo.hk_dayfile, etc). Must have apidlist file.                #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on May,7, 2008    #
##############################################################################
# main program                                                               #
##############################################################################
#set environment variables for HMI RELATED DECODE DAYFILE PROCESSING

#(0)get source argument
&check_agruments();

#(1)set paths
$hm=$ENV{'HOME'};
$script_dir=$ENV{'HK_SCRIPT_DIR'}="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";

#(2)cron setting
$ENV{'MAILTO'}="";

#(3)set debug mode
$dflg=$ENV{'DF_GDFDRMS_DEBUG'}=1;

#(4)setup log file for with name based on input argument(hmi,aia,sdo), 
$logfile=$ENV{'HK_DF_LOGFILE'}="$script_dir/log-$dspnm-ddf";

####open log and record log info####
if ($dflg) {open(LF,">>$logfile") || die "ERROR in ddf.pl:Can't Open <$logfile>: $!\n; exit;";}
if ($dflg) {print LF `date`;}
if ($dflg) {print LF "--->ddf.pl:debug:log file is <$ENV{'HK_DF_LOGFILE'}>\n";}
if ($dflg) {print LF "--->ddf.pl:debug:source argument is <$src>\n";}

#(5)set report mode on or off. if turn on, set report name want to use below.Creates gzipped report.
$rptflg=$ENV{'DF_REPORT_FLAG'}=0;

#(6)set report naming convention for using input argument(hmi,aia,sdo) if reporting is turned on
$ENV{'DF_PROJECT_NAME_FOR_REPORT'}="$dspnm";
$ENV{'DF_DATA_TYPE_NAME_FOR_REPORT'}="lev0";

#(7)set input data series to use  to check for dayfiles to decode for based on input arg(hmi,sdo,aia)
$dsnm=$ENV{'DF_SERIES_NAME'}="$dspnm\_ground.hk_dayfile"; #test LMSAL formatted dayfiles
##used to test hsb dayfiles input files from drms
#$dsnm=$ENV{'DF_SERIES_NAME'}="su_carl.hk_dayfile";   #test hsb formatted dayfiles

####log info####
if ($dflg) {print LF "--->ddf.pl:debug:input dayfile is <$dsnm>\n";}

#(8)get apid list to use
$list=&get_apid_list();

#(9)set up command for  dayfile series with source egse
# execute and retrieve dayfile for today and decode keywords.
$command="perl $script_dir/gdfdrms.pl  apidlist=$script_dir/$list  src=$src";

#Note:Used to test with known data. do:perl ddf.pl hmi egsefm
#$command="perl $script_dir/gdfdrms.pl  apidlist=$script_dir/$list  start=20080419 end=20080419 src=$src";

####log info####
if ($dflg) {print LF "--->ddf.pl:debug:Running:Command running:<$command>\n";}

#(10)execute to decode dayfiles using APIDs in apidlist
$log=`$command`;

####log info####
if ($dflg) {print LF "--->ddf.pl:debug:Finished. Command Log is <$log>\n";}
if ($dflg) {print LF `date`;}
if ($dflg) {close(LF);}


#############################################################################
# get_apid_list()
#############################################################################
sub get_apid_list
{
  #as convention keep apid list in file using this format so can easily add or delete apids
  return( "ddf_apid_list_$dspnm\_$src");
}

##############################################################################
# check_arguments()                                                          #
##############################################################################
sub check_agruments()
{
  #data series project name(i.e.,sdo,hmi,aia) for hk_dayfile data series
  $dspnm=$ARGV[0]; #data series project name(i.e.,sdo,hmi,aia)
  #source index in hk_dayfile data series
  $src=$ARGV[1];
  #check arguments
  if ($#ARGV <= 0 || $#ARGV > 2 )
  {
    print "Usage: perl ddf.pl <project name of input dayfile> <dayfile source>\nwhere project name is sdo,aia,hmi.\nwhere source is either:hsb,moc,egsefm.\n";
    exit;
  }
  elsif ("-h" eq substr($ARGV[0],0,2) )
  {
     print "Usage: perl ddf.pl <project name of input dayfile> <dayfile source>\nwhere project name is hmi,aia,sdo.\nwhere source is either:hsb,moc,egsefm.\n";
     exit;
  }
  elsif("hmi" eq substr($ARGV[0],0,3) or "aia" eq substr($ARGV[0],0,3) or "sdo" eq  substr($ARGV[0],0,3))
  {
    if("hsb" eq substr($ARGV[1],0,3) or "egsefm" eq substr($ARGV[1],0,6) or "egseem" eq  substr($ARGV[1],0,6) or "moc" eq  substr($ARGV[1],0,3) or "rtmon" eq  substr($ARGV[1],0,5) or "lmcmd" eq  substr($ARGV[1],0,5)   )
    {
      #print "okay";
    }
    else
    {
      print "ERROR: Entered incorrect source value. Use egsefm,hsb, or moc.\n";
      print "Usage: perl ddf.pl <project name of input dayfile> <dayfile source>\nwhere project name is hmi,aia,sdo.\nwhere source is either:hsb,moc,egsefm.\n";
      exit;
    }
  }
  else
  {
     print "ERROR: Entered incorrect project name for data series value. Use sdo,hmi or aia\n";
     print "Usage: perl ddf.pl <project name of input dayfile> <dayfile source>\nwhere project name is hmi,aia,sdo.\nwhere source is either:hsb,moc,egsefm.\n";
     exit;
  }
}
