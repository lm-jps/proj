#!/usr/bin/perl
##############################################################################
# Name:        gdfdrms.pl - Get  dayfiles from DRMS & send to decode_dayfile #
# Description: Used to gather Level 0 high speed bus,LMSAL,and SDO formatted #
#              dayfiles from DRMS. The script can gather a list of dayfiles  #
#              to process based on arguments used for script and then sends  #
#              dayfile to the decode_dayfile.c executable which writes       #
#              keywords to DRMS. There is  few different ways to gather up   #
#              day files. Based on apid, dates, and source. To run this the  #
#              decode_dayfile.c need to have environment variable setup with #
#              project name(hmi,aia,etc), data type name(lev0,etc), etc.     #
#              Switch on a small amount of debug using DF_GDFDRMS_DEBUG set  #
#              to 1 below.  DF_REPORT_FLAG produces a report file if set     #
#              to 1. DF_PROJECT_NAME_FOR_REPORT and                          #
#              DF_DATA_TYPE_NAME_FOR_REPORT are used to create reportname    #
#              The data series project name and data type name used by       #
#              decode_dayfile is set in SOURCE_ENV_FOR_HK_DAYFILE_DECODE file#
#              Use -h to see more information on running. The script was     #
#              create to run as cron job or a command line.                  #
# Execution:   (1)The run options are show in help listing:                  #
#                   gdfdrms  -h                                              #
#              (2)The execution options are show when run "gdfdrms  -h"      #
# Limitation:  Setup required environment variables at top of file based on  #
#              where running. Script works only using hsb dayfile file       #
#              and LMSAL formats that are in DRMS. Here are examples input   #
#              day file:                                                     #
#              -For src=hsb   : hsb_0445_2008_04_15_10_41_00.hkt             #
#              -For src=egsefm: 200804019.0x001d                             #
#              -For src=rtmon : 200804019.0x001d                             #
#              -For src=moc   :0445_2008_04_15_10_41_01.hkt                  #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on May,2, 2008    #
##############################################################################
# main program                                                               #
##############################################################################

#(1)common setting for all environmnents
$ENV{'SUMSERVER'}="k1";
$hm=$ENV{'HOME'};
$ENV{'MAILTO'}="";
$exec_dir=$ENV{'DF_EXEC_PATH'}="$hm/cvs/JSOC/bin/linux_x86_64";
$script_dir=$ENV{'HK_SCRIPT_DIR'}="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$exec_dir";

#(2) set debug flag 1 to turn on and 0 to turn off
$dflg=$ENV{'DF_GDFDRMS_DEBUG'}=0;

#(3) set report flag to turn on and off report creation and set report name 
$rptflg=$ENV{'DF_REPORT_FLAG'};
$pjn=$ENV{'DF_PROJECT_NAME_FOR_REPORT'};
if($pjn eq "")
{
  #default
  $pjn=$ENV{'DF_PROJECT_NAME_FOR_REPORT'}="su_carl";
}
$dtn=$ENV{'DF_DATA_TYPE_NAME_FOR_REPORT'};
if($dtn == 0)
{
  #default
  $dtn=$ENV{'DF_DATA_TYPE_NAME_FOR_REPORT'}="lev0";
}

#(4) check for any arguments passed in command
&check_arguments();

#(5) set choices for input dayfile
$dsnm_sdo=$ENV{'DF_SERIES_NAME'}="sdo.hk_dayfile"; 
$dsnm_hmi=$ENV{'DF_SERIES_NAME'}="hmi.hk_dayfile"; 
$dsnm_aia=$ENV{'DF_SERIES_NAME'}="aia.hk_dayfile"; 

#(6) setup log file -use logfile setting set in hmi_dfd.pl,aia_dfd.pl, sdo_dfd.pl or use setting below
$logfile=$ENV{'HK_DF_LOGFILE'};
if($logfile eq "")
{
   $logfile="$script_dir/log-gdfdrms";#send log here when run at command line
}

#(7) setup log file -use logfile setting set in hmi_dfd.pl,aia_dfd.pl, sdo_dfd.pl or use setting below
$logfile=$ENV{'HK_DF_LOGFILE'};
if($logfile eq "")
{
   $logfile="$script_dir/log-gdfdrms";#send log here when run at command line
}

#(8) open log file
open(LF,">>$logfile") || die "ERROR in gdfdrms.pl:Can't Open-0 <$logfile>: $!\n; exit;";
print LF `date`;

if($dsnm != 0)
{
  print LF "-->(0)Getting dayfiles from DRMS from data series <$dsnm>\n";
}


#(9) get list of apids to do 
&get_apid_to_do();

#(10)get list of dates to do 
&get_date_to_do();

#(11)go through day files in directory looking for all files for given apid and date range
&get_list_from_drms();

#(12)decode keywords from list of  day files and write keywords to DRMS data series
&decode_keywords_for_dayfiles();
 
#(13)close logfile
print LF `date`;
close LF;


#############################################################################
# subroutine check arguments and set flags                                  #
#############################################################################
sub check_arguments()
{
  $help_flg= "0";
  $apid_list_flg="0";
  $date_range_flg="0";
  $view_flg="0";
  
  if ($#ARGV < 0)
  {
    $help_flg="1";
  }
  if ($#ARGV >= 0)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-H")
    {
      $help_flg = "1";
    }
    elsif (substr($ARGV[0],0,2) eq "-v" )
    {
      #view only dayfiles want to load
      $view_flg="1";

      if (substr($ARGV[1],0,9) eq "apidlist=" )
      {
        #use file to get list of apids to create map files for
        $apid_list_flg = "1";
      }
      elsif (substr($ARGV[1],0,5) eq "apid=" )
      {
        #use apid following -a to create map files
        #push decimal character value in this format dddd. Example: -a 0001
        $apid_list_flg = "2";
      }
 
    }
    elsif (substr($ARGV[0],0,9) eq "apidlist=" )
    {
      #use file to get list of apids to create map files for
      $apid_list_flg = "1";
    }
    elsif (substr($ARGV[0],0,5) eq "apid=" )
    {
      #use apid following -a to create map files 
      #push decimal character value in this format dddd. Example: -a 0001
      $apid_list_flg = "2";
    }
  }
  if ($#ARGV >= 1 && $view_flg eq "0")
  {
    if ( (substr($ARGV[1],0,6) eq "start=") and (substr($ARGV[2],0,4) eq "end=") )
    {
      #use file to get list of apids to create map files for
      $date_range_flg = "1";
    }
    elsif (substr($ARGV[1],0,4) eq "src=" )
    {
      $source=substr($ARGV[1],4,6);
      if($dflg) {print LF "--->gdfdrms.pl:DEBUG:src is $source\n";}
    }
    else
    {
      print  LF "gdfdrms.pl:ERROR: Did not use src argument. Exiting script\n";
      print   "gdfdrms.pl:ERROR: Did not use src argument. Exiting script\n";
      exit;
    }
  }
  elsif ($#ARGV >= 2 && $view_flg eq "1")
  {
    if ( (substr($ARGV[2],0,6) eq "start=") and (substr($ARGV[3],0,4) eq "end=") )
    {
      #use file to get list of apids to create map files for
      $date_range_flg = "1";
    }
    elsif (substr($ARGV[2],0,4) eq "src=" )
    {
      $source=substr($ARGV[2],4,6);
      if($dflg) {print LF "gdfdrms.pl:DEBUG:src is $source\n";}
    }
    else
    {
      print LF "WARNING: Did not use src argument. Exiting script\n";
      exit;
    }
 }


  if ($#ARGV >= 3 && $view_flg eq "0")
  {
    if ( substr($ARGV[3],0,4) eq "src=") 
    {
      $source=substr($ARGV[3],4,6);
      if($dflg) {print LF "gdfdrms.pl:DEBUG:src is $source\n";}
    }
    else
    {
       print LF "gdfdrms.pl:ERROR: argument 3 is not correct: $ARGV[3]. Existing. Enter correct value. \n";
       print  "gdfdrms.pl:ERROR: argument 3 is not correct: $ARGV[3]. Existing. Enter correct value. \n";
       exit;
    }
  }
  elsif ($#ARGV >= 4 && $view_flg eq "1")
  {
    if ( substr($ARGV[4],0,4) eq "src=") 
    {
      $source=substr($ARGV[4],4,6);
      if($dflg) {print LF "gdfdrms.pl:DEBUG:src is $source\n";}
    }
    else
    {
       print  "gdfdrms.pl:ERROR: Did not use src argument. Exiting script\n";
       print LF "gdfdrms.pl:ERROR: argument 3 is not correct: $ARGV[4]\n";
       print LF "This is the source argument.Example: src=hsb. \n";
       print LF "Existing. Enter correct value. \n";
       exit;
    }
  }
  if($dflg) {print LF "gdfdrms.pl:DEBUG:SOURCE passed as argument is $source\n";}

  if ( $help_flg eq "1")
  {
     &show_help_info;
  }
}


#############################################################################
# subroutine show_help_info: show help information                          #
#############################################################################
sub show_help_info
{
  print "Help Listing\n";
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Decode Day Files using apidfile option will decode keywords todays(localtime()) dayfiles based on apids contained in file:\n
             gdfdrms.pl  apidlist=<filename containing APID List to do> src=<source of dayfile>\n";
  print "             Example: gdfdrms.pl apidlist=./df_apid_list_day_files_egsefm src=egsefm\n\n";
  print "(1b)Decode Day Files using apid option will decode keywords for today(localtime()) dayfiles for given apid  and src value:\n
             gdfdrms.pl apid=<apid-value in decimal format>  src=<source of dayfile>\n";
  print "             Example: gdfdrms.pl apid=445  src=hsb\n\n";
  print "(1c)Decode Day Files using date range with apidfile option to decode keywords for all dayfile with given apid, date, and src:\n
             gdfdrms.pl apidlist=<file> start=<yyyymmdd> end=< yyyymmdd> src=<source of dayfile>\n";
  print "             Example: gdfdrms.pl apidlist=./df_apid_list_day_files_hsb start=20070216 end=20070218 src=hsb\n\n";
  print "(1d)Decode Day Files using date range with apid option to decode keywords for dayfile with given apids,dates, and src:\n
             gdfdrms.pl  apid=<apid in decimal>  start=<yyyymmdd> end=<yyyymmdd>  src=<source of dayfile>\n";
  print "             Example: gdfdrms.pl apid=445 start=20070216  end=20070218 src=hsb\n\n";
  print "(1e)Get Help Information:\n
              gdfdrms.pl -h  or  gdfdrms.pl -help\n\n";
  print "(1f)View what going to save in data series by adding -v as first argument. This should work for 1a,1b,1c,and 1d cases.:\n
             gdfdrms -v apidlist=<file> start=<yyyymmdd> end=< yyyymmdd>  src=<source of dayfile>\n";
  print "             Example: gdfdrms.pl -v apidlist=./df_apid_list_day_files_hsb start=20070216 end=20070218  src=hsb\n\n";
  print "*Note:Propose using option for cron jobs(1a): gdfdrms.pl apidlist=<file containing APID list>  src=<Source of Dayfile>\n";
  print "(2) Requires setup and dayfile in hk dayfile data series(i.e.sdo.hk_dayfile,hmi.hk_dayfile,etc.)\n";
  print "(2a)Set choices to possible dayfile in step #5 in script: #(5) set choices for input dayfile.\n";
  print "(3) Requires setup of apid list in file when use \"apidlist\" option\n";
  print "(3a)Enter values in file in decimal format(i.e.,0001,0015,0021,0445,0475).\n";
  print "(3b)Example format and file of apid list is located in this current file: ./df_apid_list_day_files\n";
  print "(4) Requires setup of data series name and apid list in file when use \"dsname\" argument.\n";
  print "(4a)Enter values in file in decimal format for APID value and this data series name(i.e.,0445  hmi_ground.hk_dayfile).\n";
  print "(4b)Example format and file of apid list is located in this current files: df_apid_list_rtmon,df_apid_list_hsb,etc.\n";
  print "(4c)Create a file for different types of data to save files in correct series names.\n";
  print "(4d)For example, create df_apid_ds_list_for_rtmon to direct which series to save decoded keywords.View example.\n";
  print "(4e)Other Examples files:df_apid_ds_list_moc, df_apid_ds_list_hsb, df_apid_ds_list_hsb, df_apid_ds_list_egsefm, etc.\n";
  print "(5)****Limitation****: a)Works only on hsb,rtmon,egsefm,and moc dayfile formats\n";
  print "                       b)Enter arguments in specified order when running script is required.\n";
  exit;
}

#############################################################################
# subroutine get_apid_to_do: gets list of apid to create maps files for     #
#############################################################################
sub get_apid_to_do
{
  # create list of files to process using apid values in file
  if ($apid_list_flg eq "1")
  {

    if ($view_flg eq "0")
    {
      $fn= substr($ARGV[0],9);
    }
    else
    {
      $fn= substr($ARGV[1],9);
    }
    open(FILE_APID_LIST, "$fn") || die "ERROR in gdfdrms.pl:Can't Open-2: $fn file: $!. Need to create apidlist file.\n";
    while (<FILE_APID_LIST>)
    {
      push(@all_apids, int $_) ;
    }
    close FILE_APID_LIST ;
    print  LF "-->(1)Doing decimal formatted apids: @all_apids\n";
  }
  elsif ($apid_list_flg eq "2")
  {
    #push decimal character value in this format dddd.Example 0001
    if ($view_flg eq "0")
    {
      push(@all_apids, substr($ARGV[0],5));
    }
    else
    {
      push(@all_apids, substr($ARGV[1],5));

    }
    print  LF "-->(1)Doing decimal formatted apids: @all_apids\n";
  }
  else 
  {
    print LF "WARNING: Not valid apid list flag value\n";
  }
}


#############################################################################
# subroutine get_date_to_do: gets date range to do                          #
#############################################################################
sub get_date_to_do
{
  # get start and end range of dates
  if ($date_range_flg eq "1")
  {
    if ($view_flg eq "0")
    {
      $date_r1= substr($ARGV[1],6,8);
      $date_r2= substr($ARGV[2],4,8);
    }
    else
    {
      $date_r1= substr($ARGV[2],6,8);
      $date_r2= substr($ARGV[3],4,8);
    }
    #create list of dates based on range
    &get_date_list_to_do;
  }
  else
  {
     #do case with no entry for start and end time
     #for this case can make do today
     ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
     $year = 1900 + $yearOffset;
     $month= $monthOffset + 1;
     if($dflg == 1) {print LF "gdfdrms.pl:DEBUG: year this $year month is $month day is $dayOfMonth\n";}
     #create todays date format and push on list of dates to do
     $new_date= sprintf("%4d.%02.2d.%02.2d", $year,$month,$dayOfMonth);
     push(@all_dates, $new_date);
     #set values to display in log
     $date_r2=$date_r1=sprintf("%4d%02.2d%02.2d", $year,$month,$dayOfMonth);
  }
}


#############################################
# get date list to do                       #
#############################################
sub get_date_list_to_do()
{
  #parse values in date range
  $start_yr= int substr($date_r1,0,4);
  $start_mo= int substr($date_r1,4,2);
  $start_day= int substr($date_r1,6,2);
  $end_yr= int substr($date_r2,0,4);
  $end_mo= int substr($date_r2,4,2);
  $end_day= int substr($date_r2, 6,2);

  #create list of dates  to do
  if ($start_yr == $end_yr and  $start_mo == $end_mo and $start_day == $end_day)
  {
    #only need start day
    $new_date= sprintf("%4d.%02.2d.%02.2d", $start_yr,$start_mo,$start_day);
    push(@all_dates, $new_date);
  }
  elsif ($start_yr != $end_yr or $start_mo != $end_mo or  $start_day != $end_day )
  { 
    $loop_yr=$start_yr;
    $loop_mo= $start_mo;
    $loop_day=$start_day;
    $loop_yr_mo_day=sprintf("%4d%02.2d%02.2d", $loop_yr,$loop_mo,$loop_day); #setup loop start date(i.e.,20080901)
    $end_yr_mo_day=sprintf("%4d%02.2d%02.2d", $end_yr,$end_mo,$end_day);     #setup loop end date(i.e.,20081231)

    # loop thru all dates requested thru arguments and add date to list of dates to check if file in drms dayfile data series
    while ( $loop_yr_mo_day  <=  $end_yr_mo_day)
    {
      #check if month is last day-
      #if so set day to one or 
      #if end of year day then reset month and day  or
      #if not above cases, just increment day by one
      if(&check_lastday_in_month) 
      {
        #if end of month
        if($dflg) {print LF "gdfdrms.pl:DEBUG:pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        $new_date= sprintf("%4d.%02.2d.%02.2d", $loop_yr,$loop_mo,$loop_day);
        push(@all_dates, $new_date);
        $loop_day=1;
        $loop_mo++;
      }
      elsif ($loop_mo == 12 && $loop_day == 31)
      {
         #if end of year
        if($dflg) {print LF "gdfdrms.pl:DEBUG:pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        $new_date= sprintf("%4d.%02.2d.%02.2d", $loop_yr,$loop_mo,$loop_day);
        push(@all_dates, $new_date);
        $loop_day=1; 
        $loop_mo=1;
        $loop_yr++;
      }
      else
      {
        #if not end of year or month then just increment day
        if($dflg) {print LF "gdfdrms.pl:DEBUG:pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        $new_date= sprintf("%4d.%02.2d.%02.2d", $loop_yr,$loop_mo,$loop_day);
        push(@all_dates, $new_date);
        $loop_day++;
      }
      #set up next loop date number
      $loop_yr_mo_day=sprintf("%4d%02.2d%02.2d", $loop_yr,$loop_mo,$loop_day);
    }
    if($dflg) {print LF "gdfdrms.pl:DEBUG:list of dates to process:\n @all_dates\n"};


  }
  else 
  {
    print LF "gdfdrms.pl:ERROR: Probably bad entry for dates\n";
    print  "gdfdrms.pl:ERROR: Probably bad entry for dates\n";
  }
}

###############################################################################
# subroutine get_list_from_drms: get list of filenames to decode keywords for #
###############################################################################
sub check_lastday_in_month()
{
     if( (($loop_mo == 4 || $loop_mo == 6 || $loop_mo == 9 || $loop_mo == 11 ) && $loop_day == 30) || (($loop_mo == 1 || $loop_mo == 3 || $loop_mo == 5 || $loop_mo == 7 || $loop_mo == 8|| $loop_mo == 10) && $loop_day == 31))
     {
        return 1;
     }
     elsif ( $loop_mo == 2 && (($loop_day == 29 && $loop_yr == 2008) || ($loop_day == 29 && $loop_yr == 2012) || ($loop_day == 29 && $loop_yr == 2016) || ($loop_day == 29 && $loop_yr == 2020) || ($loop_day == 29 && $loop_yr == 2024)))
     {
        return 1;
     }
     elsif ( $loop_mo == 2 && $loop_day == 28 && $loop_yr != 2008 &&  $loop_yr != 2012 && $loop_yr != 2016  && $loop_yr != 2020 && $loop_yr != 2024 )
     {
        return 1;
     }
     else
     {
        return 0;
     }
}


###############################################################################
# subroutine get_list_from_drms: get list of filenames to decode keywords for #
###############################################################################
sub get_list_from_drms()
{
  #(1)Show and get all files loaded for that day
  if ($date_range_flg eq "1" or $date_range_flg eq "0")
  {
    print LF "-->(2)current dayfile date ranges are:$date_r1 to $date_r2\n";
     foreach $d (@all_dates)
     {
       #Get Dayfiles from DRMS for this date
       print LF "-->(3)Getting dayfiles from DRMS for date <$d>\n";
       #create list of apid to process in apid string variable
       $apid_str="";
       foreach $a (@all_apids)
       {
         #remove cr
         $a=~s/\n//;

         if ($apid_str eq "")
         {
           $apid_str=sprintf("%s", $a);
         }
         else
        {
          $apid_str=sprintf("%s,%s",$apid_str, $a);
        }
       }
       print LF "-->(4)Getting dayfiles from DRMS for date <$d> and apids <$apid_str>\n";

       #remove comma at end of apid string
       $apid_str=~ s/,+$//;

       # to test single date uncomment out
       #$d="2008.04.19";#for apid=19 and src=egse
       # create command to get all directory's with dayfile from given date,apid and source

       # append to date to have correct index to hk_dayfile series
       $append_to_d="_00:00:00_TAI";

       # get source of dayfile to use based on apid and source values
       if (($source eq "moc" or  $source eq "rtmon") and (int $apid_str > 96) and (int $apid_str < 400))
       {
          #if source is moc or rtmon always grab dayfile from sdo.hk_dayfile for apid between 96-399
          $dsnm=$dsnm_sdo;
       }
       elsif (($source eq "hsb" or source eq "hsb") and (int $apid_str > 32) and (int $apid_str < 64))
       {
          #if source is hsb or egsefm and apids are aia apids, then always grab dayfile from aia.hk_dayfile
          $dsnm=$dsnm_aia;
       }
       elsif (($source eq "hsb" or source eq "egse") and (int $apid_str > 500) and (int $apid_str < 600))
       {
          #if source is hsb or egsefm and apids are aia apids, then always grab dayfile from aia.hk_dayfile
          $dsnm=$dsnm_aia;
       }
       elsif (($source eq "hsb" or $source eq "egse")  and (int $apid_str > 0) and (int $apid_str < 33))
       {
          #if source is hsb or egsefm and apids are hmi apids, then always grab dayfile from hmi.hk_dayfile
          $dsnm=$dsnm_hmi;
       }
       elsif (($source eq "hsb" or $source eq "egse") and (int $apid_str > 400) and (int $apid_str < 500))
       {
          #if source is hsb or egsefm and apids are hmi apids, then always grab dayfile from hmi.hk_dayfile
          $dsnm=$dsnm_hmi;
       }
       else
       {
          print "WARNING:gdfdrms.pl:unexpected else-did not find source and apid case for setting input dayfile name!Exiting script.\n";
          exit;
       }
       if($dflg) {print "gdfdrms.pl:dsnm value is $dsnm\n"};
 
       $command=sprintf("show_keys ds=%s[%s%s][%s][%s] -P",$dsnm,$d,$append_to_d,$apid_str,$source);
       if($dflg) {print LF "gdfdrms.pl:DEBUG:command value:$command returned:$results\n"};
       $results=`$command`;
       $return=substr($results,0,5);
       if ("SUDIR" ne $return)
       {
          if($dflg) {print LF "--->gdfdrms.pl:DEBUG:warning(not an error) no record found: $command returned: $results\n"};
          #skip to next 
          print LF "-->(5)List of dayfiles for dates and apids are:\n@there_dir_list\n";
          next;
       }

       # more error checking
       $exit_value  = $? >> 8;
       $signal_num  = $? & 127;
       $dumped_core = $? & 128;
       if($dflg)  {print LF "gdfdrms.pl:DEBUG:return code from command is <$exit_value>  <$signal_num>  <$dump_core>\n";}
       if($dflg == 1) {print LF "gdfdrms.pl:DEBUG:COMMAND IS: <$command> \nresults:\n<$results>\n";}

       #returns a list of directories to get files
       if($dflg) {print LF "gdfdrms.pl:DEBUG:directory for dayfiles for $d from drms:$results;"}
       push(@dfdate, $results);

       #remove all change line, spaces, etc and put -- delimiter like 
       $results =~ s/#//g;
       $results =~ s/\/r|\n/:/g;
       $results =~ s/ //g;
       if($dflg) {print LF "gdfdrms.pl:DEBUG:directory for dayfiles cleaned:$results\n";}

       # split line into directory locations of dayfiles 
       @s_results=split(':',$results );
       foreach $item (@s_results)
       {
         if (substr($item,0,8) eq "suidback"  )
         {
           if($dflg) {print LF "gdfdrms.pl:DEBUG:Warning:no valid directory for dayfile\n";}
         }
         elsif  (substr($item,0,5) eq "SUDIR"  )
         {
           if($dflg) {print LF "gdfdrms.pl:DEBUG:Warning:no valid directory for dayfile\n";}
         }
         else
         {  
           #loop thru all files in each item directory 
           #push full path to file name on  all_there list.
           #create loop of dates ranges to do to add to allthere list

           #open directory where file is and add files name to there list.
           if($dflg) { print LF "gdfdrms.pl:DEBUG:opening directory : $item \n";}
           opendir(DFDIR, $item) || die "ERROR in gdfdrms.pl:Can't open-3 $item:$!\n"; #open subdirectory

           # read in files
           @df = readdir(DFDIR); #get a list of directory contents
           if($dflg) { print LF "gdfdrms.pl:DEBUG:files in directory are : @df \n";}

           # add file name and directory path to there_dir_list
           foreach $dird (@df)
           { 
             if ($dird eq "." or $dird eq "..")
             {
                # skip . or ..  files found
                 if($dflg == 1) {print LF "gdfdrms.pl:DEBUG:skip . and .. files found.\n";}
                next;
             }
             elsif( index($dird, 'xml') > 0 )
             {
                # skip  xml file found
                 if($dflg == 1) {print LF "gdfdrms.pl:DEBUG:skip xml file found \n";}
                next;
             }
             else
             {
               # add directory and filename path to there_dir_list
               if($dflg == 1) {print LF "gdfdrms.pl:DEBUG:Found dird value : <$dird>\n";}
               if ($source eq "egsefm")
               {
                  #parse files like this 20071108.0x0013 where 0013 is hex apid value
                  if ($dflg == 1) {print LF "gdfdrms.pl:DEBUG::egsefm:parsing based on formated files from src eq lm\n";}
                  #get apid need to convert hex apid to decimal apid
                  $apid=  hex substr($dird,11,4);
                  if($dflg) {print LF "gdfdrms.pl:DEBUG:Found <$apid> for source <$source>\n"};
               }
               elsif ($source eq "rtmon")
               {
                  #parse files like this 20071108.0x0013 where 0013 is hex apid value
                  if ($dflg == 1) {print LF "gdfdrms.pl:DEBUG:rtmon:parsing based on formated files from src eq lm\n";}
                  #get apid need to convert hex apid to decimal apid
                  $apid=  hex substr($dird,11,4);
                  if($dflg) {print LF "gdfdrms.pl:DEBUG:Found <$apid> for source <$source>\n"};
               }
               elsif ($source eq "hsb")
               {
                  #parse files like this hsb_0445_2007_11_08_16_51_31_00.hkt
                  if ($dflg == 1) {print LF "gdfdrms.pl:DEBUG:hsb:parsing based on formated files from src eq hsb\n";}
                  #get apid which is a decimal value for hsb filenames
                  $apid=   int substr($dird,4,4);
                  if($dflg) {print LF "gdfdrms.pl:DEBUG:Found <$apid> for source <$source>\n"};
               }
               elsif ($source eq "moc")
               {
                  #parse files like this 0129_2008_260_02.hkt
                  if ($dflg == 1) {print LF "gdfdrms.pl:DEBUG:moc:parsing based on formated files from src eq moc\n";}
                  #get apid which is a decimal value for hsb filenames
                  $apid=   int substr($dird,0,4);
                  if($dflg) {print LF "gdfdrms.pl:DEBUG:Found <$apid> for source <$source>\n"};
               }
               # check if apid is in apid list
               foreach $a (@all_apids)
               {
                  if($dflg) {print LF "gdfdrms.pl:DEBUG:loop thru apids:apid ==<$a>\n";}
                  if (int $a eq $apid)
                  {
                    if($dflg == 1) {print LF "gdfdrms.pl:DEBUG:Found $apid for file in the list of apids to do\n";}
                    push (@there_dir_list, "$item/$dird");
                    last;
                  }
               }
             }
           }#foreach dird 
         }#else 
         close DFDIR;
       }#foreach item
       print LF "-->(5)List of dayfiles for dates and apids are:\n@there_dir_list\n";
     }#foreach date
  }#if date
  else
  {
    print LF "Error:get_list_from_drms() function in perl script gdfdrms.pl. Exiting\n";
    exit;
  }#else
}


######################################################################################
# subroutine decode_keywords_for_dayfiles(): get data series name by looking up apid #
######################################################################################
sub decode_keywords_for_dayfiles()
{ 
 # defined list of apids to do from apid list file
 foreach $apid  (@all_apids)
 { 
    print LF "-->(6)foreach loop- doing apid:$apid\n";    
    if ($apid eq "")
    {
      next;
    }

    #get name convention of DID dayfiles file-this contains absolute path of dayfiles already processed once.
    $didfile="";
    $didfile=sprintf("%s/DID_DECODE_DF_%s.txt", "$script_dir",int $apid);
    $didfile=~ s/\n//g;

    #open DID files for each APID and get all dayfiles DID already
    print  LF "-->(7)Opening didfile is <$didfile> \n";
    @alldid_apid="";
    if (-e "$didfile")
    {
      open(DFILE, "$didfile") || die  "ERROR in gdfdrms.pl:Can't Open-4 to read:$didfile  $!\n";
      while (<DFILE>)
      {
        $_=~ s/\n//g;
        push(@alldid_apid, $_) ;
        if($dflg) {print LF "gdfdrms.pl:DEBUG: did dayfiles:<$_> \n";}
      }
    }
    else
    {
      #create new file
      open(DFILE, ">$didfile") || die  "ERROR in gdfdrms.pl:Can't Open-5 to read:$didfile  $!\n";
      if($dflg) {print  LF "gdfdrms.pl:DEBUG: didfile is <$didfile> and was created \n";}
    }
    #close after reading or creating
    close DFILE;

      
    #open DID file again for writing dayfiles did
    open(DFILE, ">>$didfile") || die "ERROR in gdfdrms.pl:Can't Open-6 to write:$didfile  $!\n";

    #loop thru each THERE file to check if in the DID file list
    foreach $there (@there_dir_list)
    {
      
      print LF "-->(8) Checking if need to decode_dayfile the dayfile:<$there>\n";
      $hex_apid = sprintf("%04x", $apid);
      if($dflg) {print LF "gdfdrms.pl:DEBUG:decimal apid is <$apid>  hex apid is <$hex_apid>\n";}

      # get filename using regular expression
      my($directory, $filename) = $there =~ m/(.*\/)(.*)$/;

      #double check if apid is in filename
      if (( $apid eq int substr($filename,4,4) and $source eq "hsb") or ($hex_apid eq substr($filename,11,4) and $source eq "egsefm") or ( $apid eq  int substr($filename,0,4) and $source eq "moc")  or ($hex_apid eq substr($filename,11,4) and $source eq "rtmon"))
      {
        if ($dflg == 1) {print LF "gdfdrms.pl:DEBUG:there file to do is <$there> hex value of apid is <$hex_apid> decimal value of apid <$apid>\n";}
        # set found flag to no 
        $found="n";

        # loop thru DID files and check if in the THERE files list
        foreach $did (@alldid_apid)
        {
          if($dfg == 1) { print LF "gdfdrms.pl:DEBUG: did is <$did> and there <$there>\n"};
          #if did equal there skip since already done
          if ( $did eq $there)
          {
            print LF "-->(9)skip doing: there:<$there> did:<$did>\n";
            $found="y";
            last;
          }
          else
          {
            if($dflg) {print LF "gdfdrms.pl:DEBUG:not equal there vs did-continue looking for:there:<$there> did:<$did>\n"};
            next;
          }
        }
        # check if found dayfile in DID list if not-decode keywords for dayfile
        if  ( $found eq "n" )
        {
           print  LF "-->(10)Found one dayfile to do $there \n"; 

           if ($rptflg == 1)
           {
             # create series name and report name
             if ($source eq "egsefm")
             {
               $series= sprintf("%s.%s_%04d", $pjn,$dtn,hex substr($filename,11,4));
             }
             elsif ($source eq "rtmon")
             {
               $series= sprintf("%s.%s_%04d", $pjn,$dtn,hex substr($filename,11,4));
             }
             elsif($source eq "hsb")
             {
               $series= sprintf("%s.%s_%04d", $pjn,$dtn, substr($filename,4,4));
             }
             elsif($source eq "moc")
             {
               $series= sprintf("%s.%s_%04d", $pjn,$dtn, substr($filename,4,4));
             }
	     print LF  "-->(11)Series name used for reportname and creating hk series is $series \n";
	     print LF  "-->(12)Report name for report is Report-$filename-$series \n";
           }

           # Check if just want to view what would be loaded without actually loading.
           if($view_flg eq "0")
           {
             if( $rptflg == 1)
             {
	       # load this there file in DRMS
               print LF "-->(13)Running script in mode for sending dayfiles to decode_dayfile exec to be decoded and written to DRMS with report file. \n";
               $log1=`decode_dayfile -p src=$source in=$there > $script_dir/Report-$filename-$series`;

               # gzip report
               $log2=`/bin/gzip -f $script_dir/Report-$filename-$series`;

               # move report to familar place or to be deterined place if needed
               # $log3=`/bin/mv $script_dir/Report-$filename-$series.gz /home/carl/REPORTS-KEYWORDS-INGESTED/ `;

               # add filename to did file list when successfully processed dayfile already
               print DFILE "$there\n";

               # print log results to log file
               print LF "-->(14)\n--->log from running decode_dayfile: $log1\n--->log from gzip report: $log2 \n";
               #print LFILE "(15)log from moving report: $log3 \n";
             }
             else
             {
	       # load this there file in DRMS
               print LF "-->(13)Running script in mode for sending dayfiles <$there> to decode_dayfile exec to be decoded and written to DRMS without creating report file\n";
               $log1=`decode_dayfile src=$source in=$there`;

               # add filename to did file list when successfully processed dayfile already
               print DFILE "$there\n";
             }
           } #if view flag is not there
           else
           {
              print LF "-->(13)Running script in viewimg mode-keywords will not be created in DRMS \n";
           }
        } #if n-then found dayfile to process 
        elsif ($view_flg eq "0")
        {
           print LF "-->(16)Did not find dayfile to be processed \n";
        }
        else
        {
           print LF "-->(16)Running script in viewimg mode-keywords will not be created in DRMS \n";

        }
      }#if outer
    }
    if(@there_dir_list == "")
    {
      print LF "-->(8)Did not find dayfile to be processed\n";
    }
    close DFILE;
  }
  print  LF "-->(17)Finished -check above on processing done\n";
}
