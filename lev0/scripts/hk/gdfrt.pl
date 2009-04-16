#!/usr/bin/perl
##############################################################################
# Name:        gdfrt.pl - Get  dayfiles from real time files, break up into  #
#                         minute files and send to decode_dayfile for writing#
#                         keyword values to DRMS data series.                #
# Description: Used to gather Level 0 SDO formatted dayfiles from low speed  #
#              bus using today date. Then using dd unix command breaks a     #
#              minute of data and creates a minute dayfile. Then the minute  #
#              dayfile is sent to decode_dayfile executable. decode_dayfile  #
#              executable decode packets in dayfile and writes values to     #
#              DRMS data series. Then the minute dayfile is deleted.         #
#              Currently reads in all data in dayfile for first minute. Then #
#              after sleeping 60 seconds reads next new set of packets for   #
#              next minute, etc. Set the debug flag to 0,1 or 2 to control   #
#              the amount of debug to see($dflg).                            #
#              To be done in future maybe:                                   #
#              Add flexibility to write report of processing of              #
#              keywords when running decode_dayfile with -p flag.            #
# Execution:   (1)The run options are show in help listing:                  #
#                   gdfdrms  -h                                              #
#              (2)The execution options are shown using "gdfdrms  -h" or see #
#                 below the execution format and an examle.                  #
#                   gdfrt.pl apid=<decimal apid value> pkt_size=<value>      #
#                             Where value is packet size in packet + 7       #
#                   Example:  gdfrt.pl apid=129  pkt_size=128                #
#                   Example:  gdfrt.pl apid=29  pkt_size=136                 #
# Limitation:  Setup required environment variables at top of file based on  #
#              where running. decode_dayfile setup is needed since this      #
#              script executes decode_dayfile. Script works only using real  #
#              time dayfile file with LMSAL file formats for source equal to #
#              rtmon. Here are examples input dayfile file formats.          #
#              -For src=rtmon for HMI ISP dayfile: 200804019.0x01d           #
#              -For src=rtmon for AIA ISP dayfile: 200804019.0x027           #
#              -For src=rtmon for SDO ASD dayfile: 200804019.0x081           #
# Author:      Carl                                                          #
# Date:        Created April, 7, 2009                                        #
##############################################################################
# main program                                                               #
##############################################################################

#(1)common setting for all environmnents
$ENV{'SUMSERVER'}="d02.Stanford.EDU";
$hm=$ENV{'HOME'};
$ENV{'MAILTO'}="";
$exec_dir=$ENV{'DF_EXEC_PATH'}="$hm/cvs/JSOC/bin/linux_x86_64";
$script_dir=$ENV{'HK_SCRIPT_DIR'}="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'HK_DF_RT_LOGFILE'}="$hm/cvs/JSOC/proj/lev0/scripts/hk/gdfrt-log-rtmon";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir:$exec_dir";

#(2) set choices for input dayfile directory
# use for testing on Carl workspace on cl1n002
# input dayfile directory location with subdirectories with apid names(0x081,0x029,etc).
$ENV{'HK_DDF_RT_DIRECTORY'}="$hm/cvs/JSOC/proj/lev0/scripts/hk/data";
#NOTE(1):Set Production Path probably somewhere on hmisdp-mon machine and comment out tmp location used for testing!

#(3) set choices for minute file directory
# use for testing on Carl workspace on cl1n002
# minute dayfile directory location - this is temp location- since all files are deleted  
$minfiles=$ENV{'HK_DDF_MINUTE_FILE_DIRECTORY'}="$hm/cvs/JSOC/proj/lev0/scripts/hk/minute-files";
#NOTE(2):Set Production Path probably to scratch and comment out temp location!

#(4) set debug flag 
$dflg=$ENV{'DF_GDFRT_DEBUG'}=1;#use 0 to display min debug, 1 to send log to file or use 2 to send max debug log to display

#(5) check for any arguments passed in command
&check_arguments();

#(6) setup log file -use logfile setting set in hmi_dfd.pl,aia_dfd.pl, sdo_dfd.pl or use setting below
$logfile=$ENV{'HK_DF_RT_LOGFILE'};
if($logfile eq "")
{
   $logfile="$script_dir/log-gdfrt";#send log here when run at command line
}

#(7) open log file
open(LF,">>$logfile") || die "ERROR in gdfrt.pl:(0):Can't Open log file <$logfile>: $!\n; exit;";
print LF `date`;

#(8) get list of apids to do 
@avlist=&get_apid_to_do();
$av= pop(@avlist);
if($dflg == 2) {print "DEBUG:MESSAGE:gdfrt: apid value is $av\n";}

#(9)loop forever thru today and next days 
while (1)
{
  #flush print buffer
  $|=1;

  #(10)get list of dates to do 
  &get_date_to_do();
 
  ### LOG ###
  if($dflg == 2) {print "DEBUG:MESSAGE:gdfrt:calling get-rt-file function with apid=$av date-range1=$date_r1 date-range2=$date_r2\n";}

  #(11)go through day files in directory looking for all files for given apid and date range
  @list_rtfiles=&get_rt_files($date_r1, $date_r2, $av);

  ### LOG ###
  if($dflg == 0) {print "-->Listing of rt file to do: <@list_rtfiles> based on date range <$date_r1> to <$date_r2>\n";}
  if($dflg == 1) {print LF "-->(3)Listing of rt file to do: <@list_rtfiles> based on date range <$date_r1> to <$date_r2>\n";}
  if($dflg == 2) {print "DEBUG:MESSAGE:gdfrt:Listing of rt file to do: @list_rtfiles\n";}

  #(12) get dayfile to read packets and write to minute files
  $f = pop(@list_rtfiles);
  if ($f eq "")
  {
    ### LOG ###
    if($dflg == 2) {print "DEBUG:MESSAGE:gdfrt:got no files to process for apid:$av . Sleep 1 minute and check again\n";}

    # sleep 1 minute and go to top of while loop and check again
    sleep 60;
    next;
  }

  #(14) get today date to handle reading next day's dayfile
  $c_date=&get_todays_date();

  #(15)get packet size which is done by just getting argument passed for pkt_size=<value>
  $packet_size= &get_packet_size();
  
  #(16)read as many packets that are in dayfile and then read next new set of packets every 1 minute
  &read_n_packets($f,$c_date, $packet_size);

  #(17) continue while loop forever reading one dayfile at a time for each of the next days to come.
}
 
#(18)close logfile
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
  $pkt_size_flg="0";
  
  if ($#ARGV < 0)
  {
    $help_flg="1";
  }
  if ($#ARGV >= 0 && $#ARGV <= 1)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-H")
    {
      $help_flg = "1";
    }
    elsif (substr($ARGV[0],0,5) eq "apid=" )
    {
      #use apid following -a to create map files 
      #push decimal character value in this format dddd. Example: -a 0001
      $apid_list_flg = "2";

      if (substr($ARGV[1],0,9) eq "pkt_size=")
      {
        $pkt_size_flg=1;#set to trigger using arg 1 later
        if($dflg == 2) { print "DEBUG:MESSAGE:check_arguments:pkt size is $pkt_size\n";}
      }
      else
      {
        print  LF "gdfrt.pl:ERROR: Did not use pkt_size argument. Exiting script\n";
        print  "gdfrtf.pl:1:ERROR: Did not use pkt_size argument. Exiting script\n";
        exit;
      }
    }
  }
  else
  {
    print  LF "gdfrt.pl:ERROR: Entered too many arguments. Exiting script\n";
    print  "\ngdfrtf.pl:3:ERROR: Entered too many arguments. Exiting script\n";
    print "NOTE:Usage:   gdfrt.pl  apid=< APID value in decimal number format>  pkt_size=<packet size in packet + 7>\n";
    print "NOTE:Example: gdfrt.pl apid=129  pkt_size=128\n\n";
    exit;
  }

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
  print "(1a)Get dayfile real time data using todays date and apid value to determine which dayfile to process.\n\n";
  print "         gdfrt.pl  apid=< APID value in decimal number format>  pkt_size=<packet size in packet + 7>\n\n";
  print "         Example: gdfrt.pl apid=129  pkt_size=128\n";
  print "         Example: gdfrt.pl apid=29   pkt_size=136 \n\n";
  print "(1b)Get Help Information:\n\n";
  print "         gdfrt.pl -h  or  gdfrt.pl -help\n\n";
  print "(2) Requires that hk level0  data series(i.e.sdo.lev0_asd_0003,hmi.lev0_isp_0021,etc.) be already created.\n";
  print "(3) Required that the environment variable are setup correctly at top of script(HK_DDF_RT_DIRECTORY, etc).\n";
  print "(4)****Limitation****:\n";
  print "(4a)Works only on rtmon dayfile formats only(20090407.0x081)\n";
  print "(4b)Enter arguments in specified order when running script is required.\n";
  print "(4c)It is required that a correct setup for decode_dayfile execution is done since this scripts uses decode_dayfile.\n";
  print "(4d)It is required that a correct value be use for packet size argument. Get packet size in data file or hk config\n";
  print "    file then add 7 to get value to use(for sdo asd packet, 121 + 7 = 128).\n\n\n";
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
    open(FILE_APID_LIST, "$fn") || die "ERROR in gdfdrms.pl:(1):Can't Open: $fn file: $!. Need to create apidlist file.\n";
    while (<FILE_APID_LIST>)
    {
      push(@all_apids, int $_) ;
    }
    close FILE_APID_LIST ;
    print  LF "-->(1)Doing decimal formatted apid: @all_apids\n";
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
   if($dflg == 2) {print "DEBUG:MESSAGE:all apids doing = @all_apids\n";}
   return(@all_apids);
}



##########################################################################
# subroutine get_todays_date()                                           #
##########################################################################
sub get_todays_date()
{
  ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  $year = 1900 + $yearOffset;
  $month= $monthOffset + 1;
  if($dflg == 2) {print  "DEBUG:MESSAGE:get_todays_date: year this $year month is $month day is $dayOfMonth\n";}
  #create todays date format and push on list of dates to do
  $new_date= sprintf("%4d%02.2d%02.2d", $year,$month,$dayOfMonth);
  if($dflg == 2) {print "DEBUG:MESSAGE:get_todays_date: new date is ::: $new_date :::\n";}
  return ( $new_date);
}


#############################################################################
# subroutine get_date_to_do: gets date range to do                          #
#############################################################################
sub get_date_to_do()
{
  # initialize all_dates to null
  @all_dates=();
  # get start and end range of dates
  if ($date_range_flg eq "1")
  {
    $date_r1= substr($ARGV[1],6,8);
    $date_r2= substr($ARGV[2],4,8);
    #create list of dates based on range
    &get_date_list_to_do();
  }
  else
  {
     #do case with no entry for start and end time
     #for this case can make do today
     ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
     $year = 1900 + $yearOffset;
     $month= $monthOffset + 1;
     if($dflg == 0) {print  "-->Today's Date is $month-$dayOfMonth-$year\n";}
     if($dflg == 1) {print LF "-->(2)Today's Date is $month-$dayOfMonth-$year\n";}
     if($dflg == 2) {print "DEBUG:MESSAGE:gdfrt: year this $year month is $month day is $dayOfMonth\n";}
     #create todays date format and push on list of dates to do
     $new_date= sprintf("%4d.%02.2d.%02.2d", $year,$month,$dayOfMonth);
     push(@all_dates, $new_date);
     #set values to display in log
     $date_r2=$date_r1=sprintf("%4d%02.2d%02.2d", $year,$month,$dayOfMonth);
  }
  if($dflg == 2) {print "DEBUG:MESSAGE:get_date_to_do: Date range 1 is $date_r1\n";}
  if($dflg == 2) {print "DEBUG:MESSAGE:get_date_to_do: Date range 2 is $date_r2\n";}
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
        if($dflg == 1) {print LF "-->Pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        if($dflg == 2) {print  "DEBUG:MESSAGE:get_date_list_to_do:Pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        $new_date= sprintf("%4d.%02.2d.%02.2d", $loop_yr,$loop_mo,$loop_day);
        push(@all_dates, $new_date);
        $loop_day=1;
        $loop_mo++;
      }
      elsif ($loop_mo == 12 && $loop_day == 31)
      {
         #if end of year
        if($dflg == 1) {print LF "-->Pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        if($dflg == 2) {print  "DEBUG:MESSAGE:get_date_list_to_do:Pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        $new_date= sprintf("%4d.%02.2d.%02.2d", $loop_yr,$loop_mo,$loop_day);
        push(@all_dates, $new_date);
        $loop_day=1; 
        $loop_mo=1;
        $loop_yr++;
      }
      else
      {
        #if not end of year or month then just increment day
        if($dflgi == 1) {print LF "gdfrt.pl:Pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        if($dflgi == 2) {print  "DEBUG:MESSAGE:get_date_list_to_do:Pushing this date: $loop_yr.$loop_mo.$loop_day\n"};
        $new_date= sprintf("%4d.%02.2d.%02.2d", $loop_yr,$loop_mo,$loop_day);
        push(@all_dates, $new_date);
        $loop_day++;
      }
      #set up next loop date number
      $loop_yr_mo_day=sprintf("%4d%02.2d%02.2d", $loop_yr,$loop_mo,$loop_day);
    }
    if($dflg == 1) {print LF "List of dates to process:\n @all_dates\n"};
    if($dflg == 2) {print "DEBUG:MESSAGE:get_date_list_to_do:List of dates to process:\n @all_dates\n"};


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



######################################################################################
# subroutine get_rt_files(): get real time dayfiles and send back file               #
######################################################################################
sub get_rt_files($$$)
{
  # setup my variables
  my($dr1,$dr2, $dn_fn);
  $dr1=$_[0];#date range 1
  $dr2=$_[1];#date range 2
  $apid_value=$_[2];#apid value

  #initialize list of file in date range to null
  @file_list_in_date_range=();

  if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files: Passed argument values:$dr1:$dr2:$apid_value:\n";}

  # get directory name to find apid directories 
  $dn = sprintf("%s",  $ENV{'HK_DDF_RT_DIRECTORY'});
  opendir(DIR_RT, $dn) || die "ERROR in gdfrt.pl:(3):Can't open rt directory <$dn>. Check setting of HK_DDF_RT_DIRECTORY environment variable:$!\n"; #open subdirectory
  if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files: Directory to find dayfile directories:$dn:\n";}

  # read in all directories  and put in list 
  @dir_list = readdir(DIR_RT); #get a list of directory contents

  #check which apid directories to process
  foreach $dir  (@dir_list)
  {
    if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files: Looking for directory for apid: $apid_value: Got directory:$dir:\n";}
    if( hex substr($dir,2,5) ne  $apid_value)
    {
       next;
    }
    else
    {
      if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files: Found apid directory:$dir:\n";}
      push(@apiddir_to_do_list, $dir);
      last;
    }
  }

  # check each apid directory for dayfiles
  foreach $apid_dir  (@apiddir_to_do_list)
  {
    #found directory for this APID
    #check day files to do in each directory
    $dn_fn = sprintf("%s/%s",  $ENV{'HK_DDF_RT_DIRECTORY'}, $apid_dir);
    if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files:Checking this directory to find this apid's dayfiles:$dn_fn:\n";}

    #open directory looking for this apids dayfiles
    opendir(DIR_RT_APID, $dn_fn) || die "ERROR in gdfrt.pl:(4):Can't open apid directory <$dn_fn>:$!\n"; #open subdirectory

   #read all files for apid directory
   @apid_filelist = readdir(DIR_RT_APID); #get a list of directory contents

   #loop thru each file looking for file with today's date
   foreach $apid_file  (@apid_filelist)
   {
     # Filter out 20090327.0x0081x 20090320.0x0081x and . and .. files
     if (substr($apid_file,15,1) eq "x" || $apid_file eq ".." | $apid_file eq ".")
     {
        next;
     }#outer if
     else
     {
         # check first 8 characters of date (20090320.0x0081) is within range
         if (int substr($apid_file,0,8) >= int $dr1 && int substr($apid_file,0,8) <= int $dr2 )
         {
            if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files:Found dayfile with good date range and apid value :::$apid_file\n";}
            $fullpathname = sprintf("%s/%s",$dn_fn, $apid_file);
            push( @file_list_in_date_range, $fullpathname);
            last;
         }#inner if
         else
         {
           ;#skip not in date range
           if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files:skip no files found in date range\n";}
         }#inner else
     }#outer else
   }# foreach file in apid directory
  }
   if($dflg == 2) {print "DEBUG:MESSAGE:get_rt_files:Returned this file for processing :@file_list_in_date_range\n";}
  return (@file_list_in_date_range )
}


#####################################################################################################
# subroutine read_n_packets(): read n number of packet per minute and send to decode dayfile execute#
#####################################################################################################
sub read_n_packets($$$)
{
   # set my variables
   my($file,$apid, $f, $d, $source,$curr_date, $x,$y,$skip,$psize);
   # passed arguement
   $file=$_[0];
   $curr_date=$_[1];
   $psize=$_[2];#packet_size

   ### LOG ###
   if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packets: Called read_n_packets with arguments file:$file and todays date:$curr_date\n";}

   # set up delay to exiting loop when today's date changes 
   $wait_clock_trigger=10;
   $wait_clock=0;

   #get directory and filename and apid value
   ($d, $f) = $file=~ m/(.*\/)(.*)$/;

   #get apid from filename
   $apid= hex substr($f,11,4);

   ### LOG ###
   if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet:apid used for display is $apid\n";}
   if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet:packet size used for dd command is $psize\n";}

   # open dayfile to process to get current size of file
   open(DAYFILE, "$d/$f") || die "ERROR in gdfrt.pl:(5):Can't Open dayfile to read:$d/$f $!\n";

   #get file size
   $size = (stat(DAYFILE))[7];
   if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet: File:$d/$f: has a file size :$size:\n";}

   #pkt count
   $current_pkts_in_file= $size/$psize;
   $pkt_count= $current_pkts_in_file;

   ### LOG ###
   if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet:packet count to be used in dd command is $pkt_count\n";}
   if($dflg == 0) {print LF "-->(5)Packet count to be used in dd command is $pkt_count\n";}

   #skip initialize to 0
   $skip=0;

   # set up parameter to do if compare in while loop
   $max= $current_pkts_in_file;
   $x=0;
   #use y to increment every minute if execute dd command
   $y=0; 

   while ( 1 )
   {
      $|=1;#added to flush buffers of print.
      if ( $x < $max)
      {
        $wait_clock=0;

        # set pkt count base on how many pkts are left to do in file.
        $pkt_count = $current_pkts_in_file;
  
        ### LOG ###
        if($dflg == 2) {printf( "DEBUG:MESSAGE:read_n_packets: Executing dd command at %s\n", get_current_time());}
        if($dflg == 1) {printf( LF "-->(4)Executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s\n",get_current_time());}
        if($dflg == 0) {printf( "-->Executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s\n",get_current_time());}

        # Execute dd command to create minue files.
        $log=`dd if=$d/$f of=$minfiles/$f-$y-minute bs=$psize count=$pkt_count skip=$skip`;

        ### LOG ###
        # check for errors when running dd command
        if ($? ne 0)
        {
          #print " $? \n";
          print "ERROR: Got error executing dd command. Check packet size passed as argument is correct. Packet size used is $psize. ";
          print "When packet size is divided by file size to process. Should get even packet count number. Getting packet count to do as $pkt_count . Exiting script.\n";

          print LF "ERROR: Got error executing dd command. Check packet size passed as argument is correct. Packet size used is $psize. ";
          print LF "When packet size is divided by file size to process. Should get even packet count number. Getting packet count to do as $pkt_count . Exiting script.\n";
          exit;
        }
        else
        {
          if($dflg == 2) {printf("DEBUG:MESSAGE:read_n_packets: Complete executing dd command at %s\n",get_current_time());}
          if($dflg == 1) {printf( LF "-->(5)Completed executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s\n",get_current_time());}
          if($dflg == 0) {printf("-->Completed executing dd command writing <$pkt_count> packets to <$f-$y-minute> file at %s\n",get_current_time());}
        }


        ### LOG ###
        if($dflg == 0) {printf("-->Executing decode_dayfile command writing <$pkt_count> packets to drms series at %s\n",get_current_time());}
        if($dflg == 1) {printf(LF "-->(6)Executing decode_dayfile command writing <$pkt_count> packets to drms series at %s\n",get_current_time());}
        if($dflg == 2) {printf("DEBUG:MESSAGE:read_n_packets:Executing decode_dayfile command writing <$pkt_count> packets to drms series at %s\n",get_current_time());}

         
        # Execute decode dayfile on minute file
        $source="rtmon"; #hardcoded since this script only used for real time dayfiles
        # decode and write to drms and create minute report in minfiles directory- COMMENT OUT and used instead command below
        #$log=`$exec_dir/decode_dayfile -p src=$source in=$minfiles/$f-$y-minute > $minfiles/Report-$f-$y-minute`;
        # decode and write to drms only - create NO report!
        $log=`$exec_dir/decode_dayfile src=$source in=$minfiles/$f-$y-minute   2>&1`;

        ### LOG ###
        # check status returned when executing decode_dayfile
        $log =~ m/(ERROR)/ ; #regular exp - look for field
        if($log eq "" )
        {
          if($dflg == 0) {printf("-->Completed executing decode_dayfile with status PASSED at %s\n",get_current_time());}
          if($dflg == 1) {printf(LF "-->(7)Completed executing decode_dayfile with status PASSED at %s\n",get_current_time());}
          if($dflg == 2) {printf("DEBUG:MESSAGE:read_n_packets:Completed executing decode_dayfile with status PASSED at %s\n",get_current_time());}
        }
        else
        {
          if($dflg == 0) {print  "ERROR:Status failed when executing decode_dayfile. Exiting script gdfrt.pl.\nLOG returned:$log:\n";}
          if($dflg == 1) {print LF "ERROR:Status failed when executing decode_dayfile. Exiting script gdfrt.pl.\nLOG returned:$log:\n";}
          if($dflg == 2) {print  "ERROR:read_n_packets:Status failed when executing decode_dayfile.Exiting script gdfrt.pl.\nLOG returned:$log:\n";}
          exit;
        }

        ### LOG ###
        if($dflg == 0) {printf("-->Executing  delete of  minute file <$f-$y-minute> at %s\n",get_current_time());}
        if($dflg == 1) {printf(LF  "-->(8)Executing  delete of minute file <$f-$y-minute> at %s\n",get_current_time());}
        if($dflg == 2) {printf("DEBUG:MESSAGE:read_n_packets:Executing  delete  of minute file <$f-$y-minute> at %s\n",get_current_time());}

        # delete minute files that was already loaded into data series
        unlink("$minfiles/$f-$y-minute");

        ### LOG ###
        if($dflg == 0) {printf("-->Completed executing  delete of  minute file <$f-$y-minute> at %s\n",get_current_time());}
        if($dflg == 1) {printf(LF  "-->(9)Completed executing  delete of minute file <$f-$y-minute> at %s\n",get_current_time());}
        if($dflg == 2) {printf("DEBUG:MESSAGE:read_n_packets:Completed executing  delete  of minute file <$f-$y-minute> at %s\n",get_current_time());}
        if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packets:sleeping 60\n\n";}

        # sleep 1 minute and then repeat above
        sleep 60;

        #Reopen and get current file size 
        close (DAYFILE);
        # reopen file and check size
        open(DAYFILE, "$d/$f") || die "ERROR in gdfdrms.pl:(6):Can't Open day file to read:$d/$f  $!\n";
        $size = (stat(DAYFILE))[7];

        ### LOG ###
        if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet: File:$d/$f: has a file size :$size:\n";}

        # get total number of packets in file
        $total_pkts_in_file= $size/$psize;

        ### LOG ###
        if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packets: After reopen, total packets in file: $total_pkts_in_file\n";}

        # Number of packets created minute files for so far is current_pkts_in_file amount just did and x amount did before
        $x=$current_pkts_in_file + $x ; #number of packets did incrementally so far
 
        # Reset max value in loop to total pkts in file
        $max=$total_pkts_in_file;

        # Update skip value to amount skipped last time with the amount did just now
        $skip= $pkt_count + $skip ; #watch pkt-count value and current-pkts-in-file are same most the time but not all the time.

        ### LOG ###
        if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet: during next dd command will skip  $skip\n";}

        # calculate number of packet left to do now
        $current_pkts_in_file = $total_pkts_in_file - $x;

        # increment y
        $y++;
      
      }# if (x < max) 
      else
      {
        # get today date and time
        $td_date=get_todays_date();

        # increment delay clock
        $wait_clock++;

        ### LOG ###
        if($dflg == 2) {print "DEBUG:MESSAGE:If today's date:$td_date > $curr_date than previous day, then break from loop.\n";}

        # check if waited enough time and if its a new day
        if(($wait_clock >= $wait_clock_trigger) && ($td_date > $curr_date))
        {
          ### LOG ###
          if($dflg == 0) {print  "-->Next day detected. Today date is now <$td_date> and previous day was <$curr_date>. Try getting new dayfile to process.\n";}
          if($dflg == 1) {print LF "-->(11)Next day detected. Today date is now <$td_date> and previous day was <$curr_date>. Try getting new dayfile to process.\n";}
          if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packets: Next day detected. Break from loop and reread new dayfile for next day.\n";}

          # reset current date
          $curr_date=$td_date;

          # break from this loop to check for dayfiles for next day
          last;
        }

        #### LOG ###
        if($dflg == 0) {printf("-->No new packet in dayfile sleeping 1 minute and will check again for new packets in dayfile at %s\n",get_current_time());}
        if($dflg == 1) {printf(LF  "-->(10)No new packet in dayfile sleeping 1 minute and will check again for new packets in dayfile at %s\n",get_current_time());}
        if($dflg == 2) {printf("DEBUG:MESSAGE:read_n_packets:No new packets in dayfile, skip writing minute file, sleep 60 seconds and wait for more at %s\n",get_current_time());}

        # sleep 1 minute and then repeat above
        sleep 60;

        #Reopen file and get current file size
        close (DAYFILE);
        open(DAYFILE, "$d/$f") || die "ERROR in gdfdrms.pl:(7):Can't Open dayfile to read file:$d/$f  $!\n";
        $size = (stat(DAYFILE))[7];

        ### LOG ###
        if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet: File:$d/$f: has a file size :$size:\n";}

        # get total number of packets in file
        $total_pkts_in_file= $size/$psize;

        ### LOG ###
        if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packets: After reopen, total packets in file: $total_pkts_in_file\n";}

        # Number of packets created minute files for so far is current_pkts_in_file amount just did and x amount did before
        $x=$current_pkts_in_file + $x ; #number of packets did incrementally

        # Reset max value in loop to total pkts in file
        $max=$total_pkts_in_file;

        ### LOG ###
        if($dflg == 2) {print "DEBUG:MESSAGE:read_n_packet: during next dd command will skip  $skip\n";}

        # calculate number of packet left to do now
        $current_pkts_in_file = $total_pkts_in_file - $x;
         
      }#else if x = max
   } # end while
   close (DAYFILE);
}

#####################################################################################################
# subroutine get_packet_size(): gets packet size argument passed                                    #
#####################################################################################################
sub get_packet_size($)
{
  my($p_size);

  if (  $pkt_size_flg eq 1 )
  {
    $p_size= substr($ARGV[1],9);#get pkt size
  }
  elsif (  $pkt_size_flg eq 2 )
  {
    $p_size= substr($ARGV[3],9);#get pkt size
  }
  else
  {
    print "ERROR: Packet size flag not set. Missing packet size value\n";
    print LF "ERROR: Packet size flag not set. Missing packet size value\n";
  }
  return ($p_size);

}

##########################################################################
# subroutine get_current_time()                                          #
##########################################################################
sub get_current_time()
{
  my($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings, $year, $month, $new_date);
  ($second, $minute, $hour, $dayOfMonth, $monthOffset, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  $year = 1900 + $yearOffset;
  $month= $monthOffset + 1;
  #create todays date and time format 
  $new_date= sprintf("%4d.%02.2d.%02.2d_%02.2d:%02.2d:%02.2d",$year,$month,$dayOfMonth,$hour,$minute,$second);
  return ( $new_date);
}
