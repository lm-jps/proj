#!/usr/bin/perl
##############################################################################
# Name:        set_mp.pl - create image location file and executes load_mp to#
#                          update or add value in master pointing series     #
# Description: Prompts users to makes image_location file and then gets data #
#              from iss and leg series and load values in master pointing    #
#              series by executing load_mp                                   #
#              table.                                                        #
# Execution:   (1)To run :   set_mp.pl                                       #
#              (2)For help:  set_mp.pl and review guide                      #
# Process:     See guide.                                                    #
# Environment Variable: Set master pointing series in script. Currently use  #
#              production name.                                              #
# Requirement: View requirements in guide information.                       #
# Author:      Carl                                                          #
# Date:        Created March, 23, 2010                                       #
##############################################################################
# main function                                                              #    
##############################################################################
# set environment variable
$hm=$ENV{'HOME'};
###$exec_dir="$hm/cvs/JSOC/proj/lev1/apps";
$exec_dir="$hm/cvs/JSOC/_linux_x86_64/proj/lev1/apps";

#set place to retrieve default data from
$mp_sn=$ENV{'HK_SET_MP_MASTER_POINTING_SERIES'}="sdo.master_pointing";
###JIM'S VERSION::$mp_sn=$ENV{'HK_SET_MP_MASTER_POINTING_SERIES'}="su_carl.master_pointing";
##Test VERSION::$mp_sn=$ENV{'HK_SET_MP_MASTER_POINTING_SERIES'}="su_carl.test99_master_pointing";

# get keywords to do in master pointing series
$metadata_kw=`show_info -j ds=$mp_sn | egrep "Keyword:T_START|Keyword:T_HKVALS|T_STOP"`;
$hmi_manual_kw=`show_info -j ds=$mp_sn        | egrep "Keyword:H_"`;
$hmi_cam1_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:H_CAM1"`;
$hmi_cam2_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:H_CAM2"`;
$aia_w094_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:A_094"`;
$aia_w131_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:A_131"`;
$aia_w171_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:A_171"`;
$aia_w193_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:A_193"`;
$aia_w211_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:A_211"`;
$aia_w304_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:A_304"`;
$aia_w335_manual_kw=`show_info  -j ds=$mp_sn  | egrep "Keyword:A_335"`;
$aia_w1600_manual_kw=`show_info -j ds=$mp_sn  | egrep "Keyword:A_1600"`;
$aia_w1700_manual_kw=`show_info -j ds=$mp_sn  | egrep "Keyword:A_1700"`;
$aia_w4500_manual_kw=`show_info -j ds=$mp_sn  | egrep "Keyword:A_4500"`;

system("clear");
&print_guide();

###########################################
# step 1 Set T_START,T_STOP,T_HKVALS time #
###########################################
#show user which values to enter next
print "Step(1)Select y(yes) to set new value for T_START, T_STOP, T_HKVALS keywords or select n(no) or return key to use current values\n";

# show current keyword values in master pointing series
@ret_kws = get_keyword_names($metadata_kw,"Current T_START Values:");

#show current values
print ":Current values:\n@ret_kws\n";

# get list of keyword-name and keyword-values pairs
%kw_val_ht = set_values(@ret_kws);

#adjust display
system("sleep 1");
system("clear");
&print_guide();

############################
# step 2 hmi image location#
############################
push(@chunks,$hmi_cam1_manual_kw);
push(@chunks,$hmi_cam2_manual_kw);
$camera=1;
foreach $kw_chunk (@chunks)
{
  #show user which values to enter next
  print "Step(2)Select y(yes) to set new values for  HMI Camera $camera Image Location keywords or select n(no) or return key to use current value\n";

  # show current keyword values in master pointing series
  @ret_kws="";
  @ret_kws = get_keyword_names($kw_chunk,"Current HMI Image Camera $camera Location Keywords Values:");

  #show current values
  print ":Current values:\n@ret_kws\n";

  # get list of keyword-name and keyword-values pairs
  %kw_val_ht = set_values(@ret_kws);
  $camera++;

  #adjust display
  system("sleep 1");
  system("clear");
  &print_guide();
}

############################
# step 3 aia image location#
############################
@chunks= ("$aia_w094_manual_kw","$aia_w131_manual_kw","$aia_w171_manual_kw","$aia_w193_manual_kw","$aia_w211_manual_kw","$aia_w304_manual_kw","$aia_w335_manual_kw","$aia_w1600_manual_kw","$aia_w1700_manual_kw","$aia_w4500_manual_kw");
@wave_label= ("094","131","171","193","211","304","335","1600","1700","4500");
$i=0;
foreach $kw_chunk (@chunks)
{

  #show user which values to enter next
  print "Step(3)Select y(yes) to set new values for AIA Wavelength $wave_label[$i] Image Location keywords or select n(no) or return key to use current value\n";

  # show current keyword values in master pointing series
  @ret_kws="";
  @ret_kws = get_keyword_names($kw_chunk,"Current AIA Wavelength $wave_label[$i] Image Location Keywords Values:");

  #show current values
  print ":Current values:\n@ret_kws\n";

  # get list of keyword-name and keyword-values pairs
  %kw_val_ht = set_values(@ret_kws);
  $i++;

  #adjust display
  system("sleep 1");
  system("clear");
  &print_guide();
}

#########################################
# step 4 save values to image loc file  #
#########################################
#print everything to image location file
print_image_loc_file( %kw_val_ht);

#########################################
# step 5 load mp data in series         #
#########################################

print "Step(5)Run load_mp C executable on latest created file.\n\n";

# manual entry of time keywords values if user selects yes
$answer = &promptUser("Do you want to load this image loc file using load_mp executable(y,n)?","n");

# check if got yes then prompt for new value
if ($answer eq "y" )
{
   print "-->starting load using: load_mp ilf=$image_loc_file\n";
   print "-->values be updated in the  master pointing series set in load_mp executable.\n";
   print "-->master pointing series name is set by load_mp executable using SOURCE_ENV_LOAD_MP file.\n";
   $log=`ls $exec_dir/load_mp`;
   print "-->using executable at:$log";
   $log=`$exec_dir/load_mp ilf=$image_loc_file`;
   print "-->results of log:$log\n";

}
else
{
   print "-->if want to run later use this command: $exec_dir/load_mp ilf=$image_loc_file\n";
   print "-->currrently- the values will be updated in the master pointing series set in load_mp executable which reads series from SOURCE_ENV_LOAD_MP file.\n";
}
print "-->Completed running set_mp.pl script. Exiting.\n";

############################
# promptUser               #
############################
sub promptUser 
{

   #-------------------------------------------------------------------#
   #  two possible input arguments - $promptString, and $defaultValue  #
   #  make the input arguments local variables.                        #
   #-------------------------------------------------------------------#

   local($promptString,$defaultValue) = @_;

   #-------------------------------------------------------------------#
   #  if there is a default value, use the first print statement; if   #
   #  no default is provided, print the second string.                 #
   #-------------------------------------------------------------------#

   if ($defaultValue) {
      print $promptString, "[", $defaultValue, "]: ";
   } else {
      print $promptString, ": ";
   }

   $| = 1;               # force a flush after our print
   $_ = <STDIN>;         # get the input from STDIN (presumably the keyboard)


   #------------------------------------------------------------------#
   # remove the newline character from the end of the input the user  #
   # gave us.                                                         #
   #------------------------------------------------------------------#

   chomp;

   #-----------------------------------------------------------------#
   #  if we had a $default value, and the user gave us input, then   #
   #  return the input; if we had a default, and they gave us no     #
   #  no input, return the $defaultValue.                            #
   #                                                                 # 
   #  if we did not have a default value, then just return whatever  #
   #  the user gave us.  if they just hit the <enter> key,           #
   #  the calling routine will have to deal with that.               #
   #-----------------------------------------------------------------#

   if ("$defaultValue") {
      return $_ ? $_ : $defaultValue;    # return $_ if it has a value
   } else {
      return $_;
   }
}

###################################
# get keyword name and value      #
###################################
sub get_keyword_names($,$)
{
  #set up local variables
  my($kwds0, $hmi_aia_str,@kw_list, @kw_names,$qsn);
  @kwvalue_list="";

  # get argument passed
  $kwds=$_[0];
  $hmi_aia_str=$_[1];

  #split into lines
  @s_kwds=split ("\n",$kwds);

  #loop thru list of keywords and print values
  $value="0.000" ; #until get real values
  $j=0;
  foreach $i ( @s_kwds)
  {
    #take out keyword field
    $s_kwds[$j++]=~ /(\w+):(\w+)/; #regular expression-fish out field
    ##put together using mp_sn: $valuepair=`show_info ds='su_carl.test99_master_pointing[\$] ' key=$2 n=1`;
    $qsn=sprintf("%s[\\\$]",$mp_sn);
    $valuepair=`show_info ds=$qsn key=$2 n=1`;
    $valuepair =~ s/\n/ = /; #replace cr with equal -regular exp
    push(@kwvalue_list,  $valuepair);
  }
  return @kwvalue_list;
}

#######################
#   set_values        #
#######################
sub set_values(@)
{
  # get list of keyword-name and keyword-values pairs
  my(@ret_kws);
  @ret_kws= @_;
  foreach $pair (@ret_kws)
  {
    $answer="";
    if ($pair eq "")
    {
      next;#do nothin'
    }

    # split two values to get name and value
    @s_pair= split(" ",  $pair);
    
    # manual entry of time keywords values if user selects yes
    $answer = &promptUser("Do you want to enter a new value for $s_pair[0](y,n)?","n");

    if ($answer eq "y" || $answer eq "n")
    {
      ;#skip-print "y-n ans:<$answer>\n";
    }
    else
    {
      print "Retry bad value entered:<$answer>\n";
      while (1)
      {
        $answer = &promptUser("Do you want to enter a new value for $s_pair[0](y,n)?","n");
        if($answer eq "y" || $answer eq "n")
        {
          last;
        }
        else
        {
          print "Retry bad value entered::<$answer>\n";

        }
      }
    }

    # check if got yes then prompt for new value
    if ($answer eq "y" )
    {
      # prompt user to set values 
      $entry_str="Enter a new value for $s_pair[0]?";
      $default_str="";
      # call have user enter values and return in updated retvalues list.
      $ret="n";
      while ($ret eq "n")
      {
        #do set_modify_kw_value in-line
        $kw_value = &promptUser($entry_str,$default_str);
        $len=length($kw_value);
        if($len == 0)
        {
          $kw_value="NO_VALUE_ENTERED";
          print "No Value Entered <> - Retry Entering.\n";
          $ret="n";
          next;
        }
        $ret= &promptUser("Are these values set correctly(y,n)[default=y]","y");

        # check if user enter bad y-n-cr answer
        if ($ret eq "y" || $ret eq "n")
        {
          ;#skip-print "y-n ans<$ret>\n";
        }
        else
        {
          print "Retry bad value entered:<$ret>\n";
          while (1)
          {
            $ret= &promptUser("Are these values set correctly(y,n)[default=y]","y");
            if($ret eq "y" || $ret eq "n")
            {
              last;
            }
            else
            {
              print "Retry bad value entered:<$ret>\n";
            }
          }
        }

        # add values to hash table if got y
        if ($ret eq "n")
        {
          next;
        }
        elsif($ret eq "y")
        {
          push( @kwvalue, $kw_value);
          push( @kwname, $s_pair[0]);
          $kw_value_ht{$s_pair[0]} =  $kw_value;
        }
        else
        {
          $kw_value="NO_VALUE_ENTERED";
          print "bad select <$ret> - Retry Entering.\n";
          $ret="n";
          next;
        }
      }#while
    }#endof big if
    elsif($answer eq "n") 
    { 
      print "Using this value: $s_pair[0] = $s_pair[2]\n" ;
      # if selected no then push default value on list and use this for value
      push(@kwname, $s_pair[0]);
      push( @kwvalue, $s_pair[2]);
      $kw_value_ht{$s_pair[0]} =  $s_pair[2];
    }
    else
    {
      print "bad select <$ret> - Retry Entering.\n";
      print "Using this value: $s_pair[0] = $s_pair[2]\n" ;
      # if selected no then push default value on list and use this for value
      push(@kwname, $s_pair[0]);
      push( @kwvalue, $s_pair[2]);
      $kw_value_ht{$s_pair[0]} =  $s_pair[2];
    }
       
  }#for each pair if keyword-Value
  return(%kw_value_ht);



}

#################################
# print_image_loc_file          #
#################################
sub print_image_loc_file(%)
{
  my($ht_values);
  %ht_values= %_;

  # print keyword and value to file
  #display and print results to image_location.txt file
  # open image location file 
  $image_loc_file="./sdo_image_location.txt";
  open(OUTFILE, ">$image_loc_file") || die "(6)Can't Open $ifn file: $!\n";
  printf( OUTFILE "# Filename : sdo_image_location.txt\n");
  $d=`date`;
  $d=~s/\n//;
  printf( OUTFILE "# Date : %s\n", $d);

  # print Times
  while(($key,$value) = each (%kw_val_ht))
  {
    if( $key eq "T_START" || $key eq "T_STOP" || $key eq "T_HKVALS")
    {
      printf( OUTFILE "%-3.3s %-20.20s %-25.25s %-7.7s\n","KWD",$key,$kw_val_ht{$key},"time");
    }
  }

  # print version
  printf( OUTFILE "%-3.3s %-20.20s %-25.25s %-7.7s\n","KWD","VERSION","0","int");

  # print hmi values
  while(($key,$value) = each (%kw_val_ht))
  {
    if( substr($key,0,6) eq "H_CAM1")
    {
      printf( OUTFILE "%-3.3s %-20.20s %-25.25s %-7.7s\n","KWD",$key,$kw_val_ht{$key},"float");
    }
  }
  while(($key,$value) = each (%kw_val_ht))
  {
    if( substr($key,0,6) eq "H_CAM2")
    {
      printf( OUTFILE "%-3.3s %-20.20s %-25.25s %-7.7s\n","KWD",$key,$kw_val_ht{$key},"float");
    }
  }

  # create wavelist order in which they will appear in image_location_file.txt
  @wavelist= ("A_094","A_131","A_171","A_193","A_211","A_304","A_335","A_1600","A_1700","A_4500");
  foreach $wave (@wavelist)
  {
    while(($key,$value) = each (%kw_val_ht))
    {
      if(substr($key,0,5) eq $wave || substr($key,0,6) eq $wave)
      {
        printf( OUTFILE "%-3.3s %-20.20s %-25.25s %-7.7s\n","KWD",$key,$kw_val_ht{$key},"float");
      }
    }
  }
  close(OUTFILE);
 
  system("clear");
  &print_guide();
  print "Step(4)Write values entered to image location file.\n\n";
  print "Completed writing keyword names and values to $image_loc_file.\n\n\n";
}

###############################
# sub print_guide()
###############################
sub print_guide()
{
# begin information on steps on updating master pointer DRMS data series.
print "************************************************************************************************************\n";
print "*   Information and Steps for setting values in master pointer series                                      *\n";
print "*   Step(1)Set T_START, T_STOP, and T_HKVALS Times                                                         *\n";
print "*       ---select default value by selecting enter key or return key or n(no).                             *\n";
print "*       ---enter new values for T_START,T_STOP, T_HKVALS keywords by entering y(yes).                      *\n";
print "*       ---T_START Time is index value used in master pointing series.                                     *\n";
print "*       ---T_STOP Time will be not be calculated automatically or checked. Enter value T_START+6 months.   *\n";
print "*       ---HKVALS Time will be used to get keyword values in drms series for image location keywords.      *\n";
print "*   Step(2)Set HMI Image Location keywords for camera 1 and 2 by entering new value or using default values*\n";
print "*       ---select default value by selecting enter key or return key or n(no).                             *\n";
print "*       ---enter values for each keyword by entering y(yes).                                               *\n";
print "*   Step(3)Set AIA Image Location keywords for each wavelength by entering new value or use default values.*\n";
print "*       ---select default value by selecting enter key or return key or n(no).                             *\n";
print "*       ---enter values for each keyword by entering y(yes).                                               *\n";
print "*       ---wavelength values are for 094,131,171,193,211,304,335,1600,1700 and 4500.                       *\n";
print "*   Step(4)Write values entered to image location file.                                                    *\n";
print "*   Step(5)Optional choice to run load_mp executable using this new image location file.                   *\n";
print "*       ---select y(yes) to run load_mp ilf=<image loc file>.                                              *\n";
print "*       ---select n(no) or hit return to not run load_mp executable and run at command line later after    *\n";
print "*       ---checking values in image location file are correct.                                             *\n";
print "*       ---set the master series value to use in SOURCE_LOAD_MP file.                                      *\n";
print "*  Limitations and Help:                                                                                   *\n";
print "*  (1)To exit script immediately use : cntl-c                                                              *\n";
print "*  (2)Does not check values entered by user at this time including T_STOP                                  *\n";
print "*  (3)The master pointing series setup in script is currently hardcoded to sdo.master_pointing             *\n";
print "*     This is where data is retrieved to display for script and set default values in image location file. *\n";
print "*  (4)Set the master seriesname to load data in SOURCE_LOAD_MP file. This is used by load_mp C exec code.  *\n";
print "*     Item (4) and (3) should have same series name.                                                       *\n";
print "*  (5)May want to add entry field in future for master pointing series name                                *\n";
print "************************************************************************************************************\n\n\n";
}
