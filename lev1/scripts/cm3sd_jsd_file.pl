#!/usr/bin/perl
##############################################################################
# Name:        cm3sd_jsd_file.pl - Make JSD File                             #
# Description: create jsd file for mean,min,max and standard deviation values#
# Process:     (1)Create instruction file.                                   #
#              (2)Run cm3sd_jsd_file.pl to make jsd file                     #
#              (3)Create jsd series by  using create_series command          #
# Execution:   (1)To run :                                                   #
#              cm3sd_jsd_file.pl inf=<durectory to instruction file>         #
#              (2)For help: cm3sd_jsd_file.pl  -h                            #
# Example:     cm3sd_jsd_file.pl inf=./apid19_temp_instruct_file             #
# Environment Variable: Set HK_CM3SD_JSD_DIRECTORY variable were want        #
#                       jsd to go to check in to cvs if needed. See below    #
# Requirement: Adjust the cvsver variable below when change script and check #
#              into cvs in order to track which version of script created    #
#              jsd files since version is written to jsd files created.      #
# Author:      Carl                                                          #
# Date:        Created August,21, 2009                                       #
##############################################################################
# main function                                                              #    
##############################################################################
# get environment variables and initialize variables.
# use perl -d <script file> to debug
$hm=$ENV{'HOME'};
$scriptname="cm3sd_jsd_file.pl";
$cvsver="1.0"; ##version needs to change as updates script##

#PRODUCTION SETTINGS to send jsd files to production directory
$ENV{'HK_CM3SD_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev1/hk_jsd_file/prod";

#USER SETTING to send jsd files to user directory 
#turn off production by commenting out, turn on what is need below
#su_carl setting for test creating instruction file.
#$ENV{'HK_CM3SD_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev1/hk_jsd_file/su_carl";
#su_rock setting for test creating instruction file.
#$ENV{'HK_CM3SD_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev1/hk_jsd_file/su_rock";

#local variables for main program
my($ins_filename,$templatefn,$seriesname,$author,$owner,$description,$interval,@kw,@all_m3sd_keywords);

# (1)check arguments used are correct 
$ins_filename=&check_arguments(@ARGV );
print ". . . . using instruction filename <$ins_filename> to create new jsd file.\n";

# (2) check instruction file does not have tabs or DOS special chararters
($curr_ifn)=&check_instruction_file($ins_filename);

# (3)read template instruction file
($templatefn,$seriesname,$author,$owner,$description,$interval, @kw)=&read_instruction_file($curr_ifn);

# (4)create mean,min,max, and standard dev. keywords for each keyword
@all_m3sd_keywords=&create_m3sd_keywords(@kw);

# (5)write jsd file
&write_m3sd_keywords($curr_ifn,$seriesname,$description,$author,$owner,$interval,@kw);




###############################################################################
# sub create_m3sd_keywords 
###############################################################################
sub create_m3sd_keywords(@)
{

  # local variable and setting variable using passed argument
  my(@keywords) = @_;

  #loop thru array of keywords to make m3sd keywords
  $j=0;
  foreach $i (@keywords)
  {
    if($i ne "")
    {
      # get keyword  short name were line:keyword:<apid>,<long name>,<short name>
      @s_kwl=split(',',$i);
      $s_kwl[2] =~ s/\n//g; #remove cr -regular exp

      # make m3sd keyword
      push(@m3sd_kw_list, sprintf("%s_max",$s_kwl[2]));
      push(@m3sd_kw_list, sprintf("%s_min",$s_kwl[2]));
      push(@m3sd_kw_list, sprintf("%s_mean",$s_kwl[2]));
      push(@m3sd_kw_list, sprintf("%s_sd",$s_kwl[2]));

      # push(@mkeywords, $new_kw )
      # return @mkeyword list
    }
    $j++;
  }
  return (@m3sd_kw_list);
}
###############################################################################
# sub  write_m3sd_keywords - write max,min,mean, and standard deviation keywords
###############################################################################
sub write_m3sd_keywords($,$,$,$,$,$,@)
{
  # local variables
  my($inst,$sn,$dc,$au,$ow,$intval,@kwds, $sn_fn,$date);

  # set local variables to passed arguments
 ($inst,$sn,$dc,$au,$ow,$intval,@kwds)=@_;

  # get date
  $date=`date`;
  $date=~ s/\n//g;

  # get directory to write jsd file
  $directory=$ENV{'HK_CM3SD_JSD_DIRECTORY'};

  # create seriesname filename
  $sn_fn=sprintf("%s%s",$sn,".jsd");

  # Create directory using environment variable, version number and filename
  open(OUTFILE, ">$directory/$sn_fn") || die "(8)Can't Open  $directory/$sn_fn: $!\n" ;
  print ". . . . creating jsd filename <$directory/$sn_fn>\n";

  #set other required values for jsd
  $unitsize="1";
  $archive= "0";
  $index="T_START";
  $tapegroup= "0";
  $retention= "0";

  # Create global section of file
  printf(OUTFILE "#================================================================================================\n");
  printf(OUTFILE "# Maximum, Minimum, Mean and Standard Deviation JSD \n");
  printf(OUTFILE "# Created by script: $scriptname  cvs version:$cvsver\n");
  printf(OUTFILE "# Instruction file used: $inst\n");
  printf(OUTFILE "# Date created: $date\n");
  printf(OUTFILE "#================================================================================================\n");
  printf(OUTFILE "#================================== Global Series Information ===================================\n");
  printf(OUTFILE  "SeriesName:  %-s\n", $sn);
  print OUTFILE   "Description:$dc \n" ;
  printf(OUTFILE  "Author:      %-s\n", $au);
  printf(OUTFILE  "Owner:       %-s\n", $ow);
  printf(OUTFILE  "Unitsize:    %-s\n", $unitsize);
  printf(OUTFILE  "Archive:     %-s\n", $archive);
  printf(OUTFILE  "Retention:   %-s\n", $retention);
  printf(OUTFILE  "Tapegroup:   %-s\n", $tapegroup);
  printf(OUTFILE  "Index:       %-s\n", $index);

  # Create keyword section of file
  printf(OUTFILE "#================================ Keywords Series Information ===================================\n");

  # Get data values for T_START keyword
  $keyword_name = "T_START";
  $type_value = "time";
  $print_value="\"2\"";#change 9-1
  $default_value = "TSEQ_EPOCH";
  $unit_value="\"UTC\"";#changed 9-1
  $comment_value = "\"T_START\"";

  # Create line for T_START
  $tmp_line = sprintf("%s%20.20s,%9.9s,%12.12s,%8.8s,%20.20s,%6.6s,%12.12s,  %-s\n",
              "Keyword:", $keyword_name, $type_value,
              "ts_eq", "record", $default_value,
              $print_value, $unit_value, $comment_value);
  printf(OUTFILE "%s", $tmp_line);

  # Get data values for T_START epoch
  $keyword_name = "T_START_epoch";
  $type_value = "time";
  $print_value="\"2\"";#changed 9-1
  $default_value = "TSEQ_EPOCH";
  $unit_value="\"UTC\"";#changed 9-1
  $comment_value = "\"T_START_epoch\"";

  # Create and write line for T_START_epoch
  $tmp_line = sprintf("%s%20.20s,%9.9s,%12.12s,%8.8s,%20.20s,%6.6s,%12.12s,  %-s\n",
              "Keyword:", $keyword_name, $type_value,
              "constant", "record", $default_value,
              $print_value, $unit_value, $comment_value);
  printf(OUTFILE "%s", $tmp_line);

  # Get data for T_START step
  $keyword_name = "T_START_step";
  $type_value = "float";
  $print_value="\"%f\"";
  $default_value = $intval;
  $unit_value="\"$intval sec\"";
  $comment_value = "\"T_START_step\"";

  #Create and write line for T_START_step
  $tmp_line = sprintf("%s%20.20s,%9.9s,%12.12s,%8.8s,%20.20s,%6.6s,%12.12s,  %-s\n",
              "Keyword:", $keyword_name, $type_value,
              "constant", "record", $default_value,
              $print_value, $unit_value, $comment_value);
  printf(OUTFILE "%s", $tmp_line);

  # Create Number of points per block or interval(T_START to next T_START)
  $keyword_name = "NUMPTS";
  $type_value = "\"int\"";
  $print_value= "\"%d\"";
  $default_value = "DRMS_MISSING_VALUE";
  $unit_value= "\"none\"";
  $comment_value = "\"NUMPTS - number of points in interval\"";
  $tmp_line = sprintf("%s%20.20s,%9.9s,%12.12s,%8.8s,%20.20s,%6.6s,%12.12s,  %-s\n",
              "Keyword:", $keyword_name, $type_value,
              "variable", "record", $default_value,
              $print_value, $unit_value, $comment_value);
  printf(OUTFILE "%s", $tmp_line);


  # create min,max,mean and sd keywords using keywords in instruction file
  foreach $kwd_line (@kwds)
  {
    if($kwd_line eq "")
    {
      next; #skip lines with no data
    }
    # Get data for each keyword line
    $type_value = "\"float\"";
    $print_value= "\"%f\"";
    $default_value ="DRMS_MISSING_VALUE";
    $unit_value= "\"none\"";

    # use long name in comments
    # get keyword  short name were line:keyword:<apid>,<long name>,<short name>
    @s_kwl=split(',',$kwd_line);
    $s_kwl[1] =~ s/\n//g; #remove cr -regular exp
    $s_kwl[1] =~ s/ //g; #remove cr -regular exp
    $s_kwl[2] =~ s/\n//g; #remove cr -regular exp
    $s_kwl[2] =~ s/ //g; #remove cr -regular exp
    $comment_value = sprintf("\"%s\"", $s_kwl[1]);

    # create line to print to jsd file
    foreach $suffix ("MAX","MIN","MEAN","SD")
    {
      # add suffix MAX,MIN,MEAN,SD etc
      $kw_long=sprintf("%s - %s",$s_kwl[1],$suffix);
      $kw_short=sprintf("%s_%s",$s_kwl[2],$suffix);

      # take out tabs and blanks
      $kw_short =~ s/^\s+//;
      $type_value =~ s/^\s+//;
      $print_value=~ s/^\s+//;
      $default_value=~ s/^\s+//;
      $unit_value=~ s/^\s+//;

      # create line for keyword in jsd file
      $tmp_line = sprintf("%s%20.20s,%9.9s,%12.12s,%8.8s,%20.20s,%6.6s,%12.12s,  \"%-s\"\n",
      "Keyword:", $kw_short, $type_value, "variable", "record", $default_value, $print_value, $unit_value, $kw_long);

      # print keyword line to jsd file
      printf(OUTFILE "%s", $tmp_line);
    }
  }
  close (OUTFILE);
  print ". . . . completed creating jsd filename <$directory/$sn_fn>\n";
  print ". . . . please create series in DRMS using: % create_series  $directory/$sn_fn\n";

}



###############################################################################
# sub  read_instruction_file
###############################################################################
sub read_instruction_file($)
{
    
    #local variables
    my($ifn);

    # set argument passed to local variable
    $ifn = $_[0];   #instruction filename

    #initialize array variable to null
    @all_kwd_lines="";

    # open one instruction file
    open(FILE, "$ifn") || die "(6)Can't Open $ifn file: $!\n";
    while (<FILE>)
    {
      @s_line= split / \s*/,  $_;

      if ( substr($s_line[0],0,1) eq "#")
      {
        ;#skip
      }
      elsif ( substr($s_line[0],0,7) eq "Keyword" || substr($s_line[0],0,7) eq "keyword")
      {
        push(@all_kwd_lines, $_) ;
      }
      elsif ( substr($s_line[0],0,12) eq "templatename" || substr($s_line[0],0,12) eq "TemplateName")
      {
        $tplname=$_ ;
      }
      elsif ( substr($s_line[0],0,10) eq "seriesname" || substr($s_line[0],0,10) eq "SeriesName")
      {
        $ser_name=$_ ;
        $ser_name=~ s/([a-z]*:)//g; #remove up to : using regular expression
        $ser_name =~ s/\n//g; #remove cr -regular exp
      }
      elsif ( substr($s_line[0],0,11) eq "description" || substr($s_line[0],0,11) eq "Description")
      {
        $desc=$_ ;
        $desc=~ s/([a-z]*:)//g;   #remove up to : using regular expression
        $desc =~ s/\n//g; #remove cr -regular exp
      }
      elsif ( substr($s_line[0],0,6) eq "author" || substr($s_line[0],0,6) eq "Author")
      {
        $auth=$_ ;
        $auth=~ s/([a-z]*:)//g;   #remove text up to : using regular expression
        $auth =~ s/\n//g; #remove cr -regular exp
      }
      elsif ( substr($s_line[0],0,5) eq "owner" || substr($s_line[0],0,5) eq "Owner")
      {
        $own=$_ ;
        $own=~ s/([a-z]*:)//g;    #remove up text to : using regular expression
        $own =~ s/\n//g; #remove cr -regular exp
      }
      elsif ( substr($s_line[0],0,8) eq "interval" || substr($s_line[0],0,8) eq "Interval")
      {
        $ivalue=$_ ;
        $ivalue=~ s/([a-z]*:)//g; #remove text up to : using regular expression
        $ivalue =~ s/\n//g; #remove cr -regular exp
      }
      else
      {
        print "Warning:skipping value of line from file:  <@s_line>\n";
      }
      next;
    }#end-while
    close( FILE);
    return ($tplname, $ser_name, $auth,$own,$desc,$ivalue, @all_kwd_lines);
}



####################################################
#  check arguments
####################################################
sub check_arguments(@)
{
  #declare local variable and set to argument passed to function
  my ( @array) = @_; 

  if (($#array != 0)) 
  {
    print "Waring: not correct usage!\nUsage: $0  isf=<absolute path to instruction-file>\nUsage: $0  -h \n";
    die "exiting $!";
  }

  if ( $array[0] eq "-h" )
  {

    print "Help Listing\n";
    print "(1)Ways to Execute Perl Script: \n";
    print "(1a)Create Min, Max, Mean and Standard deviation jsd file, \n";
    print "    where isf=<absolute path to instruction filename>:\n";
    print "            cm3sd_jsd_file.pl isf= <instruction-file>\n\n";
    print "(1b)Get Help Information:  cm3sd_jsd_file.pl -h \n";
    print "            cm3sd_jsd_file.pl -h \n";
    print "(2)Example Execution of  Perl Script: \n";
    print "(2a)% cm3sd_jsd_file.pl isf=/home/carl/cvs/TBL_JSOC/lev1/instruction_file/su_carl/hmitest1200_thermal_template.txt\n";
    exit;
  }
  else 
  {
    @s_arg=split('\=', $array[0] );
    #print "s_arg-0 is <$s_arg[0]>\n";
    #print "s_arg-1 is <$s_arg[1]>\n";

    #check argument passed looks okay
    if ($s_arg[0] eq "isf" && length($s_arg[1]))
    {
       if (-e $s_arg[1] && not -z $s_arg[1])
       {
          return($s_arg[1]);
       }
       else
       {
         print "ERROR: Got a file <$s_arg[1]> that is zero size or does not exist. Exiting.\n";
         die "Usage: $0 < isf=<absolute path to instruction-file> :$!";
       }
    }
    else
    {
       die "Usage: $0 < isf=<absolute path to instruction-file> :$!";
    }
  }
  return("");
}

####################################################
#  check instruction file
####################################################
sub check_instruction_file($)
{
    
    #local variables
    my($ifn, $redo_ifn);

    # set argument passed to local variable
    $ifn = $_[0];   #instruction filename

    #initialize array variable to null
    @all_kwd_lines="";
    $redo_ifn=0;

    # open instruction file and check if \r\c or \t in file 
    # remove these character since c code lm3s does not like these
    open(FILE, "$ifn") || die "(6)Can't Open $ifn file: $!\n";
    while (<FILE>)
    {
      @s_line= split / \s*/,  $_;
      if($_ =~ /\t/)
      {
         $fix_to_prt=$_;
         $fix_to_prt=~ s/\r\n//g;
         print ". . . . compiling Instruction file:Found Error - tab(s) in line:<$fix_to_prt>\n";
         $redo_ifn=1;
      }
      if($_ =~ /\r\n/)
      {
         $fix_to_prt=$_;
         $fix_to_prt=~ s/\r\n//g;
         print ". . . . compiling Instruction file:Found Error- DOS Return characters in line:<$fix_to_prt>\n";
         $redo_ifn=1;
      }
      next;
    }#end-while
    close( FILE);

    #remove tabs and \r from instruction files created in dos.
    if($redo_ifn)
    {
      $oldifn=sprintf("%s%s",$ifn,"-OLD");
      print ". . . . when compiling Instruction file: Found Errors with instruction file:<$ifn>\n";
      print ". . . . copying  old or orginal instruction file to filename:<$oldifn>. View your orginal file here.\n";
      $log=`cp $ifn $oldifn`;
      print ". . . . new fixed instruction file's filename  <$ifn> is used to create jsd.\n";
      open(NEWFILE, ">$ifn") || die "(6)Can't Open $newifn file: $!\n";
      # open old instruction file
      open(OLDFILE, "$oldifn") || die "(6)Can't Open $ifn file: $!\n";
      while (<OLDFILE>)
      {
        $newline=$_;
        if($_ =~ /\t/)
        {
          $newline =~ s/\t//g;
        }
       if($_ =~ /\r\n/)
       {
          $newline =~ s/\r//;
       }
       print NEWFILE $newline;
       next;
      }
      close(OLDFILE);
      close(NEWFILE);
      $current_instruction_filename=$ifn;
    }
    else
    {
      $current_instruction_filename=$ifn;
    }
    return ($current_instruction_filename);
}

