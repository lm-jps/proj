#!/usr/bin/perl
#############################################################################
# Name:        do_jsvn_map_file.pl  - Do JSOC version number map file       #
# Description: Creates map of Jsoc series version number to packet version  #
#              number file for each apid. The apid list can be inputed 3    #
#              ways. Need to create help log and checks for arguments.      #
# Execution:   (1) To run :                                                 #
#                  do_jsvn_map_file.pl -f                                   # 
#                  do_jsvn_map_file.pl -g <file version number              # 
#                  do_jsvn_map_file.pl -a <apid number>                     # 
#                  do_jsvn_map_file.pl -h                                   # 
# Logic:       (1)read and process arguments                                #
#              (2)If -g is argument read in folder and read in files in     #
#                 folder to get all apids in filenames.                     #
#              (3)Get list of all initially created JSD files               #
#              (4)Sort throught JSD files and find all with first apid in   #
#                 list of apid to do                                        #
#              (5)Create head of map file                                   #
#              (6)Loop throught list of init JSD files and do differences of#
#                 files and set JSVN_LIST with 1 or increment if different  #
#              (7)Write table, by looping through GTCIDS lines capturing the#
#                 packet-version-number, file-version-number,etc and writing#
#                 line with values from JSVN_LIST.                          #
# Execution:   (1) To run :                                                 #
#                   do_jsvn_map_file.pl                                     # 
# Author:      Carl                                                         #
# Date:        Move from EGSE to JSOC software environment on March 21, 2008#
#############################################################################
# main program                                                              #
#############################################################################
  #check for any arguments passed in command
  &check_arguments();

  ##get environment variables and initialize variables.
  $hm=$ENV{'HOME'};
  $ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
  $ENV{'HK_JSD_PRELIM_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_prelim_file/prod";
  $ENV{'HK_JSVN_MAP_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsn_map_file/prod";
  $ENV{'HK_GTCIDS_FILE'}="gtcids.txt";

  #common setting for all environments
  $ENV{'MAILTO'}="";
  $script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
  $ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";
  # set after set script_dir
  $ENV{'HK_APID_LIST_DAY_FILES'}="$script_dir/hk_apid_list_mapping_jsvn";

  #get list of initial jsoc definition files to use
  &get_list_init_jsds();

  #get list of apids to do map files for
  &get_apid_to_do();

  #go through jsd files looking for first apid
  foreach $apid (@all_apids)
  {
    $strapid= substr($apid, 0, 4);
    print "..doing apid <$strapid>\n";
    foreach $file (@jsdfiles)
    {
      if ($file eq "")
      {
        #print "breaking. final item in list \n";
        break;
      }
      #compare apid value in list with files
      $find_str= sprintf("%s%s",  $strapid, "_");
      if (( index  $file, $find_str ) != -1)
      {
        #if found one push on list of lists
        push( @list, $file);
      }
    }#foreach file
    #Check if have jsd files for this apid
    $list_count=@list;
    if ($list_count == 0)
    {
        print "WARNING: No JSD Files for apid(decimal character value) <$strapid>\n";
        print "Check if HKPDF files exist for this apid and run make_jsd_file.pl to create jsd files\n";
        break;
    }
    else
    {
      $apidfile= shift(@list);
      $firstfilename=$apidfile;
      @jsvn_name = "";
      $jsvn = 1;
      #get hash key from prelim jsd filename's packet version number
      @s_fn= split /\_/, $apidfile;#assume hmi.0001_001032
      $d = substr($s_fn[2],3,3);
      $w = substr($s_fn[2],2,1);
      $hashkey= sprintf("%s.%d",  $w, $d); #hash key is 1.32
      #set first hash key to 1.
      $jsvn_ht{$hashkey}=1;
      # loop thru all prelim jsd files in list for given APID
      while ($apidfile)
      {
        # get file to compare with
        $comparefile=shift(@list);
        #get packet version number to use as hash key in hash table
        @s_fn= split /\_/, $comparefile;#assume hmi.0001_001032
        $d = substr($s_fn[2],3,3);
        $w = substr($s_fn[2],2,1);
        $hashkey= sprintf("%s.%d",  $w, $d); #use 1.32 as hash key 
        # check if there is another file to compare with
        if ( $comparefile eq "" )
        {
          #ends setting jsvn at end of file
          $apidfile = $comparefile;
        }
        # compare files using diff command and set to 1++ is different in hash table
        $difference =`diff -qs  $dir_init_jsds/$apidfile  $dir_init_jsds/$comparefile`;
        if (( index $difference, "differ" ) != -1)
        {
           $jsvn++;
           $jsvn_ht{$hashkey}=$jsvn;
        }
        elsif( (index $difference, "identical" ) != -1)
        {
           $jsvn_ht{$hashkey}=$jsvn;
        }
        elsif( (index $difference, "No such file" ) != -1)
        {
          print "WARNING:no such file ===> $difference";
          print "WARNING:apidfile is $apidfile\n";
        }
        else
        {
          ;
          #print "WARNING:Did not find either differ or identical state ====> $difference\n";
          #print "WARNING:apidfile is $apidfile\n";
        }
        #get next file to compare
        $apidfile=$comparefile;
      }# while loop through all apidfiles in list for this apid find ones that differ      
      $strapid= substr($apid, 0, 4);
      # write header of APID's PVN-TO-JSVN file
      &create_file_header;
      # write each line in APID's PVN-TO-JSVN file for each JSVN number 
      # corresponding to packet version number 
      &get_gtcids_lines;
      #clear out hash table and do next APID's PVN-TO-JSVN file
      while (( $key, $value) = each (%jsvn_ht))
      {
        delete($jsvn_ht{$key});
      }
    }#else-end
  }#for each apid in list
#############################################################################
# subroutine check arguments and set flags                                  #
#############################################################################
sub check_arguments()
{
  $help_flg= "0";
  $apid_list_flg="0";
  if ($#ARGV >= 0)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "help")
    {
      $help_flg = "1";
    }
    elsif ($ARGV[0] eq "-f" )
    {
      #use file to get list of apids to create map files for
      $apid_list_flg = "1";
    }
    elsif ($ARGV[0] eq "-g" )
    {
      #get list of apids to create map files for based on apids in current STHA file version
      $apid_list_flg = "2";
    }
    elsif ($ARGV[0] eq "-a" )
    {
      #use apid following -a to create map files 
      #push decimal character value in this format dddd. Example: -a 0001
      $apid_list_flg = "3";
    }
  }
  else
  {
      $help_flg = "1";
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
  print "(1a)Create map files using -g option will create map files with all apids contained in folder specified :";
  print " do_jsvn_map_file.pl -g <apid-#-version# folder>\n";
  print "    Example: do_jsvn_map_file.pl -g 1.140\n"; 
  print "(1b)Create map files using argument passed for apid number will create map file for that apid only:";
  print " do_jsvn_map_file.pl -a <apid-value in dec.chars>\n";
  print "    Example: do_jsvn_map_file.pl -a 0445\n"; 
  print "(1c)Create map files using list of apids in file, where filename is set in environment variable HK_APID_LIST_DAY_FILES: do_jsvn_map_file.pl -f \n";
  print "    Example: do_jsvn_map_file.pl -f \n"; 
  print "(1d)Get Help Information: do_jsvn_map_file.pl -h  or  do_jsvn_map_file.pl -help\n";
  print "***Note:Currently using option -g\n";
  print "(2)Environment variable HK_APID_LIST_DAY_FILES\n"; 
  print "(2a)Used to list apid to process using argument -f\n";
  print "(2b)The values are in decimal character format.\n";
  print "(2c)If set to 0x1BD 0x1DB then only 2 apid will be processed\n";
  print "(3) Environment variable HK_JSVN_MAP_DIRECTORY\n";
  print "(3a) Used to store location the output files of this script.\n";
  print "(3b) The current format of files is : <APID-Value(dec.chars.)-JSVN-TO_PVN\n";
  print "(4) Environment variable HK_GTCIDS_FILE\n";
  print "(4a) Used to store filename of GROUND TO CODE IDS file to process.\n";
  print "(4b) This file is input data for this script\n";
  print "(4c) The filename format is gtcids.txt.\n";
  print "(5) Environment variable HK_CONFIG_DIRECTORY\n";
  print "(5a) Used to store location of -latest- GROUND TO CODE File(gtcids.txt) file.\n";
  print "(5b) This file is the input for this script.\n";
  print "(6) Environment variable HK_JSD_PRELIM_DIRECTORY\n";
  print "(6a) Used to store location of -initally- created JSOC Series Definition files without headers.\n";
  print "(6b) These files are input for this script to determine differences between initial JSD files.\n";
  print "(6c) Based on differences the JSOC Series Version Number is created in the JSVN Map files.\n";
  exit;
}
#############################################################################
# subroutine get_apid_to_do: gets list of apid to create maps files for     #
#############################################################################
sub get_apid_to_do
{
  # create list of files to process using apid values in file
  my($fn);
  if ($apid_list_flg eq "1")
  {
    #$dn=$ENV{'HK_CONFIG_DIRECTORY'};
    #$dn=".";#Currently, place in directory containing script.
    $fn=$ENV{'HK_APID_LIST_DAY_FILES'};
    open(FILE_APID_LIST, "$fn") || die "(1)Can't Open $fn file: $!\n";
    while (<FILE_APID_LIST>)
    {
      push(@all_apids, $_) ;
    }
    close FILE_APID_LIST ;
  }
  elsif ($apid_list_flg eq "2")
  {
     &get_list_hkpdf_filename;
 
  }
  elsif ($apid_list_flg eq "3")
  {
    #push decimal character value in this format dddd.Example 0001
    push(@all_apids, $ARGV[1]);
    print "...doing map file for apid  $ARGV[1]\n";
  }
  else 
  {
    print "WARNING: Not valid apid list flag value\n";
  }
}
#############################################################################
# subroutine get_list_list_hkpdf_filenames                                  #
#############################################################################
sub get_list_hkpdf_filename()
{
  my($dn);
  $dn= $ENV{'HK_CONFIG_DIRECTORY'};
  #open directory file handle to read in initial jsd files.
  opendir(DIR_HKPDF,"$dn/$ARGV[1]") || die "(2)Can't open folder <$dn$ARGV[1]>:$!\n";
  #get list of hkpdf files
  @hkpdf_filenames=readdir( DIR_HKPDF );
  #close directory file handle
  closedir DIR_HKPSD; 
  push(@hkpdf_filenames, "");

  # get list of apids to do
  (@apidslist)=&get_list_apids_to_do($ARGV[1]);

  # parse apids values from filename to get valid list of apids
  foreach $filename (@hkpdf_filenames)
  {

    # check apid is one we want to create jsd for.These are packets with VER_NUM keywords
    # this skips doing prelim and final jsd that are not in list
    if(substr($filename,0,5) eq "apid-" )
    {
        # get current apid
        $currapid=substr($filename,5,3); 

        # check list of apid want to create JSDs for
        ($foundflg)=&check_apid_list(@apidslist,$currapid);

        # if apid not in list skip creating jsd
        if ($foundflg == 0)
        {
          next;#skip doing
        }
    }

    if (( index  $filename, "apid-" ) != -1)
    {
       
      push(@all_apids, sprintf("%0.4d", hex  substr($filename,5,3)) );
    }
  }
}
#############################################################################
# subroutine get_list_init_jsds: get list of jsd files to check             #
#############################################################################
sub get_list_init_jsds()
{
  #Open input data files
  $dir_init_jsds=$ENV{'HK_JSD_PRELIM_DIRECTORY'};
  #open directory file handle to read in initial jsd files.
  opendir(DIR_JSD, $dir_init_jsds) || die "(3)Can't open:$!\n";
  #get list of jsd files
  @jsdfiles= sort readdir( DIR_JSD );
  #close directory file handle
  closedir DIR_JSD; 
  push(@jsdfiles, "");
}
#############################################################################
# subroutine create file header: Create header of APID-JSVN-TO-PVN file     #
#############################################################################
sub create_file_header
{
  #Open Output file
  $directory=$ENV{'HK_JSVN_MAP_DIRECTORY'};
  #create outfile file name
  $filename= sprintf("%4.4s-%s\n",  $apid, "JSVN-TO-PVN");
  #create file and start writing lines to hdr file
  open(OUTFILE, ">$directory/$filename") || die "(4)Can't Open  $directory/$filename: $!\n" ;
  &get_date;
  &get_gtcids_lookup_column;
  $hexapid = sprintf("%s%03X", "0x",$strapid);
  print OUTFILE "##########################################################################################################################\n";
  printf(OUTFILE  "# FILENAME:                  %s%s                                                                            #\n",
  $strapid,"-JSVN-TO-PVN");
  print OUTFILE "# APID(decimal char. value): $strapid                                                                                        #\n";
  print OUTFILE "# APID(hex. value):          $hexapid                                                                                       #\n";
  print OUTFILE "# PACKET DESCRIPTION:                                                                                                    #\n";
  print OUTFILE "# GTCIDS LOOKUP COLUMN:      $columnid                                                                                      #\n";
  print OUTFILE "# DATE LAST UPDATED:         $theTime                                                                   #\n";
  print OUTFILE "# DESCRIPTION OF COLUMNS:    The JSOC Series and packet version numbers are used to create series name.                  #\n"; 
  print OUTFILE "#                            The File version number,Checkin Date/Time and master file is reference information.         #\n";
  print OUTFILE "##########################################################################################################################\n";
  print OUTFILE "# JSOC Series    | Packet         | File           | Checkin Date/Time   | Master File                                   #\n";
  print OUTFILE "# Version Number | Version Number | Version Number |                     |                                               #\n";
  print OUTFILE "##########################################################################################################################\n";
} 
#############################################################################
# subroutine get_date: lookup date on your system to tag time               #
# created <APID>-JSVN-TO-PVN files                                          #
#############################################################################
sub get_date
{
  #get date created hdr file
  @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
  @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
  ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
  $year = 1900 + $yearOffset;
  $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
}
#############################################################################
# subroutine get_gtcids_lookup_column: Find column to lookup packet version #               
# number in Ground to Code IDS file.                                        #
#############################################################################
sub get_gtcids_lookup_column
{
  if($apid == 5 || $apid == 37 || $apid == 435 || $apid == 465 || $apid == 537 || $apid == 567)
  {
    $columnid ="KERNAL_ID";
  }
  elsif($apid >= 1 && $apid <= 31 )
  {
    $columnid ="HMI_ID";
  }
  elsif($apid >= 420 && $apid <= 499 )
  {
    $columnid ="HMI_ID";
  }
  elsif($apid >= 32 && $apid <= 63 )
  {
    $columnid ="AIA_ID";
  }
  elsif($apid >= 500 && $apid <= 599 )
  {
    $columnid ="AIA_ID";
  }
}
#############################################################################
# sub get_gtcids_lines: gets  lines from gtcids.txt file                    #
#############################################################################
sub get_gtcids_lines
{
  my($dn,$fn,@all_gtcids_lines);
  #open GTCIDS file
  $dn_fn = sprintf("%s/%s",  $ENV{'HK_CONFIG_DIRECTORY'}, $ENV{'HK_GTCIDS_FILE'});
  open(FILE, "$dn_fn") || die "(5)Can't Open: $!\n";
  $i=1;
  #loop throught lines in GTCIDS file
  while (<FILE>)
  {
    push(@all_gtcids_lines, $_) ;
    @s_line= split / \s*/,  $_;
    if ( substr($s_line[0],0,1) eq '#')
    {
      next; #skip
    }
    else #for HMI column id only-need to fix!
    {
      #load key parameters from each line to use in APID-JSVN-PVN map file
      $file_version_number = $s_line[13];
      $stha_ci_date = $s_line[0];
      $stha_ci_time = $s_line[1];
      $stha_gen_time = $s_line[15];
      if ( $columnid eq "HMI_ID") 
      {
        $packet_version_number = $s_line[7];
      }
      elsif ( $columnid eq "AIA_ID" ) 
      {
        $packet_version_number = $s_line[9];
      }
      elsif ( $columnid eq "KERNAL_ID")
      {
        $packet_version_number = $s_line[5];
      }
      else
      {
         print "WARNING: Did not find columnid type <$columnid>\n";
      }

      # calculate start point for jsvn
      @s_ffn= split /\_/, $firstfilename;
      $startjsvn = substr($s_ffn[3],3,3);
      #skip loading file version number below 1.32- 1.32 is based esblished.
      $str_fvn = sprintf("%s",$file_version_number);
      @s_fvn = split /\./,  $str_fvn;
      #start jsvn at starting packet version number of prelim jsd files 
      $str_fvn = sprintf("%s",$file_version_number);
      $str_pvn = sprintf("%s",$packet_version_number);
      @s_pvn = split /\./,  $str_pvn;
      # start adding lines from file version 1.32
      if ( int $s_fvn[1] < 32)
      {
        next; #skip
      }
      #if jsvn hash table has this packet version number as key then have a JSVN.
      elsif ( $jsvn_ht{$packet_version_number} )
      {
        #create line for apid-JSVN-TO-PVN file
        $jsn_map_line = sprintf(" %15.4d | %14s | %14s | %10s %8s | %47s", int $jsvn_ht{$packet_version_number},  $packet_version_number,  $file_version_number, $stha_ci_date, $stha_ci_time, $stha_gen_time);
        print OUTFILE "$jsn_map_line";
        $i++;
      }
      #else just set the JSVN to 0000
      else
      {
        #create line for apid-JSVN-TO-PVN file
        $jsn_map_line = sprintf(" %15.4d | %14s | %14s | %10s %8s | %47s", 0000,  $packet_version_number,  $file_version_number, $stha_ci_date, $stha_ci_time, $stha_gen_time);
        print OUTFILE "$jsn_map_line";
      }
      $jsn_map_line="";
    }#else
  }#end while
  close(FILE);
  close(OUTFILE);
}
###########################################
#  check apid list to do                  #
###########################################
sub check_apid_list (@$)
{
  $foundflg=0;
  foreach $a (@apidslist)
  {
    if($currapid =~ m/$a/i)
    {
      $foundflg=1;
      last;
    }
  }
 return ($foundflg);
}
###########################################
#  get list of apids to do                #
###########################################
sub get_list_apids_to_do($)
{
  # misses VER_TEMPERATURE
  # $files=`cd /home1/carl/cvs/TBL_JSOC/lev0/hk_config_file/$file_version/; grep VER_NUM apid*`;
  # get all apids that have VER_NUM or VER_TEMPERATURE Keyword
  my $files=`cd $ENV{'HK_CONFIG_DIRECTORY'}/$ARGV[1]/;  egrep '(VER_NUM|VER_TEMPERATURE)' apid*`;

  # remove apid text regular expression
  $files =~ s/(a.+?-)//g;

  # remove everything from after apid number to end of line regular expression
  $files =~ s/(-.+?\n)/ /g;

  # split apids into separate field in array
  my @apidslist = split(/ /, $files);

  # return list of apids to do
  return ( @apidslist);
}
