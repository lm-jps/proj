#!/usr/bin/perl
#############################################################################
# Name:        jsoc_do_jsvn_map_file.pl  - Do JSOC version number map file  #
# Description: Creates map of Jsoc series version number to packet version  #
#              number file for each HK apid(1-63,400's, and 500's). Also    #
#              create map of Jsoc series version number to file version     #
#              number for each SDO HK APID(129).The apid list can be        #
#              inputed 3 ways. Need to create help log and checks for       #
#              arguments.                                                   #
# Execution:   (1) To run :                                                 #
#                jsoc_do_jsvn_map_file.pl -f                                # 
#                jsoc_do_jsvn_map_file.pl -g <file version number>  <HK|SDO># 
#                jsoc_do_jsvn_map_file.pl -a <apid number>                  # 
#                jsoc_do_jsvn_map_file.pl -h                                # 
#                                                                           #
# Logic:       (1)read and process arguments                                #
#              (2)If -g is argument read in folder and read in files in     #
#                 folder to get all apids in filenames.                     #
#              (3)Get list of all initially created JSD files               #
#              (4)Sort throught JSD files and find all with first apid in   #
#                 list of apid to do                                        #
#              (5)Create head of map file                                   #
#              (6)Loop throught list of init JSD files and do differences of#
#                 files and set JSVN_LIST with 1 or increment if different  #
#              (7)Write map file, by looping through GTCIDS or SHCIDS       #
#                 files' s lines capturing the packet-version-number,       #
#                 file-version-number,etc and writing lines with values     #
#                 from JSVN_LIST.                                           #
# Example Execution:                                                        #
#               (1)jsoc_do_jsvn_map_file.pl  -a 0129                        # 
#               (2)jsoc_do_jsvn_map_file.pl  -a 0445                        # 
#               (3)jsoc_do_jsvn_map_file.pl  -f                             # 
#               (4)jsoc_do_jsvn_map_file.pl  -g 1.163 HK                    # 
#               (5)jsoc_do_jsvn_map_file.pl  -g 1.2 SDO                     # 
#                                                                           #
# Author:      Carl                                                         #
# Date:        Move from EGSE to JSOC software environment on March 21, 2008#
#############################################################################
# main program                                                              #
#############################################################################

#check for any arguments passed in command
&check_arguments();

#get environment variables and initialize variables.
$hm=$ENV{'HOME'};
$ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
$ENV{'HK_JSD_PRELIM_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_prelim_file/prod";
$ENV{'HK_JSVN_MAP_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsn_map_file/prod";
$ENV{'HK_GTCIDS_FILE'}="gtcids.txt";
$ENV{'HK_SH_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/sdo_hk_config_file";
$ENV{'HK_SHCIDS_FILE'}="shcids.txt";

#common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";

#set after set script_dir - 
#contains the listing of apids to process when use -f option.
#the values are one apid per line in file.
#the values should be in decimal format(i.e.,0445,0475,0129,etc) 
$ENV{'HK_APID_LIST_DAY_FILES'}="$script_dir/hk_apid_list_mapping_jsvn";

#get list of initial jsoc definition files to use
&get_list_init_jsds();

#get list of apids to do map files for
&get_apid_to_do();

#go through jsd files looking for first apid
foreach $apid (@all_apids)
{
  $strapid= substr($apid, 0, 4);
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
    print "WARNING:jsoc_do_jsvn_map_file.pl: No JSD Files for apid(decimal character value) <$strapid>\n";
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
    @s_fn= split /\_/, $apidfile;#ASSUME using this format hmi.lev0_0001_001032
    if ($s_fn[1] eq "ground.lev0")
    {
      print "ERROR: jsoc_do_jsvn_map_file.pl:Sorry script is assuming jsd names are hmi or aia or sdo format but might be getting <name>_ground!\n";
      print "Exiting script. Adjust code to handle this case. Note jsd name using after split is:@s_fn\n";
      exit;
    }

    # check if processing SDO or HK apids
    if( $apid > 96 &&  $apid < 400) 
    {
      #for sdo apids, use start packet date range as hash value
      $year = substr($s_fn[2],0,4);
      $month = substr($s_fn[2],4,2);
      $day = substr($s_fn[2],6,2);
      $hashkey= sprintf("%s%s%s",  $year, $month,$day); #hash key is 20080914
    }
    else
    { 
      #for hk apids, use packet version number as hash value
      $d = substr($s_fn[2],3,3);
      $w = substr($s_fn[2],2,1);
      $hashkey= sprintf("%s.%d",  $w, $d); #hash key is 1.32
    }

    #set first hash key to 1.
    $jsvn_ht{$hashkey}=1;

    # loop thru all prelim jsd files in list for given APID
    while ($apidfile)
    {
      # get file to compare with
      $comparefile=shift(@list);
    
      #get packet version number to use as hash key in hash table
      @s_fn= split /\_/, $comparefile;#assume hmi.0001_001032
      if( $apid > 96 &&  $apid < 400) 
      {
        #for sdo apids, use start packet date range as hash value
        $year = substr($s_fn[2],0,4);
        $month = substr($s_fn[2],4,2);
        $day = substr($s_fn[2],6,2);
        $hashkey= sprintf("%s%s%s",  $year, $month,$day); #hash key is 20080914
      }
      else
      {
        $d = substr($s_fn[2],3,3);
        $w = substr($s_fn[2],2,1);
        $hashkey= sprintf("%s.%d",  $w, $d); #use 1.32 as hash key 
      }
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
        print "WARNING:jsoc_do_jsvn_map_file.pl: no such file ===> $difference";
        print "WARNING:jsoc_do_jsvn_map_file.pl: apidfile is $apidfile\n";
      }
      else
      {
        ;
        #print "WARNING:jsoc_do_jsvn_map_file.pl: Did not find either differ or identical state ====> $difference\n";
        #print "WARNING:jsoc_do_jsvn_map_file.pl: apidfile is $apidfile\n";
      }
      #get next file to compare
      $apidfile=$comparefile;
    }# while loop through all apidfiles in list for this apid find ones that differ      

    # write header of APID's PVN-TO-JSVN file
    $strapid= substr($apid, 0, 4);
    &create_file_header($strapid);

    # write each line in APID's PVN-TO-JSVN file for each JSVN number 
    # corresponding to packet version number 
    $apid =~ s/\n//g; #regular exp rm cr 
    if( $apid > 96 &&  $apid < 400) 
    {
      print "-->creating map file for apid:<$apid>\n";
      &get_shcids_lines;
    }
    elsif( $apid < 64 ||  $apid > 399)
    {
      print "-->creating map file for apid:<$apid>\n";
      &get_gtcids_lines;
    }
    else
    {
      $ha= hex $apid;
      print "ERROR:jsoc_do_jsvn_map_file.pl: apid decimal value unknown:<$apid>. The apid hex value is:<$ha>. Exiting script\n";
      exit;
    }

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
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help")
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
      if ($#ARGV != 2 or $ARGV[1] eq ""  or $ARGV[2] eq "" or ($ARGV[2] ne "HK" and  $ARGV[2] ne "SDO"))
      {
        print "\nERROR:jsoc_jsvn_map_file.pl: Did not enter correct arguments for -g flag:<$ARGV[0] $ARGV[1] $ARGV[2]>\n\n";
        &show_help_info("USAGE");
      }
    }
    elsif ($ARGV[0] eq "-a" )
    {
      #use apid following -a to create map files 
      #push decimal character value in this format dddd. Example: -a 0001
      $apid_list_flg = "3";
      if ($ARGV[1] eq "" or $ARGV[1] eq "0")
      {
        print "\nERROR:jsoc_jsvn_map_file.pl: Did not enter apid value! Enter 4 digit decimal value for apids.\n\n";
        &show_help_info("USAGE");
      }
      else
      {
        ;#print "got value:$ARGV[1]\n";
      }
      # check argument passed for -a is 4 digits
      # regular expression -check got exactly 4 digits
      if($ARGV[1] =~ m/^\d{4}/)
      {
          ;#match-got 4 digits
      }
      else
      {
        print "\nERROR:jsoc_jsvn_map_file.pl: Did not enter correct apid value:<$ARGV[1]>! Enter exactly 4 digit decimal value for apids including leading zeros(i.e.,0129, 0445 or 0001).\n\n";
        &show_help_info("USAGE");
      }
    }
  }
  else
  {
      $help_flg = "1";
  }
  if ( $help_flg eq "1")
  {
     &show_help_info("HELP");
  }
}

#############################################################################
# subroutine show_help_info: show help information                          #
#############################################################################
sub show_help_info($)
{
  # get arg passed
  $argpassed= $_[0];

  if ($argpassed eq "HELP")
  {
    print "Help Listing\n";
  }
  else
  {
    print "Usage:\n";
  }

  if ($argpassed eq "USAGE" or $argpassed eq "HELP")
  {
    print "(1)Ways to Execute Perl Script: \n";
    print "(1a)Create map files using -g option will create map files with all apids contained in folder specified :\n";
    print "    jsoc_do_jsvn_map_file.pl -g <apid-#-version# folder> < SDO | HK >\n";
    print "    Example: jsoc_do_jsvn_map_file.pl -g 1.140 HK\n\n"; 
    print "(1b)Create map files using argument passed for apid number will create map file for that apid only:\n";
    print "    jsoc_do_jsvn_map_file.pl -a <apid-value in dec.chars>\n";
    print "    Example: jsoc_do_jsvn_map_file.pl -a 0445\n\n"; 
    print "(1c)Create map files using list of apids in file:\n    jsoc_do_jsvn_map_file.pl -f \n";
    print "    Where filename is set in environment variable HK_APID_LIST_DAY_FILES in this script. \n";
    print "    Example: jsoc_do_jsvn_map_file.pl -f \n\n"; 
    print "(1d)Get Help Information: jsoc_do_jsvn_map_file.pl -h  or  jsoc_do_jsvn_map_file.pl -help\n\n";
  }
  if ($argpassed eq "HELP")
  {
    print "***Note:Currently using option -g. This option does many at once.This option filters out creating\n";
    print "***map file for HK config packets(not sdo hk packets) that do not have a VER keyword.\n\n";
    print "More information on Environment variables\n"; 
    print "(2)Environment variable HK_APID_LIST_DAY_FILES\n"; 
    print "(2a)Used to list apid to process using argument -f\n";
    print "(2b)The values are in decimal character format.\n";
    print "(2c)If set to 0445 and  0475 then only 2 apids will be processed\n";
    print "(3) Environment variable HK_JSVN_MAP_DIRECTORY\n";
    print "(3a) Used to store location to put output map files of this script.\n";
    print "(3b) The current format of files is : <APID-Value(dec.chars.)-JSVN-TO_PVN\n";
    print "(4) Environment variable HK_GTCIDS_FILE and HK_SHCIDS_FILE\n";
    print "(4a) Used to store filename of GROUND TO CODE IDS file to process.\n";
    print "(4b) Used to store filename of SDO HK CODE IDS file to process.\n";
    print "(4c) These files are used as input data for this script\n";
    print "(4d) The filename formats are currently gtcids.txt and shcids.txt.\n";
    print "(5) Environment variable HK_CONFIG_DIRECTORY and HK_SH_CONFIG_DIRECTORY\n";
    print "(5a) Used to store location of -latest- GROUND TO CODE File(gtcids.txt) and shcids.txt files.\n";
    print "(5b) These files are the input for this script.\n";
    print "(6) Environment variable HK_JSD_PRELIM_DIRECTORY\n";
    print "(6a) Used to store location of -initally- created JSOC Series Definition files without headers.\n";
    print "(6b) These files are input for this script to determine differences between initial JSD files.\n";
    print "(6c) Based on differences the JSOC Series Version Number is created in the JSVN Map files.\n";
  }
  exit;
}

#############################################################################
# subroutine get_apid_to_do: gets list of apid to create maps files for     #
#############################################################################
sub get_apid_to_do
{
  # create list of files to process using apid values in file
  my($fn,$found_dup_flg);
  if ($apid_list_flg eq "1") # -f flag
  {
    #$dn=$ENV{'HK_CONFIG_DIRECTORY'};

    #$dn=".";#Currently, place in directory containing script.
    $fn=$ENV{'HK_APID_LIST_DAY_FILES'};
    open(FILE_APID_LIST, "$fn") || die "(1)Can't Open $fn file: $!\n";
    while (<FILE_APID_LIST>)
    {
      #convert apid string decimal value to hex string for check_dup_merged_apid function
      $hexstr= sprintf("%x", $_);
      ($found_dup_flg)= &check_dup_merged_apid( $hexstr);        
      if ($found_dup_flg)
      { 
        next;#skip doing
      }
      push(@all_apids, $_) ;
    }
    close FILE_APID_LIST ;
  }
  elsif ($apid_list_flg eq "2") # -g flag
  {
     &get_list_hkpdf_filename($ARGV[1],$ARGV[2]);
  }
  elsif ($apid_list_flg eq "3") # -a flag
  {
    #convert apid string decimal value to hex string for check_dup_merged_apid function
    $hexstr= sprintf("%x", $ARGV[1]);
    #skip doing if duplicate of merged apid 445,529,etc
    ($found_dup_flg)= &check_dup_merged_apid( $hexstr);        
    if ($found_dup_flg)
    { 
      #skip doing
      print "-->skipping doing map file for apid  $ARGV[1]\n";
    }
    else
    {
      #push decimal character value in this format dddd.Example 0001
      push(@all_apids, $ARGV[1]);
      print "-->doing map file for apid  $ARGV[1]\n";
    }
  }
  else 
  {
    print "WARNING: jsoc_do_jsvn_map_file.pl: Not valid apid list flag value\n";
  }
}

#############################################################################
# subroutine get_list_list_hkpdf_filenames                                  #
#############################################################################
sub get_list_hkpdf_filename($$)
{
  # declare locals and set args
  my($dn,$arg_fvn,$arg_processing_type, $found_dup_flg);

  # aguments passed
  $arg_fvn=$_[0]; #first argument passed file version number
  $arg_processing_type=$_[1];#second argument passed is process type flag

  #set config directory based on processing type
  if($arg_processing_type eq "SDO" )
  {
    $dn= $ENV{'HK_SH_CONFIG_DIRECTORY'};
  }
  elsif ($arg_processing_type eq "HK" )
  {
    $dn= $ENV{'HK_CONFIG_DIRECTORY'};
  }

  #open directory file handle to read in initial jsd files.
  opendir(DIR_HKPDF,"$dn/$arg_fvn") || die "(2)Can't open folder <$dn/$arg_fvn>:$!\n";

  #get list of hkpdf files
  @hkpdf_filenames=readdir( DIR_HKPDF );

  #close directory file handle
  closedir DIR_HKPSD; 
  push(@hkpdf_filenames, "");

  # get list of apids to do
  (@apidslist)=&get_list_apids_to_do($arg_fvn,$arg_processing_type);

  # parse apids values from filename to get valid list of apids 
  foreach $filename (@hkpdf_filenames)
  {

    # check apid is one we want to create jsd for.These are packets with VER_NUM keywords
    # this skips doing prelim and final jsd that are not in list
    if(substr($filename,0,5) eq "apid-")
    {
        # get current apid
        $currapid=substr($filename,5,3); 

        # check list of apid want to create JSDs for
        ($foundflg)=&check_apid_list($currapid);

        # if apid not in list skip creating jsd
        if ($foundflg == 0)
        {
          next;#skip doing
        }

        # check if not doing duplicate for merged data series
        ($found_dup_flg)= &check_dup_merged_apid($currapid);        
        if ($found_dup_flg)
        { 
          next;#skip doing
        }
    }
    elsif( substr($filename,0,4) eq "AIA-" || substr($filename,0,4) eq "HMI-" || substr($filename,0,4) eq "SDO-")
    {
        # get current apid
        if(substr($filename,0,7) eq "AIA-ISP")
        {
          $currapid="211"; 
        }
        elsif(substr($filename,0,7) eq "AIA-SEQ")
        {
          $currapid="218"; 
        }
        elsif(substr($filename,0,7) eq "AIA-OBT")
        {
          $currapid="21C"; 
        }
        elsif(substr($filename,0,7) eq "HMI-ISP")
        {
          $currapid="1BD"; 
        }
        elsif(substr($filename,0,7) eq "HMI-SEQ")
        {
          $currapid="1C3"; 
        }
        elsif(substr($filename,0,7) eq "HMI-OBT")
        {
          $currapid="1C0"; 
        }
        elsif(substr($filename,0,7) eq "SDO-ASD")
        {
          $currapid="081"; 
        }
        else
        {
          $currapid=substr($filename,4,3); 
        }


        # check list of apid want to create JSDs for
        ($foundflg)=&check_apid_list($currapid);

        # if apid not in list skip creating jsd
        if ($foundflg == 0)
        {
          next;#skip doing
        }

        # check if not doing duplicate for merged data series
        ($found_dup_flg)= &check_dup_merged_apid($currapid);        
        if ($found_dup_flg)
        { 
          next;#skip doing
        }
    }
    else
    {
        next;
    }

    # passed check if duplicate merged apid and passed being in apids want to do so add to list of apids to do
    if (( index  $filename, "apid-" ) != -1 )
    {
      push(@all_apids, sprintf("%0.4d", hex  substr($filename,5,3)) );
    }
    elsif (( index  $filename, "HMI-" ) != -1 )
    {
      if (substr($filename,0,7) eq "HMI-ISP")
      {
        push(@all_apids, "0445" );
      }
      elsif (substr($filename,0,7) eq "HMI-SEQ")
      {
        push(@all_apids, "0451" );
      }
      elsif (substr($filename,0,7) eq "HMI-OBT")
      {
        push(@all_apids, "0448" );
      }
      else
      {
        push(@all_apids, sprintf("%0.4d", hex  substr($filename,4,3)) );
      }
    }
    elsif (( index  $filename, "AIA-" ) != -1 )
    {
      if (substr($filename,0,7) eq "AIA-ISP")
      {
        push(@all_apids, "0529" );
      }
      elsif (substr($filename,0,7) eq "AIA-SEQ")
      {
        push(@all_apids, "0536" );
      }
      elsif (substr($filename,0,7) eq "AIA-OBT")
      {
        push(@all_apids, "0540" );
      }
      else
      {
        push(@all_apids, sprintf("%0.4d", hex  substr($filename,4,3)) );
      }

    }
    elsif (( index  $filename, "SDO-" ) != -1 )
    {
      if (substr($filename,0,7) eq "SDO-ASD")
      {
        push(@all_apids, "0129" );
      }
      else
      {
        push(@all_apids, sprintf("%0.4d", hex  substr($filename,4,3)) );
      }
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
sub create_file_header($)
{
  #declare and set local variables
  my ($string_apid,$hexapid, $fn_prefix);
  $string_apid= $_[0];#passed argument in function
  

  #Open Output file
  $directory=$ENV{'HK_JSVN_MAP_DIRECTORY'};

  #create outfile file name
  $fn_prefix = &get_filename_prefix($apid);
  $filename= sprintf("%s-%s\n",  $fn_prefix, "JSVN-TO-PVN");

  #create file and start writing lines to hdr file
  open(OUTFILE, ">$directory/$filename") || die "(4)Can't Open  $directory/$filename: $!\n" ;

  # get date to tag creation time of this map file
  &get_date;

  # check if processing SDO-HK apid or HK apid- if HK apid get column to lookup in gtcids file
  if(hex $apid > 96 && hex $apid < 400)
  {
    ;#print "skip calling get_gtcids_lookup_column()\n";
  }
  else
  { 
    ;#print " calling get_gtcids_lookup_column()\n";
    &get_gtcids_lookup_column;
  }
  
  # begin outputing lines for JSD
  $hexapid = sprintf("%s%03X", "0x",$string_apid);
 $filename  =~ s/\n//g;
  print OUTFILE "##############################################################################################################################\n";
  printf(OUTFILE  "# FILENAME:                  %s                                                                             #\n", $filename);

  if($string_apid eq "0445")
  {
  print OUTFILE "# APID(decimal char. value): $string_apid, 475, 29                                                                                   #\n";
  print OUTFILE "# APID(hex. value):          $hexapid, 0x1DB, 0x1D                                                                              #\n";
  print OUTFILE "# PACKET DESCRIPTION:        HMI Image Status Packet                                                                         #\n";
  }
  elsif($string_apid eq "0529")
  {
  print OUTFILE "# APID(decimal char. value): $string_apid, 569, 39                                                                                   #\n";
  print OUTFILE "# APID(hex. value):          $hexapid, 0x239, 0x27                                                                              #\n";
  print OUTFILE "# PACKET DESCRIPTION:        AIA Image Status Packet                                                                         #\n";
  }
  elsif($string_apid eq "0451")
  {
  print OUTFILE "# APID(decimal char. value): $string_apid, 481, 21                                                                                   #\n";
  print OUTFILE "# APID(hex. value):          $hexapid, 0x1E1, 0x15                                                                              #\n";
  print OUTFILE "# PACKET DESCRIPTION:        HMI Sequencer Packet                                                                            #\n";
  }
  elsif($string_apid eq "0536")
  {
  print OUTFILE "# APID(decimal char. value): $string_apid, 576, 46                                                                                   #\n";
  print OUTFILE "# APID(hex. value):          $hexapid, 0x240, 0x2E                                                                              #\n";
  print OUTFILE "# PACKET DESCRIPTION:        AIA Sequencer Packet                                                                            #\n";
  }
  elsif($string_apid eq "0448")
  {
  print OUTFILE "# APID(decimal char. value): $string_apid, 478, 18                                                                                   #\n";
  print OUTFILE "# APID(hex. value):          $hexapid, 0x1DE, 0x12                                                                              #\n";
  print OUTFILE "# PACKET DESCRIPTION:        HMI OBT Packet                                                                                 #\n";
  }
  elsif($string_apid eq "0540")
  {
  print OUTFILE "# APID(decimal char. value): $string_apid, 580, 50                                                                                   #\n";
  print OUTFILE "# APID(hex. value):          $hexapid, 0x244, 0x32                                                                              #\n";
  print OUTFILE "# PACKET DESCRIPTION:        AIA OBT Packet                                                                                 #\n";
  }
  elsif ($string_apid eq "0129")
  {
  print OUTFILE "# APID(decimal char. value): $string_apid                                                                                            #\n";
  print OUTFILE "# APID(hex. value):          $hexapid                                                                                           #\n";
  print OUTFILE "# PACKET DESCRIPTION:        Ancillory Science Packet                                                                        #\n";
  }
  else
  {
  print OUTFILE "# APID(decimal char. value): $string_apid                                                                                            #\n";
  print OUTFILE "# APID(hex. value):          $hexapid                                                                                           #\n";
  print OUTFILE "# PACKET DESCRIPTION:                                                                                                        #\n";
  }

  if(hex $apid < 96 || hex $apid > 400)
  {
    print OUTFILE "# GTCIDS LOOKUP COLUMN:      $columnid                                                                                          #\n";
  }
  print OUTFILE "# DATE LAST UPDATED:         $theTime                                                                      #\n";
  if(hex $apid > 96 && hex $apid < 400)
  {
    print OUTFILE "# DESCRIPTION OF COLUMNS:    The file version numbers are used to lookup JSOC series version number to create series         #\n"; 
    print OUTFILE "#                            name for this SDO HK APID. The Start Packet Date/Time and master file is reference              #\n";
    print OUTFILE "#                            information from the shcids.txt file. The Packet Version is not used for SDO HK APID            #\n";
    print OUTFILE "#                            packets(i.e., 129), but are kept there to be consistent with HK APID map file.                  #\n";
  }
  else
  {
    print OUTFILE "# DESCRIPTION OF COLUMNS:    The JSOC Series and packet version numbers are used to create series name.                      #\n"; 
    print OUTFILE "#                            The File version number,Checkin Date/Time and master file is reference information.             #\n";
  }
  print OUTFILE "##############################################################################################################################\n";

  #check if processing SDO-HK apid(96-399) or HK apid(1-63,400's, or 500's)
  if(hex $apid > 96 && hex $apid < 400)
  {
    print OUTFILE "# JSOC Series    | Packet         | File           | Start Packet Date/Time  | Master File                                   #\n";
  }
  else
  {
    print OUTFILE "# JSOC Series    | Packet         | File           | Checkin Date/Time       | Master File                                   #\n";
  }
  print OUTFILE "# Version Number | Version Number | Version Number |                         |                                               #\n";
  print OUTFILE "##############################################################################################################################\n";
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
  $theTime = sprintf("%-02.2d:%-2.2d:%-2.2d %4s %4s %-02.2d %-4.4d",$hour,$minute,$second, $weekDays[$dayOfWeek], $months[$month], $dayOfMonth, $year);
}
#############################################################################
# subroutine get_gtcids_lookup_column: Find column to lookup packet version #               
# number in Ground to Code IDS file.                                        #
#############################################################################
sub get_gtcids_lookup_column
{
  if($apid == 5 || $apid == 37 || $apid == 448 || $apid == 478 || $apid == 540 || $apid == 580 ||  $apid == 18 ||  $apid == 50 ||  $apid == 580 || $apid == 537  || $apid == 550 || $apid == 567 || $apid == 527 ||  $apid == 435 ||  $apid == 465)
  {
    $columnid ="KERNAL_ID";
    print "WARNING:jsoc_do_jsvn_map_file.pl: using columnid KERNAL_ID in GROUND file for apid:$apid\n";
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
  open(FILE, "$dn_fn") || die "(5)Can't Open file <$dn_fn>: $!\n";
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
         print "WARNING:jsoc_do_jsvn_map_file.pl: Did not find columnid type <$columnid>\n";
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
      elsif ( $jsvn_ht{$packet_version_number} )
      {
      #if jsvn hash table has this packet version number as key then have a JSVN.
        #create line for apid-JSVN-TO-PVN file
        $jsn_map_line = sprintf(" %15.4d | %14s | %14s | %10s %12s | %47s", int $jsvn_ht{$packet_version_number},  $packet_version_number,  $file_version_number, $stha_ci_date, $stha_ci_time, $stha_gen_time);
        print OUTFILE "$jsn_map_line";
      }
      else
      {
        #else just set the JSVN to 0000
        #create line for apid-JSVN-TO-PVN file
        $jsn_map_line = sprintf(" %15.4d | %14s | %14s | %10s %12s | %47s", 0000,  $packet_version_number,  $file_version_number, $stha_ci_date, $stha_ci_time, $stha_gen_time);
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
sub check_apid_list ($)
{
  # apid passed.  Note apidslist is global variable
  $capid = $_[0];
  $foundflg=0;
  foreach $a (@apidslist)
  {
    if($capid =~ m/$a/i)
    {
      $foundflg=1;
      last;
    }
  }
 return ($foundflg);
}

###############################################################################################################
#  get list of apids to do : get list from HKPDF folder and filters out doing HKPDF files with no VER keyword #
###############################################################################################################
sub get_list_apids_to_do($$)
{
  # local variables
  my($argver,$argproctype,$files, @s_fvn);
  # arguments passed
  $argver=$_[0];
  $argproctype=$_[1];

  # misses VER_TEMPERATURE
  # $files=`cd /home1/carl/cvs/TBL_JSOC/lev0/hk_config_file/$file_version/; grep VER_NUM apid*`;
  # get all apids that have VER_NUM or VER_TEMPERATURE Keyword
  @s_fvn=split('\.', $argver );
  if($argproctype eq "HK" &&  int $s_fvn[1] >= 160 &&  int $s_fvn[0] == 1)
  {
    # get hmi files with VER
    $hmifiles =`cd $ENV{'HK_CONFIG_DIRECTORY'}/$argver/;  egrep '(VER_NUM|VER_TEMPERATURE)' HMI-*`;
      
    # remove HMI text regular expression
    $hmifiles =~ s/(H.+?-)//g;
     
    # remove everything from after apid number to end of line regular expression
    $hmifiles =~ s/(-.+?\n)/ /g;
    
    # fix for dealing with apid list -substitute ISP for hmi with 1BD and SEQ for hmi with 1E1
    $hmifiles =~ s/ISP/1BD/g;
    $hmifiles =~ s/SEQ/1C3/g;
   
    # get aia files with VER
    $aiafiles =`cd $ENV{'HK_CONFIG_DIRECTORY'}/$argver/;  egrep '(VER_NUM|VER_TEMPERATURE)' AIA-*`;

    # remove AIA text regular expression
    $aiafiles =~ s/(AIA-)//g;

    # remove everything from after apid number to end of line regular expression
    $aiafiles =~ s/(-.+?\n)/ /g;

    # fix for dealing with apid list -substitute ISP for hmi with 1BD and SEQ for hmi with 1E1
    $aiafiles =~ s/ISP/211/g;
    $aiafiles =~ s/SEQ/218/g;

    #put both aia and hmi lists together
    $files=sprintf("%s%s", $hmifiles,$aiafiles);

  }
  elsif($argproctype eq "HK" &&  int $s_fvn[1] < 160 && int $s_fvn[0] == 1)
  {
    $files=`cd $ENV{'HK_CONFIG_DIRECTORY'}/$argver/;  egrep '(VER_NUM|VER_TEMPERATURE)' apid* `;
    # remove apid text regular expression
    $files =~ s/(a.+?-)//g;
    # remove everything from after apid number to end of line regular expression
    $files =~ s/(-.+?\n)/ /g;
  }
  elsif ($argproctype eq "SDO")
  {
    $sdofiles=`cd $ENV{'HK_SH_CONFIG_DIRECTORY'}/$argver/; grep -H "^APID" SDO*`;

    # remove AIA text regular expression
    $sdofiles =~ s/(SDO-)//g;

    # remove everything from after apid number to end of line regular expression
    $sdofiles =~ s/(-.+?\n)/ /g;

    # fix for dealing with apid list -substitute ISP for hmi with 1BD and SEQ for hmi with 1E1
    $sdofiles =~ s/ASD/081/g;

    #put both aia and hmi lists together
    $files=sprintf("%s", $sdofiles);
  }
  else
  {
    print "ERROR:jsoc_do_jsvn_map_file.pl: Could not understand processing type passed:$argproctype\n";
  }

  # split apids into separate field in array
  #my @apidslist = split(/ /, $files);
  @apidslist = split(/ /, $files);

  # return list of apids to do
  return ( @apidslist);
}

#############################################################################
# sub get_shcids_lines: gets  lines from shcids.txt file                    #
#############################################################################
sub get_shcids_lines
{
  # declare locals
  my($dn,$fn,@all_shcids_lines);

  #open SHCIDS file
  $dn_fn = sprintf("%s/%s",  $ENV{'HK_SH_CONFIG_DIRECTORY'}, $ENV{'HK_SHCIDS_FILE'});
  open(FILE, "$dn_fn") || die "(5)Can't Open: $dn_fn  $!\n";

  #loop throught lines in SHCIDS file
  while (<FILE>)
  {                     
    push(@all_shcids_lines, $_) ;
    @s_line= split / \s*/,  $_;
    if ( substr($s_line[0],0,1) eq '#')
    {
      next; #skip
    }
    else 
    {
      #load key parameters from each line to use in APID-JSVN-PVN map file
      $file_version_number = $s_line[5];
      $gts_ci_date = $s_line[0];
      $gts_ci_time = $s_line[1];
      $gts_name = $s_line[7];
      $gts_ver = $s_line[8];

      # check adminstrator for shcids.txt and SDO_HK_CODE_ids.txt file entered version and filename correctly
      if($gts_ver eq "")
      {
         $gts_ver =~ s/\n//g; #regular exp rm cr 
         $gts_name  =~ s/\n//g; #regular exp rm cr
         print "ERROR:jsoc_do_jsvn_map_file.pl: The version number is not correct in shcids.txt file. Filename and Version value is < $gts_name  $gts_ver>.\n";
         print "Fix this by updating last two columns in file like this: <GODDARD_TLM_SDO.txt   1.1>. Exiting script...\n";
         exit;
      }

      # for SDO-HK apid the packet version number is not applicable so set it to 0.
      $packet_version_number = "000.000";

      # get hashkey based on date
      $hasharraykey =  $gts_ci_date;
      $hasharraykey =~ s/\n|\///g; #regular exp rm cr and slash

      # load string values
      $str_fvn = sprintf("%s",$file_version_number);
      $str_pvn = sprintf("%s",$packet_version_number);

      if ( $jsvn_ht{$hasharraykey} )
      {
        #create line for apid-JSVN-TO-PVN file
        $jsn_map_line = sprintf(" %15.4d | %14s | %14s | %10s %8s | %19s %8.8s", int $jsvn_ht{$hasharraykey},  $packet_version_number,  $file_version_number, $gts_ci_date,$gts_ci_time,$gts_name, $gts_ver);
        print OUTFILE "$jsn_map_line";
      }
      else
      {
        #else just set the JSVN to 0000
        #create line for apid-JSVN-TO-PVN file
        $jsn_map_line = sprintf(" %15.4d | %14s | %14s | %10s %8s | %19s %7s", 0000,  $packet_version_number,  $file_version_number, $gts_ci_date, $gts_ci_time, $gts_name,$gts_ver);
        print OUTFILE "$jsn_map_line";
      }
      $jsn_map_line="";
    }#else
  }#end while
  close(FILE);
  close(OUTFILE);
}

#############################################################################
# subroutine get_filename_prefix:get new prefix name of map filenames       #
#############################################################################
sub get_filename_prefix($)
{
  my($apid_str,$apid);
  $apid_str = $_[0];
  
  #convert to integer apid value
  $apid= int $apid_str;

  #check for merged apids,non-merged apids, or sdo apids and pass back map filename prefix
  if ($apid == 445 )
  {
    return("HMI-ISP");
  }
  elsif ($apid == 529)
  {
    return("AIA-ISP");
  }
  elsif ($apid == 451)
  {
    return("HMI-SEQ");
  }
  elsif ($apid == 536)
  {
    return("AIA-SEQ");
  }
  elsif ($apid == 448)
  {
    return("HMI-OBT");
  }
  elsif ($apid == 540)
  {
    return("AIA-OBT");
  }
  elsif ($apid == 129)
  {
    return("SDO-ASD");
  }
  elsif (($apid < 31 && $apid > 0) || ($apid < 500 && $apid > 400))
  {
     return(sprintf("HMI-%03d",$apid));
  }
  elsif (($apid < 64 && $apid > 30) || ($apid < 600 && $apid > 500))
  {
     return(sprintf("AIA-%03d",$apid));
  }
  else
  {
     print "Got unexpected apid value : $apid\n";
     return ("UNKNOWN");
  }

}
#############################################################################
# subroutine check_dup_merge_apid: check for duplicate of merged apid       #
#############################################################################
sub check_dup_merged_apid($)        
{
  my($apid_str,$apid);
  $apid_str = $_[0];

  #convert to integer apid value
  $apid= hex $apid_str;

  if ($apid == 475  || $apid == 29) #dup HMI ISPs
  {
    return(1);
  }
  elsif ($apid == 481  || $apid == 21) #dup HMI SEQs
  {
    return(1);
  }
  elsif ($apid == 478  || $apid == 18) #dup HMI OBT
  {
    return(1);
  }
  elsif ($apid == 569  || $apid == 39) #dup AIA ISPs
  {
    return(1);
  }
  elsif ($apid == 576  || $apid == 46) #dup AIA ISPs
  {
    return(1);
  }
  elsif ($apid == 580  || $apid == 50) #dup AIA OBTs
  {
    return(1);
  }
  else
  {
    return(0);
  }

} 
