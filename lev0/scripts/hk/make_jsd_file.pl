#!/usr/local/bin/perl
##############################################################################
# Name:        make_jsd_file.pl - Make JSD File                              #
# Description: Strip data from HKPDF files and create global and keyword     #
#              section for the JSOC Series Definition file. Uses GTCIDS file #
#              to get packet version number data. Use the <APID>-JSVN-TO-PVN #
#              map files to get jsoc series version number. Run first to make#
#              prelim jsd files. Then run a second time to make final JSDs.  #
#              The prelim files need to be there for all packet versions to  #
#              successfully create final jsd file.                           #
# Process:     (1)Run first to make prelim jsd files for needed APIDs.       #
#              (2)Run do_jsvn_map_file.pl to make JSVN-TO-PVN map files for  #
#                 each APID.                                                 #
#              (3)Run make final jsd file to make the final files to use to  #
#                 create.                                                    #
# Execution:   (1)To run :                                                   #
#              make_jsd_files.pl <prelim | final > <File-Version-Number>     #
#              (2)For help: make_jsd_file.pl  -h                             #
# Example:     make_jsd_file.pl prelim 1.79                                  #
# Example:     make_jsd_file.pl final  1.79                                  #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on March 21, 2008 #
##############################################################################
# main function                                                              #    
##############################################################################
##get environment variables and initialize variables.
$hm=$ENV{'HOME'};
$ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
$ENV{'HK_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_file";
$ENV{'HK_JSD_PRELIM_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_prelim_file";
$ENV{'HK_JSVN_MAP_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsn_map_file";
$ENV{'HK_GTCIDS_FILE'}="gtcids.txt";
$lnjsds=$ENV{'HK_LIST_NEW_JSDS'}="$hm/cvs/JSOC/proj/lev0/scripts/hk/new_jsd_files.txt";

if ($#ARGV < 0 || $#ARGV > 1) 
{
  die "Usage: $0 < prelim | final > <Packet Version Number | -h > :$!";
}

if ( $ARGV[0] eq "-h" )
{
  print "Help Listing\n";
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Create Preliminary JSD Files: make_jsd_file.pl prelim <File Version Number>\n";
  print "(1a)Create Final JSD Files:       make_jsd_file.pl final <File Version Number>\n";
  print "(1b)Get Help Information:          make_jsd_file.pl -h \n";
  print "(2)Environment variable HK_JSD_PRELIM_DIRECTORY and HK_JSD_DIRECTORY\n"; 
  print "(2a)Used to control where to put output final JSD file\n";
  print "(3) Environment variable HK_GTCIDS_FILE\n";
  print "(3a) Used to store filename of gtcids file which is use to look up packet version number.\n";
  print "(4) Environment variable HK_CONFIG_DIRECTORY\n";
  print "(4a) Used to store location of hk config (HKPDF) files for each apid\n";
  print "(4b) The filename format is apid-<apid#>-version-<version#>.\n";
  print "(4c) The HKPDF files are used to get print values, keyword and long names, etc >.\n";
  exit;
}
elsif ( $ARGV[0] eq "prelim")
{
  $file_version_number = $ARGV[1];
  $ffp= 0; #flag for final or prelim value;
}
elsif ( $ARGV[0] eq "final")
{
  $file_version_number = $ARGV[1];
  $ffp= 1;
}
else
{
  print "Error: use either prelim or final for first argument(i.e.,prelim)\n";
  print "Error: use packet version number for second argument(i.e.,1.132)\n";
}

# (1) get all file in folder based on file version number
&get_all_HKPDF_files();

##open New JSD File for holding new jsd files created
if($ffp == 1) 
{
  open(NJF, ">$lnjsds") || die "(1)Can't Open  $lnjsds: $!\n" ;
}

# (2) loop thru all HKPDF file and create a Final JSD file for each apid in filename
foreach $file ( @fn_list)
{
   if ( $file eq "." || $file eq ".." || $file eq "CVS")
   {
      next;#skip doing
   }

# (3) get all lines in HKPDF file
   &get_all_lines();

#(4) get packet version number based on file version number used
   &get_packet_version_number();
   if ($found_fvn == 0)
   {
     # This occurs when there is a STHA file version but is not mapped in GTCIDS file
     # to a packet version number. This are no file versions 1.53,1.66,1.68,1.70.
     print "NOTE:For Apid <$apid_value>:Skipping creating preliminary jsd file.There is no packet version number ";
     print "found corresponding to file version number<$file_version_number> in GTCIDS file.\n";
     next;
   }
#(5) open <APID>-JSVN-TO-PVN file and locate packet-version# and jsvn#
   if($ffp == 1) 
   {
     &get_jsvn(); 
   }
#(6) create jsd file based on data in each HKPDF file
   &create_jsd_file(); 
}

#close list in New JSD File 
if ($ffp == 1)
{
  close NJF;
}

###############################################################################
# sub get_jsvn: gets jsoc series version number 
###############################################################################
sub get_jsvn()
{
# setup my variables
  my($dn,$fn,$mfn_list,$mapfile,$found_flag);
# get directory name to find jsvn map files
  $dn = sprintf("%s",  $ENV{'HK_JSVN_MAP_DIRECTORY'});
  opendir(DIR_JSVN_MAP, $dn) || die ":(2)Can't open:$!\n"; #open subdirectory
# read in all map files and put in list 
  @mfn_list = readdir(DIR_JSVN_MAP); #get a list of directory contents
  $found_flag = "0";
  foreach $mapfile  (@mfn_list)
  {
    if(substr( int  $mapfile,0,4) ne  hex substr($apid_value,2,4))
    {
       next;
    }
    else
    {
      #found map file for this APID
      $found_flag = "1";
      $dn_fn = sprintf("%s/%s",  $ENV{'HK_JSVN_MAP_DIRECTORY'}, $mapfile);
      #open map file
      open(FILE1, "$dn_fn") || die "(3)Can't Open: $!\n";
      #loop thru map file looking for file version number in 3nd column or 5th field 
      #and compare with file version number passed in as argument 2.
      @all_map_file_lines="";
      while (<FILE1>)
      {
        push(@all_map_file_lines, $_) ;
        @s_line= split / \s*/,  $_;
        if ( substr($s_line[0],0,1) eq '#')
        {
          next; #skip
        }
        elsif( $s_line[5] eq $file_version_number)
        {
          $jsvn= $s_line[1];
        }
        else 
        {
         #print "skipping\n";
         #print "fvn  is  $s_line[5]\n";
        }
      }#while
      close FILE1;
      closedir DIR_JSVN_MAP; #close directory
      return;
    }#else
  }#foreach
  closedir DIR1; #close directory
  if ($found_flag eq "0")
  {
    $jsvn = "999";#no JSVN Map files exists
  }
}
  
###############################################################################
# sub get_file_version_number: gets file version number based on packet version#
###############################################################################
sub get_packet_version_number
{
  my($dn,$fn,@all_gtcids_lines);
# open GTCIDS file
  $dn_fn = sprintf("%s/%s",  $ENV{'HK_CONFIG_DIRECTORY'}, $ENV{'HK_GTCIDS_FILE'});
  open(FILE, "$dn_fn") || die "(4)Can't Open: $!\n";
  $found_fvn = 0;
  while (<FILE>)
  {
#     push(@all_gtcids_lines, $_) ;
      @s_line= split / \s*/,  $_;
      if ( substr($s_line[0],0,1) eq '#')
      {
        next; #skip
      }
      else
      {
        if( $sys_value eq "hmi" ) 
        {
          if ( $s_line[13] eq $file_version_number)
          {
            $packet_version_number = $s_line[7];
            $found_fvn = 1;
            $stha_ci_date = $s_line[0];
            $stha_ci_time = $s_line[1];
            $stha_gen_time = $s_line[15];
            $found_fvn = 1;
            break;
          }
        }
        elsif( $sys_value eq "aia" ) 
        {
          if ( $s_line[13] eq $file_version_number)
          {
            $packet_version_number = $s_line[9];
            $found_fvn = 1;
            $stha_ci_date = $s_line[0];
            $stha_ci_time = $s_line[1];
            $stha_gen_time = $s_line[15];
            $found_fvn = 1;
            break;
          }
        }
        elsif( $sys_value eq "ssim" ) 
        {
          if ( $s_line[13] eq $file_version_number)
          {
            $packet_version_number = $s_line[7]; #use hmi column as default value
            $found_fvn = 1;
            $stha_ci_date = $s_line[0];
            $stha_ci_time = $s_line[1];
            $stha_gen_time = $s_line[15];
            $found_fvn = 1;
            break;
          }
        }
        else
        {
           print "Warning: $sys_value is not hmi,aia, or ssim. Check setting.\n";
           $found_fvn = 0;
        }
      }#big else
  }
  close(FILE);
}
###############################################################################
# sub get_all_hkpdf_files: gets all files in folder specified by file version#
###############################################################################
sub get_all_HKPDF_files
{
# make directory name to find hkpdf files
  $dn = sprintf("%s/%s",  $ENV{'HK_CONFIG_DIRECTORY'}, $file_version_number);
  opendir(DIR, $dn) || die "(5)Can't open:$!\n"; #open subdirectory
  @fn_list = readdir(DIR); #get a list of directory contents
  closedir(DIR); #close directory
} 
###############################################################################
# sub get_all_lines: gets all lines from HKPDF file
###############################################################################
sub get_all_lines
{
    
    my($i);
    @all_kwd_lines="";
    @all_dcon_lines="";
    @all_acon_lines="";
#   skip over processing . or CVS, etc files
    if(substr($file,0,4) ne "apid" && substr($file,9,7) ne "version" )
    {
#       print "skipping reading in lines for  file <$file>\n";
        return;
    }
#   open one apid-#-version-# file
    open(FILE, "$dn/$file") || die "(6)Can't Open $file file: $!\n";
    while (<FILE>)
    {
      push(@all_lines, $_) ;
       @s_line= split / \s*/,  $_;
      if ( $s_line[0] eq "FILE")
      {
        $file_line = $_ ;
        $file_version= $s_line[2] ;
        next;
      }
      if ( $s_line[0] eq "APID" )
      {
        $apid_lines= $_ ;
        $apid_value = $s_line[1] ;
        $sys_value = lc $s_line[3] ;
        $i=0;
        @desc="";
        while( $s_line[4 + $i] ne "") 
        {
           push( @desc, $s_line[4 + $i]);
           $i++;
           next;
        }
        pop(@desc); #pop off last element-timestamp
        next;
      }
      if ( $s_line[0] eq "KWD" )
      {
        push(@all_kwd_lines, $_) ;
        next;
      }
      if ( $s_line[0] eq "DCON" )
      {
        push(@all_dcon_lines, $_) ;
        next;
      }
      if ( $s_line[0] eq "ACON" )
      {
        push(@all_acon_lines, $_) ;
        next;
      }
    }#end-while
    close( FILE);
    $all_kwd_count=@all_kwd_lines;
    $all_acon_count=@all_acon_lines;
    $all_dcon_count=@all_dcon_lines;
}
###############################################################################
# sub create_jsd: create jsd based on data in HKPDF file
###############################################################################
sub create_jsd_file
{
  if ($ffp == 1 )
  {
    if ($jsvn eq "999")
    {
#       print "skip doing jsn is 999\n";
       return;
    }
  }
  my($default_value);
# use environment to set or use arguments(-n<namespace> -i<identifer> ) to set
  $jsd_namespace = sprintf("%s%s", $sys_value,"_ground");
  $jsd_identifier = "lev0";
  $author="production";
  $owner="production";
  $unitsize="1";
  $archive= "1";
  $retention= "30";
  $index="PACKET_TIME";

# set tape ground based on if hmi,aia or ssim 
  if ( $sys_value eq "hmi")
  {
    $tapegroup= "2";
    $retention= "60";
  }
  elsif ( $sys_value eq "aia")
  {
    $tapegroup= "3";
    $retention= "60";
  }
  elsif ( $sys_value eq "ssim")
  {
    $tapegroup= "1";
    $retention= "60";
  }
  elsif ( $sys_value eq "unkn")
  {
    $tapegroup= "1";
  }
  else
  {
    print "Warning: Do not know sys value <$sys_value>\n";
  }

# format packet version number in MMMmmm decimal format 
  $key=".";
  $pos= index  $packet_version_number ,  $key;
  if ( $packet_version_number < 10)
  {
    $pvn_wn= substr($packet_version_number,0,1);
    $pvn_dn= substr($packet_version_number, 2);
  }
  elsif(  $packet_version_number > 10 && $packet_version_number < 100) 
  {
    $pvn_wn= substr($packet_version_number,0,2);
    $pvn_dn= substr($packet_version_number, 3);
  }
  elsif (  $packet_version_number > 100 && $packet_version_number < 1000 )
  {
    $pvn_wn= substr($packet_version_number,0,3);
    $pvn_dn= substr($packet_version_number, 4);
  }
  else
  {
     print "**WARNING: The Packet Version number for jsd file name is not set correctly.\n";
     print "Check code and updating perl script is probably needed.\n";
  }

  if ($ffp == 1)
  {
    # create file name for final jsd files
    $jsd_fn= sprintf("%s.%s_%-04.4d_%-04.4d",$jsd_namespace,$jsd_identifier,hex  substr($apid_value,2,4), $jsvn);
    # get directory to write jsd file
    $directory=$ENV{'HK_JSD_DIRECTORY'};
  }
  else ## ffp==0
  {
    # create file name for preliminary jsd files used for creating jsn map files.
    $jsd_fn= sprintf("%s.%s_%-04.4d_%-03.3d%-03.3d",$jsd_namespace,$jsd_identifier,hex  substr($apid_value,2,4), $pvn_wn,$pvn_dn);
    # get directory to write jsd file
    $directory=$ENV{'HK_JSD_PRELIM_DIRECTORY'};
  }

# check if version_number folder exists
  opendir(DIR,$directory) || die "(7)Can't open: $!\n";

# check if file exists
  @items=readdir( DIR );
  closedir(DIR);
  $found_file=0;
  foreach $current_file ( @items )
  { 
    if ($current_file eq $jsd_fn)
    {
      $found_file=1;
      last;
    }
  }
  if ( $found_file eq 1 )
  {
     if($ffp == 0)    
    {
      print "-->updating existing preliminary jsd file <$jsd_fn>.\n";
    }
    else
    {
      print "-->updating existing final jsd file <$jsd_fn>.\n";
    }
  }
  else
  {
     if($ffp == 0)    
    {
      print "-->created new preliminary jsd file <$jsd_fn>.\n";
    }
    else
    {
      print "-->create new final jsd file <$jsd_fn>.\n";
      print NJF "$jsd_fn\n";
    }
  }


# Create directory using environment variable, version number and filename
  open(OUTFILE, ">$directory/$jsd_fn") || die "(8)Can't Open  $directory/$jsd_fn: $!\n" ;
# Create global section of file
  if ($ffp == 1)
  {
    printf(OUTFILE "#======================== Global Series Information =========================\n");
    printf(OUTFILE  "SeriesName:  %-s\n", $jsd_fn);
    print OUTFILE   "Description:@desc \n" ;
    printf(OUTFILE  "Author:      \"%-s\"\n", $author);
    printf(OUTFILE  "Owner:       \"%-s\"\n", $owner);
    printf(OUTFILE  "Unitsize:    %-s\n", $unitsize);
    printf(OUTFILE  "Archive:     %-s\n", $archive);
    printf(OUTFILE  "Retention:   %-s\n", $retention);
    printf(OUTFILE  "Tapegroup:   %-s\n", $tapegroup);
    printf(OUTFILE  "Index:       %-s\n", $index);
  }

# Create keyword section of file
  printf(OUTFILE "#======================= Keywords Series Information ========================\n");

# create and print to out file index keyword
  $keyword_name = "PACKET_TIME";
  $type_value = "time";
  $print_value="UTC";
  $default_value = "1977.01.01";

  $unit_value="\"UTC\"";
  $comment_value = "\"PACKET_TIME, TIMECODE time in UTC Format\"";

# Note: Set inialially default value to TBD. Then reset and
# Generate UTC format YYYY.MM.DD_HH:MM:SS.ss from TAI time
  $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
              "Keyword:", $keyword_name, $type_value,
              "variable", "record", $default_value, 
              $print_value, $unit_value, $comment_value); 
  printf(OUTFILE "%s", $tmp_line);

# create and print to out file PACKET VERSION NUMBER keyword
  $keyword_name = "PACKET_VERSION_NUMBER";
  $type_value = "string";
  $print_value="%s";
  $default_value ="0";
  $unit_value="\"none\"";
  $comment_value = "\"PACKET_VERSION_NUMBER, decimal value with format:MMM.mmm\"";
  $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
              "Keyword:", $keyword_name, $type_value,
              "variable", "record", $default_value, 
              $print_value, $unit_value, $comment_value); 
  printf(OUTFILE "%s", $tmp_line);

# print out line for each keyword
  foreach $kwd_line (@all_kwd_lines)
  {
    if($kwd_line eq "") 
    { 
      next; #skip lines with no data
    }

# get keyword name
  @s_line = split / \s*/, $kwd_line;
  $keyword_name = $s_line[1];

# get variable type and  printing format to use
  if ($s_line[6] eq "UB" && $s_line[7] eq "R")
  {
    $type_value = "short";
    $print_value=" %hd";
    $default_value ="0";
  }
  elsif($s_line[6] eq "SB" && $s_line[7] eq "R") 
  {
    $type_value = "char";
    $print_value=" %hhd";
    $default_value ="0";
  }
  elsif($s_line[6] eq "IU1" && $s_line[7] eq "R") 
  {
    $type_value = "int";
    $print_value=" %d";
    $default_value ="0";
  }
  elsif($s_line[6] eq "IS1" && $s_line[7] eq "R") 
  {
    $type_value = "short";
    $print_value=" %hd";
    $default_value ="0";
  }
  elsif($s_line[6] eq "IL1" && $s_line[7] eq "R") 
  {
    $type_value = "int";
    $print_value=" %d";
    $default_value ="0";
  }
  elsif($s_line[6] eq "UL1" && $s_line[7] eq "R") 
  {
    $type_value = "longlong";
    $print_value=" %ld";
    $default_value ="0";
  }
  elsif($s_line[7] eq "D") 
  {
    $type_value = "string";
    $print_value=" %s";
    $default_value ="0";
  }
  elsif($s_line[7] eq "A") 
  {
    $type_value = "double";
    $print_value=" %lf";
    $default_value ="0.0";
  }
  else
  {
    $type_value = "int";
    $print_value=" %d";
    $default_value ="0";
    print "**WARNING: Apid: $apid_value Keyword: $s_line[1] type value not found. NOT VALID TYPE! Skipping using keyword.\n";
    next;#skip writing to jsd file
  }

# use hardcoded value for units
  $unit_value="\"none\"";

# use long name in comments
  $comment_value = sprintf("\"%s\"", $keyword_name);

# create line to print to jsd file
  $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
              "Keyword:", $s_line[2], $type_value,
              "variable", "record", $default_value, 
              $print_value, $unit_value, $comment_value); 
  printf(OUTFILE "%s", $tmp_line);

  }#end-for-each

# Create Link section in jsd file
  printf(OUTFILE "#========================= Link Series Information ========================\n");
  printf(OUTFILE "#Not applicable\n");

# Create Segment section in jsd file
  printf(OUTFILE "#======================== Segment Series Information ========================\n");
  printf(OUTFILE "#Not applicable\n");

# Close output file
  close( OUTFILE);

}

