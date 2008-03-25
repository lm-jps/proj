#!/usr/bin/perl
#############################################################################
#NAME: make_hkdpf.pl
#Author: Carl
#Purpose:Strip data from STANFORD_TLM_HMI_AIA.txt file and create HKPDF files
#        using the filename format apid-<#>-version-<#>. These HKPDF files 
#        contain the required data for processing housekeeping data.
#Execution for help: make_hkpdf.pl -h 
#Execution to create hkpdf files: make_hkpdf.pl 4 
#Create Date:3/15/2008 
#############################################################################
# main function
#############################################################################
#get environment variables and initialize variables.
$hm=$ENV{'HOME'};
$ENV{'HK_APID_LIST'} = "ALL";
$ENV{'HK_STHA_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/fromlmsal";
$ENV{'HK_STHA_FILENAME'}="STANFORD_TLM_HMI_AIA.txt";
$ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file/";
$ENV{'HK_GTCIDS_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/fromlmsal";
$ENV{'HK_GTCIDS_FILE'}="gtcids.txt";

#common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";
$logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-clmuq";

#display environment variables
#print "--->make_hkpdf enviroment variables:\n";
#print `env `;

if ( $ARGV[0] eq "-h" )
{
  print "Help Listing\n";
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Create Files: make_hkdpf.pl \n";
  print "(1b)Create Files: make_hkdpf.pl <sort-choice= 1, 2, 3, 4>\n";
  print "where sort-choice equal to 0 or no argument is no sort \n";
  print "where sort-choice equal to 1 is sort by keyword\n";
  print "where sort-choice equal to 2 is sort by byte position\n";
  print "where sort-choice equal to 3 is sort by long telem name\n";
  print "where sort-choice equal to 4 is sort by byte and bit position\n";
  print "************************************\n";
  print "** Note:Currently using option 4! **\n";
  print "************************************\n";
  print "(1c)Get Help Information: make_hkdpf.pl -h \n";
  print "(2)Environment variable HK_APID_LIST\n"; 
  print "(2a)Used to control which apid to process\n";
  print "(2b)If set to ALL all apid will be processed\n";
  print "(2c)If set to 0x1BD 0x1DB then only 2 apid will be processed\n";
  print "(3) Environment variable HK_STHA_DIRECTORY\n";
  print "(3a) Used to store location of STANFORD_TLM_HMI_AIA.txt file.\n";
  print "(4) Environment variable HK_STHA_FILE\n";
  print "(4a) Used to store filename of STANFORD_TLM_HMI_AIA.txt file to process.\n";
  print "(4b) This file is input data for this script\n";
  print "(5) Environment variable HK_CONFIG_DIRECTORY\n";
  print "(5a) Used to store location of hk config files for each apid\n";
  print "(5b) These files are the output for this script.\n";
  print "(5c) The filename format is apid-<apid#>-version-<version#>.\n";
  exit;
}
# Set sort option to default and check for argument passed
$sort_option="0";
if ( $ARGV[0] eq "0" )
{
  $sort_option="0";
}
elsif ( $ARGV[0] eq "1" )
{
  $sort_option="1";
}
elsif ( $ARGV[0] eq "2" )
{
  $sort_option= "2";
}
elsif ( $ARGV[0] eq "3" )
{
  $sort_option="3";
}
elsif ( $ARGV[0] eq "4" )
{
  $sort_option="4";
}
elsif ( $ARGV[0] eq "" )
{
  $sort_option="0";
}
else
{
  print "**make_hk.pl:Error bad argument: choice 0, 1, 2, 3, 4 or no argument \n";
  print "where sort-choice equal to 0 or no argument is no sort \n";
  print "where sort-choice equal to 1 is sort by keyword\n";
  print "where sort-choice equal to 2 is sort by byte position\n";
  print "where sort-choice equal to 3 is sort by long telem name\n";
  exit;
}
# (1) get all lines in STHA file
&get_all_lines();

# (2) get all apids to process based on environment variable 
&get_apids_to_do;

#print "Please wait creating  < $all_apids_count > hkdpf files\n";
# (3) create HKPDF filenames
&create_hkpdf_filenames;

# (4) create HKPDF files
&create_hkpdf_files;

# (5) print stats:total of apid files created
#print ".. completed creating < $all_apids_count > hkdpf files\n";

# (6) create workaround file apid-1AF-version-<current version number>
#&create_hkpdf_for_1AF; not used anymore

#LOG
print "--->Finished running make_hkpdf \n";


###############################################################################
# sub get_all_lines: gets all lines in STANFORD_TLM_HMI_AIA.txt file
###############################################################################
sub get_all_lines
{
  $dn=$ENV{'HK_STHA_DIRECTORY'};
#print "$dn\n";
  $fn=$ENV{'HK_STHA_FILENAME'};
#print "$fn\n";
  open(FILE, "$dn/$fn") || die "Can't Open $fn file: $!\n";
  while (<FILE>)
  {
    push(@all_lines, $_) ;
    @s_line=split(',',$_ );
    if ( $s_line[0] eq "PKT")
    {
      push(@all_pkt_lines, $_) ;
      push(@all_apids, $s_line[1]) ;
      next;
    }
    if ( $s_line[0] eq "FITS" )
    {
      push(@all_fits_lines, $_) ;
      next;
    }
    if ( $s_line[0] eq "TLM" )
    {
      push(@all_tlm_lines, $_) ;

      @mnm=split('_', $s_line[1]);
      if ( $mnm[0] eq "T" )
      {
        push(@all_tlm_t, $_);
      }
      if ( $mnm[0] eq "E" )
      {
        push(@all_tlm_e, $_);
      }
      next;
    }
    if ( $s_line[0] eq "ALG" )
    {
      push(@all_alg_lines, $_) ;
      next;
    }
    if ( $s_line[0] eq "DSC" )
    {
      push(@all_dsc_lines, $_) ;
      next;
    }
  }
  close( FILE);
  $line_count=@all_lines;
  $pkt_line_count=@all_pkt_lines;
  $all_apids_count=@all_apids;
  $all_fits_count=@all_fits_lines;
  $all_tlm_count=@all_tlm_lines;
  $all_tlm_tcount=@all_tlm_t;
  $all_tlm_ecount=@all_tlm_e;
  $all_alg_count=@all_alg_lines;
  $all_dsc_count=@all_dsc_lines;
##print "LINES=$line_count\nPKTS=$pkt_line_count\nAPIDS=$all_apids_count\n";
#print "FITS=$all_fits_count\nTLMS=$all_tlm_count\nALGS=$all_alg_count\nDSCS=$all_dsc_count\n";
#print "TLMS_T=$all_tlm_tcount\nTLMS_E=$all_tlm_ecount\n";
}

###############################################################################
# sub get_all_pkt_lines: gets all pkt lines and all apids in STHA file
###############################################################################
sub get_all_pkt_lines  
{
  my($i);
  $i=0;
  while ( $all_lines[$i] )
  {
    @apid_line=split(',',$all_lines[$i] );
    if ( $apid_line[0] eq "PKT" )
    {
      push(@all_pkt_lines, $all_lines[$i]) ;
      push(@all_apids, $apid_line[1]) ;
    }
    $i++;
  }
  $pkt_line_count=@all_pkt_lines;
  $all_apids_count=@all_apids;
}

###############################################################################
# sub check_apid_to_do: check which apid to do, if ALL do all in STHA file
###############################################################################
sub get_apids_to_do
{
  $arg=$ENV{'HK_APID_LIST'};
  if ( $arg eq "ALL" || $arg eq "")
  {
     $apids_flag=ALL;
  }
  else 
  {
     $apid_flag=SOME;
  }

  if ( $apids_flag ne "ALL")
  {
    while ( @all_apids )
    {
      shift @all_apids;
    }
    @arg_items=split(' ',$arg );
    while (@arg_items)
    {
      push(@all_apids,$arg_items[0] );
      shift @arg_items;
    }
  }
$all_apids_count=@all_apids;
print "APID List to process is as follows: @all_apids \n";
}

###############################################################################
# sub create_hkpdf_filenames: create filenames for each apid.
###############################################################################
sub create_hkpdf_filenames
{
  my($i, $j);
##get version number
  $i=0;
  while ($all_lines[$i])
  {
    @file_line=split(' ', $all_lines[$i]);
    if ($file_line[1] eq "File" && $file_line[2] eq "RCS")
    {
      $version_number=$file_line[6];
#print "VERSION NUMBER is: $version_number \n";
      last;
    }
    $i++;
  }

##get apid number
  $i=0;
  while ( $all_apids[$i] )
  {
    @apid_line=split(',',$all_pkt_lines[$j] );
    if ( hex $apid_line[1] ne  hex $all_apids[$i] )
    {
      $j++;
    }
    else
    {
      @apid_line_hex=split('x',$apid_line[1]);
      push(@filenames, sprintf("apid-%3.3s-version-%s\n",  $apid_line_hex[1], $version_number));
      $i++;
    }
  }
##print "filename list:\n  @filenames \n";
}
###############################################################################
# sub create_hkpdf_files: create files for each apid to do.
###############################################################################
sub create_hkpdf_files
{
  my($i);
  $directory=$ENV{'HK_CONFIG_DIRECTORY'};
# Check if version_number folder exists
  opendir(DIR,"$directory") || die "Can't open: $!\n";
  @items=readdir(DIR);
  closedir(DIR);
  $found_dir=0;
  foreach $folder ( @items )
  { 
    if ($folder eq $version_number)
    {
      $found_dir=1;
      last;
    }
  }
  if ( $found_dir ne 1 )
  {
      mkdir( "$directory$version_number", 0755) || die "mkdir:$version_number: $!\n";
  }
  print "Writing new files to directory: $directory$version_number \n";
  $i=0;
  while ($all_apids[$i])
  {
# Create directory using environment variable, version number and filename
    open(OUTFILE, ">$directory$version_number/$filenames[$i]") || die "Can't Open  $directory/$filenames[$i]: $!\n" ;
    &create_file_line;
    &create_apid_line( $all_apids[$i]) ;
    &create_kwd_lines( $all_apids[$i]) ;
    &create_acon_lines() ;
    &create_dcon_lines() ;
    close( OUTFILE);
    $i++;
    print "Finish file $i\n";
  }
}

###############################################################################
# sub create_pkt_line: create FILE line for hkpdf file.
###############################################################################
sub create_file_line
{
  my($i,$f_line,$fileline);
  $i=0;
  while ($all_lines[$i])
  {
    @f_line=split(' ', $all_lines[$i]);
    if ($f_line[1] eq "File" && $f_line[2] eq "RCS")
    {
      last;
    }
    $i++;
  }        
  $fileline=sprintf("%4.4s %24.24s  %-7s %-10s  %-8s\n","FILE",$f_line[5],$f_line[6], $f_line[7], $f_line[8]); 
  print OUTFILE "# hkdpf file created by make_hkdpf.pl script \n";
  print OUTFILE "$fileline";
}
###############################################################################
# sub create_apid_line: create apid line for hkpdf file.
###############################################################################
sub create_apid_line
{
  my $find_apid =shift ;
  my($i, $p_line, $apidline);
  $i=0;
  while ( $all_pkt_lines[$i] )
  {
    @p_line=split(',',$all_pkt_lines[$i] );
    if ( $p_line[1] eq $find_apid )
    {
      last;
    }
    $i++;
  }
  if ( hex $p_line[1] >= 1 && hex $p_line[1] <= 31)
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "HMI", $p_line[5], $p_line[6]); 
  }
  elsif ( hex $p_line[1] >= 400 && hex $p_line[1] <= 486)
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "HMI", $p_line[5], $p_line[6]); 
  }
  elsif ( hex $p_line[1] > 31 && hex $p_line[1] <= 63)
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "AIA", $p_line[5], $p_line[6]); 
  }
  elsif ( hex $p_line[1] >= 500 && hex $p_line[1] <= 580)
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "AIA", $p_line[5], $p_line[6]); 
  }
  elsif ( hex $p_line[1] >= 2002 && hex $p_line[1] <= 2047)
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "SSIM", $p_line[5], $p_line[6]); 
  }
  else
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "UNKN", $p_line[5], $p_line[6]); 
  }


##$apidline=sprintf("%4.4s %5.5s  %3.3s  %8.8s\n","APID",$p_line[1],$p_line[3], $p_line[6]); 
  print OUTFILE "$apidline";
}
###############################################################################
# sub create_kwd_line: create keyword line for hkpdf file. 
###############################################################################
sub create_kwd_lines
{
  
  my $find_apid =shift ;
  my($i, $p_line, $kwdline, $mnm);
  $i=0;
  while (@apid_kwd_lines)
  {
    shift @apid_kwd_lines;
  }
  $i=0;
  while ( $all_lines[$i] )
  {
    @line=split(',',$all_lines[$i] );
    if ( $line[0] eq "PKT" )
    {
      
      if ( $line[1] eq $find_apid )
      {
        $i++;
        $j=0;
        while ($all_lines[$i])
        {
          @line=split(',',$all_lines[$i] );
          if ( $line[0] eq "SRC" )
          {
            @mnm=split('_', $line[1]);


            if ( $mnm[0] ne "T" && $mnm[0] ne "E" && $mnm[0] ne "Y") 
            {
              &get_conversion_flag ($line[1]);
              &get_keyword_value ($line[1]);
              $kwdline[$j]=sprintf("%-3s %-8s  %-40s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s","KWD",$kwd_value,$line[1],$line[4],$line[5],$line[6],$line[7],$conv_value,$line[10]); 
              push(@apid_kwd_lines, $kwdline[$j]);
#              print OUTFILE "$kwdline[$j++]";
            }
          }
          else 
          {
            last;
          }
          $i++;
        }
        last;
      }
    }
    $i++;
  }
################################
# SORT                         #
################################
  &sort_kwd_lines ;
   
} 

###############################################################################
# sub sort_kwd_lines: sort keyword lines by byte-position, long and short names
###############################################################################
sub sort_kwd_lines
{
  my ($j, @s_line, @tmp_lines, @s_oline);
#################################
# No Sort                       #
#################################
  if ($sort_option eq "0")
  {
    print "Sorting keyword lines by first in first out. ";
    $j=0;
    while( $apid_kwd_lines[$j] )
    {
      @s_line = split / \s*/,$apid_kwd_lines[$j];
      printf(OUTFILE "%-3.3s  %-8.8s  %-35.35s  %-3.3s  %-1.1s  %-2.2s  %-3.3s %-1.1s  %-8.8s \n",$s_line[0],$s_line[1],$s_line[2],$s_line[3],$s_line[4],$s_line[5],$s_line[6], $s_line[7],$s_line[8]);
      $j++;
    }
  }
#################################
# Sort by Short Keyword Names   #
#################################
  if( $sort_option eq "1")
  {
    print "Sorting keyword lines by short keyword name. ";
    $j=0;
    while ( $apid_kwd_lines[$j] )
    {
      @s_line = split / \s*/,$apid_kwd_lines[$j];
      $tmp_lines[$j]=sprintf("%-8.8s  %-3.3s  %-8.8s  %-35.35s  %-3.3s  %-1.1s  %-2.2s  %-3.3s %-1.1s  %-8.8s \n",$s_line[1],$s_line[0],$s_line[1],$s_line[2],$s_line[3],$s_line[4],$s_line[5],$s_line[6], $s_line[7],$s_line[8]); 
      $j++;
    }
  }
#################################
# Sort by byte position         #
#################################
  if( $sort_option eq "2")
  {
    print "Sorting keyword lines by byte positions. ";
    $j=0;
    while ( $apid_kwd_lines[$j] )
    {
      @s_line = split / \s*/,$apid_kwd_lines[$j];
      $tmp_lines[$j]=sprintf("%-3.3d  %-3.3s  %-8.8s  %-35.35s  %-3.3s  %-1.1s  %-2.2s  %-3.3s %-1.1s  %-8.8s \n",$s_line[3],$s_line[0],$s_line[1],$s_line[2],$s_line[3],$s_line[4],$s_line[5],$s_line[6], $s_line[7],$s_line[8]); 
      $j++;
    }
  }
###################################
# Sort by Long Telem Mnemonic Name# 
###################################
  if( $sort_option eq "3")
  {
    print "Sorting keyword lines by long telemetry names. ";
    $j=0;
    while($apid_kwd_lines[$j])
    {
      @s_line = split / \s*/,$apid_kwd_lines[$j];
      $tmp_lines[$j]=sprintf("%-35.35s  %-3.3s  %-8.8s  %-35.35s  %-3.3s  %-1.1s  %-2.2s  %-3.3s %-1.1s  %-8.8s \n",$s_line[2],$s_line[0],$s_line[1],$s_line[2],$s_line[3],$s_line[4],$s_line[5],$s_line[6], $s_line[7],$s_line[8]); 
      $j++;
    }
  }
#########################################
# Sort by Byte Position and Bit Position# 
#########################################
  if( $sort_option eq "4")
  {
    print "Sorting keyword lines by byte and bit positions. ";
    $j=0;
    while($apid_kwd_lines[$j])
    {
      @s_line = split / \s*/,$apid_kwd_lines[$j];
      $tmp_lines[$j]=sprintf("%-3.3d%-1.1d  %-3.3s  %-8.8s  %-35.35s  %-3.3s  %-1.1s  %-2.2s  %-3.3s %-1.1s  %-8.8s \n",$s_line[3],$s_line[4],$s_line[0],$s_line[1],$s_line[2],$s_line[3],$s_line[4],$s_line[5],$s_line[6], $s_line[7],$s_line[8]); 
      $j++;
    }
  }
###################################
# Sort List                       #
###################################
  if ( $sort_option ne "0")
  {
    @ord_lines= sort@tmp_lines;
    $j = 0;
    while($ord_lines[$j])
    {
      @s_oline = split / \s*/,$ord_lines[$j];
      shift @s_oline;
      printf(OUTFILE "%-3.3s  %-8.8s  %-35.35s  %-3.3s  %-1.1s  %-2.2s  %-3.3s %-1.1s  %-8.8s \n",$s_oline[0],$s_oline[1],$s_oline[2],$s_oline[3],$s_oline[4],$s_oline[5],$s_oline[6], $s_oline[7],$s_oline[8]);
      $j++;
    }
  }
}
 
###############################################################################
# sub get_conversion_flag: get conversion flag value for keyword lines. 
###############################################################################
sub get_conversion_flag
{
  my $mnm= shift;
  my ($i, $found_conv);
  $found_conv=0;
  $i=0;
  while ($all_tlm_lines[$i])
  {
    @s_line=split(',',$all_tlm_lines[$i] );
    if ( $s_line[1] eq $mnm )
    {
      $conv_value= $s_line[5];
      last;
    }
    if ( $i == $all_tlm_count - 1 )
    {
#######print "make_hkpdf.pl:**Error: no conversion flag for $mnm in TLM lines, so will set value to <?> \n";
      $conv_value="?";
      last;
    }
    $i++;
  }
  $i=0;
  while ($all_tlm_e[$i])
  {
    @s_line=split(',',$all_tlm_e[$i] );
    if ( $s_line[1] eq "E_$mnm" )
    {
#print "found E_\n";
      $conv_value= $s_line[5];
      $found_conv=1;
      last;
    }
    $i++;
  }
  if ( $found_conv eq "0")
  {
    $i=0;
    while ($all_tlm_t[$i] )
    {
      @s_line=split(',',$all_tlm_t[$i] );
#print "s-line1 is $s_line[1] and T_mnm is T_$mnm \n";
      if ( $s_line[1] eq "T_$mnm" )
      {
#print "found T_\n";
        $conv_value= $s_line[5];
        last;
      }
      $i++;
    }
  }
}
###############################################################################
# sub get_keyword_value: get keyword value for keyword lines. 
###############################################################################
sub get_keyword_value
{
  my $mnm= shift;
  my ($i);
  $i=0;
  while ($all_fits_lines[$i])
  {
    @s_line=split(',',$all_fits_lines[$i] );
    if ( $s_line[1] eq $mnm || $s_line[1] eq "E_$mnm" )
    {
      $kwd_value= $s_line[2];
      last;
    }
    if ( $i == $all_fits_count - 1 )
    {
######print "make_hkpdf.pl:**Error: no keyword value for $mnm in FITS lines-so will set value to <????> \n";
      $kwd_value="????";
      last;
    }
    $i++;
  }
}
###############################################################################
# sub create_acon_lines: add ACON lines to file
###############################################################################
sub create_acon_lines
{
  my ($i, $j, $s_kwdline,$s_algline, $aline);
  $i=0;
  while( $apid_kwd_lines[$i] )
  {
      @s_kwdline = split / \s*/,$apid_kwd_lines[$i];
      $j=0;
      while ( $all_alg_lines[$j] )
      {
         @s_algline = split(",", $all_alg_lines[$j] );
         if ( "E_$s_kwdline[2]" eq  $s_algline[1] )
         {
           if ( $s_algline[2] eq 5) 
           {
               $s_algline[8]= 0;
           }
           if ( $s_algline[2] eq 4 )
           {
               $s_algline[8] =0;
               $s_algline[7] =0;
           }
            
           $aline[$j]=sprintf("%-4s %-30s %-1s %-2s %-2s %-2s %-2s %-2s %-2s %-8s","ACON",$s_kwdline[2],$s_algline[2],$s_algline[3],$s_algline[4],$s_algline[5],$s_algline[6],$s_algline[7],$s_algline[8],$s_algline[10]);
           print OUTFILE "$aline[$j]";
           last;
         }
         $j++;
      }
      $i++;
  }
}

###############################################################################
# sub create_dcon_lines: add DCON lines to file
###############################################################################
sub create_dcon_lines
{
 ####(1) Compare- if equal then
 my ($i, $j,$k, $s_kwdline,$s_dscline, $dline,$s_tlmtline);
 $i=0;
 while( $apid_kwd_lines[$i] )
 {
  @s_kwdline = split / \s*/,$apid_kwd_lines[$i];
  $j=0;
  while ( $all_tlm_t[$j])
  {
    @s_tlmtline=split(',',$all_tlm_t[$j] );
    if ( $s_kwdline[2] eq substr($s_tlmtline[1],2))
    {
      ## Get value in sixth field
      $lookup_value= $s_tlmtline[6];
      $k=0;
      while ($all_dsc_lines[$k] )
      {
        @s_dscline = split(",", $all_dsc_lines[$k] );
        # Use value in sixth field to compare with DSC's 1th field then if equal
        if( $s_tlmtline[6] eq  $s_dscline[1])
        {
          ## Use DSC's 2, 3, and 4.
          $dline[$k] = sprintf("%-4s %-30s %-5s %-5s %-35s %-8s", "DCON", $s_kwdline[2], $s_dscline[2] ,$s_dscline[3], $s_dscline[4], $s_dscline[6] );
          print OUTFILE "$dline[$k]";
        }
        $k++;
      }#loop thru dsc lines
      last;
    }# if found kwd match to tlm line
    else
    {
      $j++;
    }
  } #loop thru TLM,T_ lines for keyword match
  $i++;
 }# loop thru keywords in KWD line        
}

###############################################################################
#sub create_hkpdf_for_1AF:Create config format for 1AF apid since does not exist
###############################################################################
sub create_hkpdf_for_1AF
{
  my($directory, $fileline,$apidline,$kwdline1,$kwdline2);
  $directory=$ENV{'HK_CONFIG_DIRECTORY'};
# Check if version_number folder exists
  opendir(DIR,"$directory$version_number") || die "Can't open: $!\n";
  @items=readdir(DIR);
  closedir(DIR);
  $found_file=0;
# check if file already there
  foreach $file ( @items )
  {
    if ($file eq "apid-1AF-version-$version_number")
    {
      $found_file=1;
      last;
    }
  }
# If file does not exist create it & also if does exist create with latest code
  if ( $found_file ne 1  ||  $found_file eq 1 )
  {
#   create file line    
    $fileline=sprintf("%-4s %-16s  %-7s  %-10s %-8s\n","FILE","MADE_UP_BY_CARL",$version_number, "2006/99/99", "99:99:99");

#   create apid line    
    $apidline=sprintf("%4.4s %5.5s  %3.3s  %7.7s \"%-s\" %8.8s\n","APID","0x1AF","4", "HMI", "PACKET CREATE FOR WORKAROUND", "20069999");

#   create keyword line 1
    $kwdline1=sprintf("%-3s  %-8s %-20s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s\n","KWD","HTCS1AF","HTCS1AF","0","0","32","UL1","R","20069999"); 

#   create keyword line 2
    $kwdline2=sprintf("%-3s  %-8s %-20s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s\n","KWD","HTCSS1AF","HTCSS1AF","4","0","32","UL1","R","20069999"); 

#   Create file name
    $filename="apid-1AF-version-$version_number";

#   Create directory name using environment variable, version number and filename
    open(OUTFILE, ">$directory$version_number/$filename") || die "Can't Open  $directory/$filename: $!\n" ;
    print OUTFILE "$fileline";
    print OUTFILE "$apidline";
    print OUTFILE "$kwdline1";
    print OUTFILE "$kwdline2";
    close (OUTFILE);
  }
}

