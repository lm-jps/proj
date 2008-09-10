#!/usr/bin/perl
#############################################################################
#NAME:   make_hkdpf.pl - Make Housekeeping Data Packet format files
#Author: Carl
#Purpose:Strip data from STANFORD_TLM_HMI_AIA.txt file or/and
#        GODDARD_TLM_SDO.txt file and create HKPDF files using the filename 
#        format apid-<#>-version-<#>. These HKPDF files contain the required 
#        and minimum needed data for processing housekeeping data packets by 
#        the Level 0 Software Applications.
#Execution for help: make_hkpdf.pl -h 
#Execution to create hkpdf files: 
#  make_hkpdf.pl sort=<sort-num> apidlist=<list|keyword>
#  make_hkpdf.pl sort=<sort-num> apidfile=<full path to file containing apids>
#Example Execution commands:
#  make_hkpdf.pl sort=4 apidlist="0x081"
#  make_hkpdf.pl sort=4 apidlist="0x1BD 0x1DB"
#  make_hkpdf.pl sort=4 apidlist="0x1BD 0x1DB 0x081"
#  make_hkpdf.pl sort=4 apidlist="ALL_SDO"
#  make_hkpdf.pl sort=4 apidlist="ALL_HK
#  make_hkpdf.pl sort=4 apidlist="ALL"
#Create Date:3/15/2008 
#############################################################################
# main function
#############################################################################
#set environment variables and initialize variables.
$hm=$ENV{'HOME'};
$script_version=1.2;

# HK environment variable for processing STANFORD file
$ENV{'HK_STHA_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/fromlmsal";
$ENV{'HK_STHA_FILENAME'}="STANFORD_TLM_HMI_AIA.txt";
$ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file/";

# SDO environment variables for processing GODDARD file
$ENV{'HK_SH_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/sdo_hk_config_file/";
$ENV{'HK_GTS_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/fromgoddard";
$ENV{'HK_GTS_FILENAME'}="GODDARD_TLM_SDO.txt";

#common setting for all environments
$ENV{'MAILTO'}="";
$script_dir="$hm/cvs/JSOC/proj/lev0/scripts/hk";
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:$script_dir";

#logfile
$logfile="$hm/cvs/JSOC/proj/lev0/scripts/hk/log-clmuq";

#setup contants
use constant NO_PROCESSING =>  0;
use constant HK_PROCESSING =>  1;
use constant SDO_PROCESSING => 2;

#display environment variables
#print "--->make_hkpdf enviroment variables:\n";
#print `env `;

# check arguments passed
&check_arguments();

# set apid values to process based on arguments passed
if ($apid_list_flg == "1")
{
  ($apid_arg)=&set_apidlist_value($ARGV[1]);
}
elsif ($apid_file_flg == "1")
{
  ($apid_arg)=&set_apidfile_value($ARGV[1]);
}
else
{
   print "make_hkpdf.pl:ERROR: argument missing. Need to use apidlist or apidfile arguments!\n";
}

# check if need to process STANFORD or/and GODDARD file based on APIDS set in arguments
# passed in apidlist or apidfile.  
($hk_type_flag,$sdo_type_flag)=&check_types_to_do($apid_arg);
if ($ret=&check_processing_flag($hk_type_flag) == HK_PROCESSING)
{
  #execute processing STANFORD file
  print "Doing Processing of STANFORD file\n";
  #make files for HK_PROCESSING
  &make_files(HK_PROCESSING);
  #LOG
  print "Finished processing STANFORD file using make_hkpdf.pl \n";
}
if ($ret=&check_processing_flag($sdo_type_flag) == SDO_PROCESSING)
{
  #execute processing GODDARD file
  print "Doing Processing of GODDARD file\n";
  #make files for SDO_PROCESSING
  &make_files(SDO_PROCESSING);
  #LOG
  print "Finished processing GODDARD file using make_hkpdf \n";

}

###############################################################################
# sub make files:  
###############################################################################
sub make_files($)
{
  
  # local variables and intialize variables
  my($processing_type);
  $processing_type=$_[0];

  # (1) get all lines in STHA file
  &get_all_lines($processing_type);

  # (2) get all apids to process based on environment variable 
  &get_apids_to_do();
  #print "Please wait creating  < $all_apids_count > hkdpf files\n";

  # (3) create HKPDF filenames
  &create_hkpdf_filenames($processing_type);

  # (4) create HKPDF files
  &create_hkpdf_files;

  # (5) print stats:total of apid files created
  #print ".. completed creating < $all_apids_count > hkdpf files\n";

  # (6) create workaround file apid-1AF-version-<current version number>
  #&create_hkpdf_for_1AF; not used anymore

}


###############################################################################
# sub check type flag: check processing flag 
###############################################################################
sub check_processing_flag($flag)
{
  if ($_[0] == HK_PROCESSING )
  {
     #then bit one is set so return setting of HK_PROCESSING
     return(HK_PROCESSING);
  }
  elsif ($_[0] == SDO_PROCESSING )
  {
     #then bit two is set so return setting of HK_PROCESSING
     return(SDO_PROCESSING);
  }
  else
  {
     return($NO_PROCESSING);
  }
}
###############################################################################
# sub get_all_lines: gets all lines in STANFORD_TLM_HMI_AIA.txt file
###############################################################################
sub get_all_lines($)
{
  # initalize arrays
  @all_lines=();
  @all_pkt_lines=();
  @all_apids=();
  @all_fits_lines=();
  @all_tlm_lines=();
  @all_tlm_e=();
  @all_tlm_t=();
  @all_alg_lines=();
  @all_dsc_lines=();
  
  # determine directory and filename to use based on processing  either GODDARD file or STANFORD file
  my($processing_type);
  $processing_type=$_[0];
  if($processing_type == SDO_PROCESSING)
  {
    $dn=$ENV{'HK_GTS_DIRECTORY'};
    $fn=$ENV{'HK_GTS_FILENAME'};
  }
  elsif ($processing_type == HK_PROCESSING)
  {
    $dn=$ENV{'HK_STHA_DIRECTORY'};
    $fn=$ENV{'HK_STHA_FILENAME'};
  }
  else
  {

    print "make_hkpdf.pl:ERROR: Did not find processing type. Check values.\n";
    exit;
  }

  # open either STANFORD file or GODDARD file for processing
  open(FILE, "$dn/$fn") || die "Can't Open $fn file: $!\n";
  while (<FILE>)
  {
    push(@all_lines, $_) ;
    @s_line=split(',',$_ );
    if ( $s_line[0] eq "PKT")
    {
      push(@all_pkt_lines, $_) ;
      if($processing_type == HK_PROCESSING && hex $s_line[1] < 63 ||  hex $s_line[1] > 399)
      {
        # if processing STANFORD file push only these apids
        push(@all_apids, $s_line[1]) ;
      }
      elsif ($processing_type == SDO_PROCESSING && hex $s_line[1] > 95 &&  hex $s_line[1] < 400)
      {
        # if processing GODDARD file push only these apids
        push(@all_apids, $s_line[1]) ;
      }
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
  #statistics
  #print "LINES=$line_count\nPKTS=$pkt_line_count\nAPIDS=$all_apids_count\n";
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

###########################################
#  check types to do                      #
###########################################
sub check_types_to_do($apid_arg)
{
  $sdo_type_flag=NO_PROCESSING;
  $hk_type_flag=NO_PROCESSING;
  if ( $apid_arg eq "ALL" || $apid_arg eq "")
  {
     $sdo_type_flag=SDO_PROCESSING;#set bit 2
     $hk_type_flag=HK_PROCESSING;  #set bit 1
  }
  elsif ( $apid_arg eq "ALL_HK")
  {
     $sdo_type_flag=NO_PROCESSING; #set bit 2
     $hk_type_flag=HK_PROCESSING;  #set bit 1
  }
  elsif ( $apid_arg eq "ALL_SDO")
  {
     $sdo_type_flag=SDO_PROCESSING; #set bit 2
     $hk_type_flag=NO_PROCESSING;  #set bit 1
  }
  else
  {
     #check list of apids to deterine whether how to set type flag
    @arg_items=split(' ',$apid_arg );
    while (@arg_items)
    {
      $num= hex $arg_items[0];
      if ($num > 95 and $num < 400)
      {
         $sdo_type_flag=SDO_PROCESSING;
      }
      elsif ($num > 0 and $num < 64)
      {
         $hk_type_flag=HK_PROCESSING;

      }
      elsif ($num > 399 and $num < 600)
      {
         $hk_type_flag=HK_PROCESSING;

      }
      shift @arg_items;
    }
   }
   return ($hk_type_flag,$sdo_type_flag);
}
###############################################################################
# sub check_apid_to_do: check which apid to do, if ALL do all in STHA file
###############################################################################
sub get_apids_to_do
{
  #$arg=$ENV{'HK_APID_LIST'};
  $arg=$apid_arg;
  if ( $arg eq "ALL" || $arg eq "ALL_HK" || $arg eq "ALL_SDO")
  {
     $apids_flag=ALL;
  }
  else 
  {
     $apids_flag=SOME;
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
  #print "APID List to process is as follows: @all_apids \n";
}

###############################################################################
# sub create_hkpdf_filenames: create filenames for each apid.
###############################################################################
sub create_hkpdf_filenames($)
{
  # local variables and initialized variables
  my($i, $j, $processing_type);
  @filenames=();
  $processing_type= $_[0];

  # get version number
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
  $j=0;
  $pkt_count=@all_pkt_lines;
  while ( $all_apids[$i])
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
    if ( $j eq $pkt_count)
    {
      push(@filenames, sprintf("not-used"));
        $i++;
        $j=0;
    }
  }
}
###############################################################################
# sub create_hkpdf_files: create files for each apid to do.
###############################################################################
sub create_hkpdf_files($)
{
  # set local variables and initialize variables
  my($i,$processing_type);
  $processing_type= $_[0];

  # Determine directory to put HKPDF files
  if($processing_type == SDO_PROCESSING)
  {
    $directory=$ENV{'HK_SH_CONFIG_DIRECTORY'};
  }
  elsif ($processing_type == HK_PROCESSING)
  {
    $directory=$ENV{'HK_CONFIG_DIRECTORY'};
  }
  else
  {
    print "make_hkpdf.pl:ERROR: Did not find processing type. Exiting.\n";
    exit;
  }

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

  $i=0;
  while ($all_apids[$i])
  {
    if($processing_type == HK_PROCESSING && hex $all_apids[$i] > 95 && hex $all_apids[$i] < 400)
    {
        $i++;
    }
    elsif($processing_type == SDO_PROCESSING && (hex $all_apids[$i] < 64 || hex $all_apids[$i] > 399))
    {
        $i++;
    }
    else
    {
      # Create directory using environment variable, version number and filename
      open(OUTFILE, ">$directory$version_number/$filenames[$i]") || die "Can't Open  $directory/$filenames[$i]: $!\n" ;
      &create_file_line($processing_type);
      &create_apid_line( $all_apids[$i]) ;
      &create_kwd_lines( $all_apids[$i]) ;
      &create_acon_lines() ;
      &create_dcon_lines() ;
      close( OUTFILE);
      print "Finish file $i. APID is $all_apids[$i]\n";
      $i++;
    }
  }
}

###############################################################################
# sub create_pkt_line: create FILE line for hkpdf file.
###############################################################################
sub create_file_line($)
{
  # set local variables and initialize variables
  my($i,$f_line,$fileline,$processing_type);
  $i=0;
  $processing_type= $_[0];

  # create file line
  while ($all_lines[$i])
  {
    @f_line=split(' ', $all_lines[$i]);
    if ($f_line[1] eq "File" && $f_line[2] eq "RCS")
    {
      last;
    }
    $i++;
  }        
  # Determine if getting GODDARD or STANFORD filename inside file
  if($processing_type == SDO_PROCESSING)
  {
    $fileline=sprintf("%4.4s %19.19s  %-7s %-10s  %-8s\n","FILE",$f_line[5],$f_line[6], $f_line[7], $f_line[8]); 
  }
  elsif ($processing_type == HK_PROCESSING)
  {
    $fileline=sprintf("%4.4s %24.24s  %-7s %-10s  %-8s\n","FILE",$f_line[5],$f_line[6], $f_line[7], $f_line[8]); 
  }

  # print time created file and version of this script used to create file
  $datestamp=`date`;
  $datestamp =~ s/\n//g;
  print OUTFILE "# hkdpf file created by make_hkdpf.pl script.\n# date created file: $datestamp\n# script version: $script_version\n";
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
    if ( hex $p_line[1] eq hex $find_apid )
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
  elsif ( hex $p_line[1] >= 96 && hex $p_line[1] <= 399)
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "SDO", $p_line[5], $p_line[6]); 
  }
  else
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], "UNKN", $p_line[5], $p_line[6]); 
  }


  #$apidline=sprintf("%4.4s %5.5s  %3.3s  %8.8s\n","APID",$p_line[1],$p_line[3], $p_line[6]); 
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
###############################################################################
#sub check_arguments
###############################################################################
sub check_arguments
{

  $help_flg= "0";
  $sort_flg= "0";
  $apid_list_flg="0";
  $apid_file_flg="0";
  
  if ($#ARGV < 1 || $#ARGV > 1)
  {
    $help_flg="1";
  }
    
  if ($#ARGV >= 0)
  {
    if ($ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "-H")
    {
      $help_flg = "1";
    }
    elsif (substr($ARGV[0],0,5) eq "sort=" )
    {
      #use file to get list of apids to create map files for
      $sort_flg = "1";
      &set_sort_value($ARGV[0]);
    }
    else
    {
       print  "ERROR: Did not use -h or sort arguments. Exiting script\n\n";
       &show_help_info("USAGE");
    }
  }
  if ($#ARGV >= 1)
  {
    if (substr($ARGV[1],0,9) eq "apidlist=" )
    {
      #use file to get list of apids to create map files for
      $apid_list_flg = "1";
      #($apid_arg)=&set_apidlist_value($ARGV[1]);
    }
    elsif (substr($ARGV[1],0,9) eq "apidfile=" )
    {
      #use apid following -a to create map files 
      #push decimal character value in this format dddd. Example: -a 0001
      $apid_file_flg = "1";
      #&set_apidfile_value($ARGV[1]);
    }
    else
    {
       print  "ERROR: Did not use apidfile or apidlist arguments. Exiting script\n\n";
       &show_help_info("USAGE");
    }
  }

  if ( $help_flg eq "1")
  {
     &show_help_info("HELP");
  }
}
###############################################################################
#sub set_apidlist_value
###############################################################################
sub set_apidlist_value($)
{
  my($arg,$s_arg);
  $arg=$_[0];
  @s_arg=split("=",$arg );
  if ( $s_arg[1] eq "")
  {
    print "make_hkpdf.pl:ERROR:exiting bad apidlist argument:<$s_arg[1]>\n";
    exit;
  } 
  elsif($s_arg[1] eq "ALL") 
  {
    $apidlist_value=$s_arg[1];
    print "Doing all apids in STANFORD and GODDARD files\n";
  }
  elsif($s_arg[1] eq "ALL_HK") 
  {
    $apidlist_value=$s_arg[1];
    print "Doing all apids in STANFORD file only\n";
  }
  elsif($s_arg[1] eq "ALL_SDO") 
  {
    $apidlist_value=$s_arg[1];
    print "Doing all apids in GODDARD file only\n";
  }
  else 
  {
    $apidlist_value=$s_arg[1];
  }
  #print "apidlist value is <$apidlist_value>\n";
  return($apidlist_value);
}
###############################################################################
#sub set_apidlist_value
###############################################################################
sub set_apidfile_value($)
{
  my($arg,$s_arg);
  $arg=$_[0];
  @s_arg=split("=",$arg );
  if ( $s_arg[1] eq "")
  {
    print "make_hkpdf.pl:ERROR:exiting bad apidfile argument:<$s_arg[1]>\n";
    exit;
  } 
  else 
  {
    $apidfile_value=$s_arg[1];
    # open apidfile passed in as argument
    open(FILEAPID, "$s_arg[1]") || die "Can't Open $s_arg[1] file: $!\n";
    while (<FILEAPID>)
    {
      #read in apids 
      $_ =~ s/\n//g; #regular exp rm cr
      push(@all_apid_from_file, $_) ;
    }
    close(FILEAPID);
    #print "apids from file: @all_apid_from_file\n";
    # put in format for processing... 0x1BD 0x081
    $apidfile_value=join(' ',@all_apid_from_file);
    #return back values
  }
  return  $apidfile_value;
}
###############################################################################
#sub set_sort_value
###############################################################################
sub set_sort_value($)
{
  my($arg);
  $arg=$_[0];
  # Set sort option to default and check for argument passed
  $sort_option="0";
  if ( $arg eq "sort=0" )
  {
    $sort_option="0";
  }
  elsif ( $arg eq "sort=1" )
  {
    $sort_option="1";
  }
  elsif ( $arg eq "sort=2" )
  {
    $sort_option= "2";
  }
  elsif ( $arg eq "sort=3" )
  {
    $sort_option="3";
  }
  elsif ( $arg eq "sort=4" )
  {
    $sort_option="4";
  }
  elsif ( $arg eq "" )
  {
    $sort_option="0";
  }
  else
  {
    print "**make_hk.pl:Error bad argument for sort parameter: choice 0, 1, 2, 3, 4 or no argument \n";
    print "where sort-choice equal to 0 is no sort \n";
    print "where sort-choice equal to 1 is sort by keyword\n";
    print "where sort-choice equal to 2 is sort by byte position\n";
    print "where sort-choice equal to 3 is sort by long telem name\n";
    print "where sort-choice equal to 4 is sort by byte position and bit start position\n";
    exit;
  }
}
###############################################################################
#sub show_help_info
###############################################################################
sub show_help_info($)
{
  my ($help_arg);
  $help_arg=$_[0];
  
  
  if($help_arg eq "HELP")
  {
    print "Help Listing\n";
  }
  else
  {
    print "Usage:\n";
  }
  print "(1)Ways to Execute Perl Script: \n";
  print "(1a)Help Information : make_hkdpf.pl -h\n";
  print "(1b)Create Files: make_hkdpf.pl sort=<0,1,2,3 or 4> apidfile=<full path to file>\n";
  print "(1c)Create Files: make_hkdpf.pl sort=<0,1,2,3 or 4> apidlist=<list|keyword>\n";
  print "where sort-choice equal to 0 is no sort \n";
  print "where sort-choice equal to 1 is sort by keyword\n";
  print "where sort-choice equal to 2 is sort by byte position\n";
  print "where sort-choice equal to 3 is sort by long telem name\n";
  print "where sort-choice equal to 4 is sort by byte and bit position\n";
  print "************************************\n";
  print "** Note:Currently using option 4! **\n";
  print "************************************\n";
  print "where apidfile is full path and filename containing list of apids to do\n";
  print "where apidlist contains list of apids to do(i.e.,0x1BD 0x081 ).\n";
  print "where apidlist contains keyword (i.e., ALL, ALL_HK, or ALL_SDO)on apids to do\n";
  if($help_arg eq "HELP")
  {
    print "(2) Environment variable HK_STHA_DIRECTORY\n";
    print "(2a) Used to store location of STANFORD_TLM_HMI_AIA.txt file.\n";
    print "(3) Environment variable HK_STHA_FILE\n";
    print "(3a) Used to store filename of STANFORD_TLM_HMI_AIA.txt file to process.\n";
    print "(3b) This file is input data for this script\n";
    print "(4) Environment variable HK_CONFIG_DIRECTORY\n";
    print "(4a) Used to store location of hk config files for each apid\n";
    print "(4b) These files are the output for this script.\n";
    print "(4c) The filename format is apid-<apid#>-version-<version#>.\n";
    print "(5) Environment variable HK_GTS_DIRECTORY\n";
    print "(5a) Used to store location of GODDARD_TLM_SDO.txt file.\n";
    print "(6) Environment variable HK_GTS_FILE\n";
    print "(6a) Used to store filename of GODDARD_TLM_SDO.txt file to process.\n";
    print "(6b) This file is input data for this script\n";
    print "(7) Environment variable HK_SH_CONFIG_DIRECTORY\n";
    print "(7a) Used to store location of sdo-hk config files for each apid\n";
    print "(7b) These files are the output for this script.\n";
    print "(7c) The filename format is apid-<apid#>-version-<version#>.\n";
  }
  exit;
}
