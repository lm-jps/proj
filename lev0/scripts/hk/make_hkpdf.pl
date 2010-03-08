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
#  make_hkpdf.pl sort=4 apidlist="ALL_SDO" (Normally use)
#  make_hkpdf.pl sort=4 apidlist="ALL_HK"  (Normally use)
#  make_hkpdf.pl sort=4 apidlist="ALL"
#Create Date:3/15/2008 
#############################################################################
# main function
#############################################################################
#set environment variables and initialize variables.
$hm=$ENV{'HOME'};
$script_version=1.6;

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
use constant HMI_AIA_MERGED_SERIES_FVN_DN => 160;
use constant HMI_AIA_MERGED_SERIES_FVN_WN => 1;
use constant KERNAL_MERGED_SERIES_FVN_DN  => 161;
use constant KERNAL_MERGED_SERIES_FVN_WN  => 1;
use constant MERGED_SERIES => 1;
use constant NOT_MERGED_SERIES => 0;
use constant KERNAL_APID => 1;
use constant NOT_KERNAL_APID => 0;

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
  @column_table={};
  
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
  open(FILE, "$dn/$fn") || die "1-Can't Open $fn file: $!\n";
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
    if ( $s_line[0] eq "SRC")
    {
      # strip VER_NUM and KER_VER_NUM lines 
      # set column_table to use later to set column id in get_column function
      if((index $_, "KER_VER_NUM") != -1 ) 
      {
        @line=split(',',$_);
        $tmpapid=$line[3];
        $numapid= hex $tmpapid;
        if(($numapid <= 31 && $numapid > 0) || ($numapid >= 400 && $numapid < 500))
        { 
           $column_table{$numapid}="HMIKER";
        }
        elsif(($numapid > 31 && $aid < 63) || ($numapid >= 500 && $numapid < 600))
        {
           $column_table{$numapid}="AIAKER";
        }

      }
      elsif((index $_, "VER_NUM") != -1 ) 
      {
        @line=split(',',$_);
        $tmpapid=$line[3];
        $numapid= hex $tmpapid;
        if(($numapid <= 31 && $numapid > 0) || ($numapid >= 400 && $numapid < 500))
        {  
          $column_table{$numapid}="HMI";
        }
        elsif(($numapid > 31 && $numapid < 63) || ($numapid >= 500 && $numapid < 600))
        {
           $column_table{$numapid}="AIA";
        }
        elsif($numapid >= 96 && $numapid < 400)
        {
           $column_table{$numapid}="SDO";
        }
        elsif($numapid >= 2002 && $numapid <= 2047)
        {
           $column_table{$numapid}="SSIM";
        }
        else 
        {
           $column_table{$numapid}="UNKN";
        }
      }
    }
    if ( $s_line[0] eq "FITS" )
    {
      push(@all_fits_lines, $_) ;
      next;
    }
    if ( $s_line[0] eq "TLM" )
    {
      push(@all_tlm_lines, $_) ;

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
      #check size of variable typed in
      $sv=length($arg_items[0]);
      if($sv > 5) #new check added 6-23-2009:cc
      {
        print "ERROR:Retype in argument for hex apid value for $arg_items[0] using 5 characters! ";
        print "Exiting script. Retry.\n";
        exit;
      }
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
  my($i, $j, $processing_type, $tmpapid);
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
      print "DOING VERSION NUMBER : $version_number \n";
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
      $j++; #if pkt line does not contain current apid doing go to next pkt line
    }
    else
    {
      if(SDO_PROCESSING == $processing_type)
      {
         # if SDO type- create files like new way 
        $a= hex $apid_line[1];
        push(@filenames, sprintf("SDO-%3.3s-version-%s\n",  "ASD", $version_number));
        $i++;
      }
      elsif(&check_for_merged_version($version_number,hex $apid_line[1])) #CC:Add block to handle merge series requirement
      { #checks only for HK version numbers
        #create name of filename for merged series requirements (HMI-ISP-version-1.163,AIA-ISP-version-1.163)
        $a= hex $apid_line[1];
        if (hex $apid_line[1] == 445 )
        {
          # push this filename to array. will use this config to do all hmi isps
          $a= hex $apid_line[1];
          push(@filenames, sprintf("HMI-%3.3s-version-%s\n",  "ISP", $version_number));
        }
        elsif ( hex $apid_line[1] == 529)
        {
          # push filename for this apid to array. will use this config to do all aia isps.
          $a= hex $apid_line[1];
          push(@filenames, sprintf("AIA-%3.3s-version-%s\n",  "ISP", $version_number));
        }
        elsif ( hex $apid_line[1] == 451)
        {
          # push these filename to array
          $a= hex $apid_line[1];
          push(@filenames, sprintf("HMI-%3.3s-version-%s\n",  "SEQ", $version_number));
        }
        elsif ( hex $apid_line[1] == 536)
        {
          # push these filename to array
          $a= hex $apid_line[1];
          push(@filenames, sprintf("AIA-%3.3s-version-%s\n",  "SEQ", $version_number));
        }
        elsif ( hex $apid_line[1] == 448)
        {
          # push these filename to array
          $a= hex $apid_line[1];
          push(@filenames, sprintf("HMI-%3.3s-version-%s\n",  "OBT", $version_number));
        }
        elsif ( hex $apid_line[1] == 129)
        {
          # push these filename to array
          $a= hex $apid_line[1];
          push(@filenames, sprintf("SDO-%3.3s-version-%s\n",  "ASD", $version_number));
        }
        elsif ( hex $apid_line[1] == 540)
        {
          # push these filename to array
          $a= hex $apid_line[1];
          push(@filenames, sprintf("AIA-%3.3s-version-%s\n",  "OBT", $version_number));
        }
        elsif (hex $apid_line[1] == 569 || hex $apid_line[1] == 475 || hex $apid_line[1] == 39 ||  hex $apid_line[1] == 29 ||
               hex $apid_line[1] == 481 || hex $apid_line[1] == 576 || hex $apid_line[1] == 21 ||  hex $apid_line[1] == 46 ||
               hex $apid_line[1] == 478 || hex $apid_line[1] == 580 || hex $apid_line[1] == 18 ||  hex $apid_line[1] == 50)
        {
          # push null for these filenames to array since already did primary apid for isp, seq, etc.
          $tmpapid= hex $apid_line[1];
          push(@filenames, "");
          print "make_hkpdf.pl:WARNING-Skipping creating HKPDF files for apid. Probably got ISP, SEQ or OBT apids that are duplicates of apid 0x1BD or 0x211, etc . apid is $tmpapid\n";
        }
        elsif((hex $apid_line[1] <= 31  && hex $apid_line[1] > 0) || (hex $apid_line[1] >= 400 && hex $apid_line[1] < 500)) 
        {
          # do others hmi apids except isp and seq. using new file format using "HMI"  instead of "apid"
          @apid_line_hex=split('x',$apid_line[1]);
          push(@filenames, sprintf("HMI-%3.3s-version-%s\n",  $apid_line_hex[1], $version_number));
        }
        elsif((hex $apid_line[1] > 31 && hex $apid_line[1] < 65 ) || (hex $apid_line[1] >= 500 && hex $apid_line[1] < 600)) 
        {
          # do others aia apids except isp and seq. using new file format using "AIA" instead of "apid"
          @apid_line_hex=split('x',$apid_line[1]);
          push(@filenames, sprintf("AIA-%3.3s-version-%s\n",  $apid_line_hex[1], $version_number));
        }
        else
        {
           $tmpapid=hex $apid_line[1];
           print "make_hkpdf.pl:WARNING--Skipping creating HKPDF file for apid. Probably got sim apids or ADP apid. apid is $tmpapid\n";
        }
        $i++;
      }#end if merged case
      else
      {
        # if not merged case use apid-#-version-# filename format
        @apid_line_hex=split('x',$apid_line[1]);
        push(@filenames, sprintf("apid-%3.3s-version-%s\n",  $apid_line_hex[1], $version_number));
        $i++;
      }
    }
    if ( $j eq $pkt_count)
    {
      push(@filenames, sprintf("not-used"));
      print "WARNING:Entered possible bad value for apid:$all_apids[$i].\n1-Redo execution of script without this apid if apid was incorrectly entered\n2-Remove not-used file from folder which got created.\n3-If good apid, try rerun with this \"apid only\" via make_hkpdf sort=4 apidlist=\"<hex apid value for one apid>\"\n";
      $i++;
      $j=0;
    }
  }#while loop
}
###############################################################################
# sub create_hkpdf_files: create files for each apid to do.
###############################################################################
sub create_hkpdf_files($)
{

  # initilize warning notification flag
  $warning_flag=-1;

  # set local variables and initialize variables
  my($i,$processing_type, $files_did);
  $processing_type= $_[0];
  $files_did=1;

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
    # filter out apids that are not HK or SDO-HK apids
    if($processing_type == HK_PROCESSING && hex $all_apids[$i] > 95 && hex $all_apids[$i] < 400)
    {
        $i++;
    }
    elsif($processing_type == SDO_PROCESSING && (hex $all_apids[$i] < 64 || hex $all_apids[$i] > 399))
    {
        $i++;
    }
    elsif($filenames[$i] eq "")  #CC:Added-case-1-28-2009
    {   
        $i++;#skip values not set because not needed.do apid 445 but not 475.
    }
    else
    {
      # Create directory using environment variable, version number and filename
      open(OUTFILE, ">$directory$version_number/$filenames[$i]") || die "2-Can't Open  $directory/$filenames[$i]: $!\n" ;
      &create_file_line($processing_type);
      &create_apid_line( $all_apids[$i]) ;
      &create_kwd_lines( $all_apids[$i]) ;
      #check if bit position is within range- if not exit.
      $cbp_status=check_bit_positions($processing_type,@apid_kwd_lines);#added check on 10-12-2009
      if($cbp_status == -1)
      {
        close(OUTFILE);
        $delfile="$directory$version_number/$filenames[$i]";
        $delfile=~ s/\n//;
        unlink($delfile) || die "After got error can't delete file <$delfile>: $!\n";
        exit;
      }
      &create_acon_lines() ;
      &create_dcon_lines() ;
      close( OUTFILE);
      #print "Finish file $i. APID is $all_apids[$i]\n";
      print "Finished file where APID is $all_apids[$i]. Files did:$files_did\n";

      #check if need to notify user to run patch01.pl script.
      $warning_flag=check_patch01($version_number,$all_apids[$i]);
      &patch_notice($warning_flag,$version_number,$all_apids[$i]);

      #do next apid
      $files_did++;
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
  my $find_apid =shift ;#passed in apid array
  my($i, $p_line, $apidline, $colid);
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
  # get column id field to use
  $colid=get_ground_column_value(hex $p_line[1]);
  $len_colid=length($colid) + 1;
  $variable_setting=sprintf("%s%d.%d%s","%",$len_colid,$len_colid,"s");

  if ( hex $p_line[1] >= 1 && hex $p_line[1] <= 31)
  {
      $apidline=sprintf("%4.4s %5.5s %3.3s  $variable_setting \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], $colid, $p_line[5], $p_line[6]); 
  }
  elsif ( hex $p_line[1] >= 400 && hex $p_line[1] <= 486)
  {
    #fix APID lines APID value for apid 445,529,448 only if equal or greater than the merged name threshold
    if((hex $find_apid == 445 || hex $find_apid == 451 || hex $find_apid == 448) && (check_for_merged_version($version_number,hex $find_apid) == MERGED_SERIES ))
    {
      $apidline=sprintf("%4.4s %3.3s %3.3s  $variable_setting  \"%-s\" %8.8s\n","APID",&fix_apid_value(hex  $find_apid),$p_line[3], $colid, $p_line[5], $p_line[6]); 
    }
    else
    {
      $apidline=sprintf("%4.4s %5.5s %3.3s  $variable_setting  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], $colid, $p_line[5], $p_line[6]); 
    }
  }
  elsif ( hex $p_line[1] > 31 && hex $p_line[1] <= 63)
  {
    $apidline=sprintf("%4.4s %5.5s %3.3s  $variable_setting  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], $colid, $p_line[5], $p_line[6]); 
  }
  elsif ( hex $p_line[1] >= 500 && hex $p_line[1] <= 580)
  {
    #fix APID lines APID value for apid 529,536,540 only if equal to or greater than  merged name threshold
    if((hex $find_apid == 529 || hex $find_apid == 536 || hex $find_apid == 540) && (check_for_merged_version($version_number,hex $find_apid) == MERGED_SERIES ))
    {
      $apidline=sprintf("%4.4s %3.3s %3.3s  $variable_setting  \"%-s\" %8.8s\n","APID", &fix_apid_value(hex  $find_apid),$p_line[3], $colid, $p_line[5], $p_line[6]); 
    }
    else
    {
      $apidline=sprintf("%4.4s %5.5s %3.3s  $variable_setting  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], $colid, $p_line[5], $p_line[6]); 
    }
  }
  else
  {
    # gets either SDO, SSIM or UNKN
    $apidline=sprintf("%4.4s %5.5s %3.3s  %4.4s  \"%-s\" %8.8s\n","APID",$p_line[1],$p_line[3], $colid, $p_line[5], $p_line[6]); 
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
          #look for keyword lines
          @line=split(',',$all_lines[$i] );
          if ( $line[0] eq "SRC" )
          {
            @mnm=split('_', $line[1]);


            if ( $mnm[0] ne "T" && $mnm[0] ne "E" && $mnm[0] ne "Y") 
            {
              &get_conversion_flag ($line[1]);
              &get_keyword_value ($line[1]);
              #handle merge series requirement
              if(&check_for_merged_version($version_number, hex $find_apid ) || hex $find_apid == 129)
              {
                #fix timecode and checksum values 
                if(substr($line[1],8,8) eq "TIMECODE")
                {
                   if( hex $find_apid == 445 || hex $find_apid == 529 || hex $find_apid == 536 ||  hex $find_apid == 451 || hex $find_apid == 448  || hex $find_apid == 540  || hex $find_apid == 129)
                   {
                     $kwdline[$j]=sprintf("%-3s %-8s  %-40s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s","KWD",&fix_timecode(hex $find_apid,$kwd_value,"SHORT"),&fix_timecode(hex $find_apid,$line[1],"LONG"),$line[4],$line[5],$line[6],$line[7],$conv_value,$line[10]); 
                   }
                   else
                   {
                     $kwdline[$j]=sprintf("%-3s %-8s  %-40s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s","KWD",$kwd_value,$line[1],$line[4],$line[5],$line[6],$line[7],$conv_value,$line[10]); 
                   }
                    
                }
                elsif(substr($line[1],8,8) eq "CHECKSUM")
                {
                   if( hex $find_apid == 445 || hex $find_apid == 529 || hex $find_apid == 536 ||  hex $find_apid == 451 || hex $find_apid == 448 || hex $find_apid == 540)
                   {
                     $kwdline[$j]=sprintf("%-3s %-8s  %-40s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s","KWD",&fix_checksum_keywords(hex $find_apid,$kwd_value,"SHORT"),&fix_checksum_keywords(hex $find_apid,$line[1],"LONG"),$line[4],$line[5],$line[6],$line[7],$conv_value,$line[10]); 
                   }
                   else
                   {
                     #if  merged version but not a timecode and checksum name then create kwd line using values in Stanford file. 
                     $kwdline[$j]=sprintf("%-3s %-8s  %-40s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s","KWD",$kwd_value,$line[1],$line[4],$line[5],$line[6],$line[7],$conv_value,$line[10]); 
                   }
                }
                else 
                {
                  $kwdline[$j]=sprintf("%-3s %-8s  %-40s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s","KWD",$kwd_value,$line[1],$line[4],$line[5],$line[6],$line[7],$conv_value,$line[10]); 
                }
              }
              else
              { #if not merged version then do all kwd lines without fixes to timecode and checksum names
                $kwdline[$j]=sprintf("%-3s %-8s  %-40s  %-3s  %-1s  %-2s  %-3s %-1s  %-8s","KWD",$kwd_value,$line[1],$line[4],$line[5],$line[6],$line[7],$conv_value,$line[10]); 
              }
              push(@apid_kwd_lines, $kwdline[$j]);
#             print OUTFILE "$kwdline[$j++]";
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
      #print "make_hkpdf.pl:**Error: no conversion flag for $mnm in TLM lines, so will set value to <?> \n";
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
      @s_line=split(',',$all_tlm_t[$i] ); #print "s-line1 is $s_line[1] and T_mnm is T_$mnm \n";
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
# sub get_keyword_value: get short keyword value for keyword lines. 
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
      #print "make_hkpdf.pl:**Error: no keyword value for $mnm in FITS lines-so will set value to <????> \n";
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
    open(OUTFILE, ">$directory$version_number/$filename") || die "3-Can't Open  $directory/$filename: $!\n" ;
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
    open(FILEAPID, "$s_arg[1]") || die "4-Can't Open $s_arg[1] file: $!\n";
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
  print "****************************************\n";
  print "** Note:Currently using option sort=4 **\n";
  print "****************************************\n";
  print "where apidfile is full path and filename containing list of apids to do\n";
  print "where apidlist contains list of apids to do(i.e.,0x1BD 0x081 ).\n";
  print "where apidlist contains keyword (i.e., \"ALL\", \"ALL_HK\", or \"ALL_SDO\")on apids to do\n";
  print "**********************************************************************\n";
  print "** Note:Currently using apidlist=\"ALL_HK\" or apidlist=\"ALL_SDO\"    **\n";
  print "**********************************************************************\n";
  if($help_arg eq "HELP")
  {
    print "(2)  Environment variable HK_STHA_DIRECTORY\n";
    print "(2a) Used to store location of STANFORD_TLM_HMI_AIA.txt file.\n";
    print "(3)  Environment variable HK_STHA_FILE\n";
    print "(3a) Used to store filename of STANFORD_TLM_HMI_AIA.txt file to process.\n";
    print "(3b) This file is input data for this script\n";
    print "(4)  Environment variable HK_CONFIG_DIRECTORY\n";
    print "(4a) Used to store location of hk config files for each apid\n";
    print "(4b) These files are the output for this script.\n";
    print "(4c) The filename format is apid-<apid#>-version-<version#>.\n";
    print "(4d) Or filename format is HMI-<apid#|name>-version-<version#>.\n";
    print "(4e) Or filename format is AIA-<apid#|name>-version-<version#>.\n";
    print "(5)  Environment variable HK_GTS_DIRECTORY\n";
    print "(5a) Used to store location of GODDARD_TLM_SDO.txt file.\n";
    print "(6)  Environment variable HK_GTS_FILE\n";
    print "(6a) Used to store filename of GODDARD_TLM_SDO.txt file to process.\n";
    print "(6b) This file is input data for this script\n";
    print "(7)  Environment variable HK_SH_CONFIG_DIRECTORY\n";
    print "(7a) Used to store location of sdo-hk config files for each apid\n";
    print "(7b) These files are the output for this script.\n";
    print "(7c) The filename format is apid-<apid#>-version-<version#>.\n";
    print "(8)  Examples:\n";
    print "(8a) make_hkpdf.pl sort=4 apidlist=\"0x032 0x1C0 0x1DE 0x21C 0x244\"\n";
    print "(8b) make_hkpdf.pl sort=4 apidlist=\"ALL_HK\"\n";
    print "(8c) make_hkpdf.pl sort=4 apidlist=\"0x032\"\n";
    print "(8d) make_hkpdf.pl sort=4 apidlist=\"ALL\"\n";
    print "(8e) make_hkpdf.pl sort=4 apidlist=\"ALL_SDO\"\n";

  }
  exit;
}

###############################################################################
# sub fix_timecode
###############################################################################
sub fix_timecode($,$,$)
{
  # pass arguments
  my $apid= $_[0];
  my $timecode= $_[1];
  my $shortlongflag= $_[2];

  # if 445 or 529 convert keyword strings
  if ($apid == 445 ||  $apid == 529)
  {
    $packetname="ISP";
    if($shortlongflag eq "LONG")
    {
      # change APID1BD_TIMECODE to APIDISP_TIMECODE
      @s_timecode=split("_",$timecode );
      $new_str= sprintf("%s%s_%s_%s", "APID", $packetname, $s_timecode[1],$s_timecode[2] );
    }
    else ##if else ($shortlongflag eq "SHORT")
    {
      # change HTCS1BD to HTCSISP or HTCSS1BD to HTCSSISP
      if ($timecode eq "HTCS1BD" || $timecode eq "ATCS211" )
      {
        $new_str= sprintf("%s%s", substr($timecode,0,4), $packetname );
      }
      else
      {
        $new_str= sprintf("%s%s", substr($timecode,0,5), $packetname );
      }
    }
  }
  elsif ($apid == 451 ||  $apid == 536)
  {
    $packetname="SEQ";
    if($shortlongflag eq "LONG")
    {
      @s_timecode=split("_",$timecode );
      $new_str= sprintf("%s%s_%s_%s", "APID", $packetname, $s_timecode[1],$s_timecode[2] );
    }
    else ##if else ($shortlongflag eq "SHORT")
    {
      # change HTCS1C3 to HTCSSEQ or HTCSS1C3 to HTCSSSEQ
      if ($timecode eq "HTCS1C3" || $timecode eq "ATCS02E")
      {
        $new_str= sprintf("%s%s", substr($timecode,0,4), $packetname );
      }
      else
      {
        $new_str= sprintf("%s%s", substr($timecode,0,5), $packetname );
      }
    }
  }
  elsif ($apid == 448 ||  $apid == 540)
  {
    $packetname="OBT";
    if($shortlongflag eq "LONG")
    {
      @s_timecode=split("_",$timecode );
      $new_str= sprintf("%s%s_%s_%s", "APID", $packetname, $s_timecode[1],$s_timecode[2] );
    }
    else ##if else ($shortlongflag eq "SHORT")
    {
      # change HTCS1C0 to HTCSOBT or HTCSS1C0 to HTCSSOBT
      if ($timecode eq "HTCS1C0" || $timecode eq "ATCS1C0")
      {
        $new_str= sprintf("%s%s", substr($timecode,0,4), $packetname );
      }
      else
      {
        $new_str= sprintf("%s%s", substr($timecode,0,5), $packetname );
      }
    }
  }
  elsif ($apid == 129)
  {
    $packetname="ASD";
    if($shortlongflag eq "LONG")
    {
      @s_timecode=split("_",$timecode );
      $new_str= sprintf("%s%s_%s_%s", "APID", $packetname, $s_timecode[1],$s_timecode[2] );
    }
    else ##if else ($shortlongflag eq "SHORT")
    {
      # change OTCS081 to OTCSASD or OTCSS081 to OTCSSASD
      if ($timecode eq "OTCS081")
      {
        $new_str= sprintf("%s%s", substr($timecode,0,4), $packetname );
      }
      else
      {
        $new_str= sprintf("%s%s", substr($timecode,0,5), $packetname );
      }
    }
  }
  return ($new_str);
}


###############################################################################
# sub fix_checksum_keywords
###############################################################################
sub fix_checksum_keywords($,$,$)
{
  # pass arguments
  my $apid= $_[0];
  my $checksum_kw= $_[1];
  my $shortlongflag= $_[2];

  # local variables 
  my $packetname;
  my $s_cs_kw;
  my $new_str;

  # if 445 or 529 convert keyword strings
  if ($apid == 445 ||  $apid == 529)
  {
    $packetname="ISP";
    if($shortlongflag eq "LONG")
    {
      # change APID1BD_CHECKSUM to APIDISP_CHECKSUM
      @s_cs_kw=split("_",$checksum_kw );
      $new_str= sprintf("%s%s_%s", "APID", $packetname, $s_cs_kw[1] );
    }
    else ##if else ($shortlongflag eq "SHORT")
    {
      # change HCSUM1BD to HCSUMISP
      $new_str= sprintf("%s%s", substr($checksum_kw,0,5), $packetname );
    }
  }
  elsif ($apid == 451 ||  $apid == 536)
  {
    $packetname="SEQ";
    if($shortlongflag eq "LONG")
    {
      @s_cs_kw=split("_",$checksum_kw );
      $new_str= sprintf("%s%s_%s", "APID", $packetname, $s_cs_kw[1] );
    }
    else ##if else ($shortlongflag eq "SHORT")
    {
      $new_str= sprintf("%s%s", substr($checksum_kw,0,5), $packetname );
    }
  }
  elsif ($apid == 448 ||  $apid == 540)
  {
    $packetname="OBT";
    if($shortlongflag eq "LONG")
    {
      @s_cs_kw=split("_",$checksum_kw );
      $new_str= sprintf("%s%s_%s", "APID", $packetname, $s_cs_kw[1] );
    }
    else ##if else ($shortlongflag eq "SHORT")
    {
      $new_str= sprintf("%s%s", substr($checksum_kw,0,5), $packetname );
    }
  }
  else
  {
     print "WARNING: This argument passed is not correct. Should not be calling this function.\n";
  }
  return ($new_str);
}

###############################################################################
# sub fix_apid_value
###############################################################################
sub fix_apid_value($,$)
{
  # pass arguments
  my $apid= $_[0];
  # local variables
  my $new_apid_value;

  # if 445 or 529 convert keyword strings
  if ($apid == 445 ||  $apid == 529)
  {
    $new_apid_value="ISP";
  }
  elsif ($apid == 451 ||  $apid == 536)
  {
    $new_apid_value="SEQ";
  }
  elsif ($apid == 448 ||  $apid == 540)
  {
    $new_apid_value="OBT";
  }
  else
  {
    print "WARNING:This case should not call this function.\n";
  }
  return ($new_apid_value);
}
###############################################################################
# sub check_for_merged_version
###############################################################################
sub check_for_merged_version($,$)
{
  # pass arguments
  my $ver= $_[0];
  my $aid= $_[1];

  #split file version number into whole number(WN) and decimal number(DN) for comparing below
  if (length($ver) > 4)
  {
    $format="%1.3f";
  }
  else
  {
    $format="%1.2f";
  }
  $str_ver= sprintf("$format",$ver);
  #print "$str_ver\n";
  @s_str_ver=split('\.', $str_ver);

  # if apid packet is KERNAL type which is has KER_VER_NUM value in keywords then do this
  # check using this constant if merge-series
  if (check_for_kernal_apid($aid) == KERNAL_APID)
  {
    if( $s_str_ver[1] >= KERNAL_MERGED_SERIES_FVN_DN)
    {
      return(MERGED_SERIES); #return(MERGE-FLAG)
    }
    else
    {
      return(NOT_MERGED_SERIES)
    }
  }
  else
  {

    # this apid packet is probably a HK HMI or HK AIA type so check using this merged constant value
    if(  $s_str_ver[1] >= HMI_AIA_MERGED_SERIES_FVN_DN )
    {
      return(MERGED_SERIES); #return(MERGE-FLAG)
    }
    else
    {
      return(NOT_MERGED_SERIES)
    }
  }
}
###############################################################################
# sub check_for_kernal_apid
###############################################################################
sub check_for_kernal_apid($)
{
  # pass arguments
  my $aid= $_[0];

  if($aid == 5 || $aid == 18 || $aid == 37 || $aid == 50 ||
     $aid == 448 || $aid == 478 || $aid == 540 || $aid == 580 )
  {
    return(KERNAL_APID); #return(MERGE-FLAG)
  }
  else
  {
    return(NOT_KERAL_APID)
  }
}
###############################################################################
# sub get_ground_column_value
###############################################################################
sub get_ground_column_value($)
{
  # pass arguments
  my $aid= $_[0];
  my $lookup_str=$column_table{$aid};

  if($lookup_str eq "HMI" || $lookup_str eq "AIA")
  {
    #use table to filter out KER packets and set column_id to KER
    return($lookup_str); 
  }
  elsif($lookup_str eq "AIAKER" || $lookup_str eq "HMIKER")
  {
    #use table to filter out KER packets and set column_id to KER
    if(($aid <= 31 && $aid > 0) || ($aid >= 400 && $aid < 500))
    {
      return("HMIKER");
    }
    elsif(($aid > 31 && $aid < 63) || ($aid >= 500 && $aid < 600))
    {
      return("AIAKER");
    }

    return($lookup_str); 
  }
  else
  {
    if(($aid <= 31 && $aid > 0) || ($aid >= 400 && $aid < 500))
    {
      return("HMI");
    }
    elsif(($aid > 31 && $aid < 63) || ($aid >= 500 && $aid < 600))
    {
      return("AIA");
    }
    elsif($aid >= 96 && $aid < 400)
    {
      return("SDO");
    }
    elsif($aid >= 2002 && $aid <= 2047)
    {
      return("SSIM");
    }
    else 
    {
      return("UNKN");
    }
  }
}
###############################################################################
# sub  check_bit_positions() 
###############################################################################
sub check_bit_positions($,@)
{

  # define local variables
  my(@keywordlines,@s_keywordlines,$pt); 

  # get variables passed: processing-type flag and array of keyword lines
  ($pt, @keywordlines) = @_; 

  # loop thru keyword lines and check if bit position is between 0-7 
  # note:LMSAL Users Guide defines bit position should be 0-7.
  $j=0;
  while( $keywordlines[$j] )
  {
    # split line at spaces and get bit position value
    @s_keywordlines=split(/ \s*/, $keywordlines[$j]);

    # check bit position is in range 0-7
    if ($s_keywordlines[4] > 7 || $s_keywordlines[4] < 0)
    {
      print "\nERROR: Found bit position in file not within range of 0-7 bits:<$s_keywordlines[4]>. ";
      print "Keyword with incorrect bit position are long name:<$s_keywordlines[2]> and short name:<$s_keywordlines[1]>. ";
      if($pt == SDO_PROCESSING)
      {
        print "Correct the keywords values in the input file<$ENV{'HK_GTS_DIRECTORY'}/$ENV{'HK_GTS_FILENAME'}> and then rerun. ";
      }
      elsif ($pt == HK_PROCESSING)
      {
        print "Correct the keywords values in the input file<$ENV{'HK_STHA_DIRECTORY'}/$ENV{'HK_STHA_FILENAME'} > and then rerun. ";
      }
      else
      {
        print "Correct the keywords values in the input file<Unknown> and then rerun. ";
      }
      print "Exiting make_hkpdf.pl script.\n";
      return -1;
    }
    $j++;
  }
  return 1;
}


###############################################################################
# sub  check_bit_positions() 
###############################################################################
sub check_patch01($,$)
{
   my($fvn, $id);
   $fvn=$_[0];
   $id=$_[1];
 #  print "fvn is <$fvn>\n";
 #  print "id is <$id>\n";
   @s_fvn=split('\.',$fvn);
   $fvn_wn= $s_fvn[0];
   $fvn_dn= $s_fvn[1];
   #print "fvn_dn is $fvn_dn\n";
   #print "fvn_wn is $fvn_wn\n";
   # patch warning for only these apids
   if ($id eq "0x02C" || $id eq "0x010")
   {
     if($s_fvn[0] == 1 && ($s_fvn[1] == 161 || $s_fvn[1] == 162 || $s_fvn[1] == 163 ))
     { 
       # if 1.161 to 1.163 then tell user to do patch
       return (1);
     }
     elsif($s_fvn[0] == 1 && $s_fvn[1] < 160 )
     { 
        # if 1.159 or less- tell user nothing- do NOT need patch
        return (0);
     }
     elsif($s_fvn[0] == 1 && $s_fvn[1] > 163 )
     { 
        # if greater then 1.163 file version then tell user may need to patch
        return (2);
     }
     elsif($s_fvn[0] > 1  )
     { 
       # if 2.x then tell user may need to patch
       return(3);
     }
   }
   return(0);
}

###############################################################################
# sub  displaynotice() 
###############################################################################
sub patch_notice($,$,$)
{
   my($msgnum,$vernum, $id);

   # get arguments passed in
   $msgnum=$_[0];
   $vernum=$_[1];
   $id=$_[2];

   if ($msgnum == 0)
   {
     return;
   }
   elsif($msgnum == 1)
   {
      print "\*\*\*Important Notification\*\*\*\n";
      print "-->For updates of HK configuration file with apid <$id> and file version <$vernum>, you need to run a patch 01 script.\n";
      print "-->Here is process steps. Run script. Check into cvs hk config files shown in log of patch01.pl script. Do cvs update of files on production.\n";
      print "-->This patch 01 fixes TRAC bug #261\n";
      print "-->Example run of script: patch01.pl  or patch01.pl > LOG-PATCH01\n";
   }
   elsif($msgnum == 2 || $msgnum == 3 )
   {
      print "\*\*\*Important Notification\*\*\*\n";
      print "-->For updates of HK configuration file versions <$vernum>, you \"MAY\" need to run a patch script.\n";
      print "-->If the STANFORD file was not fixed to make unsigned variables for Analog conversion keywords for packet 0x10 and 0x1c.\n";
      print "-->Then update patch01.pl to handle this new version number  <$vernum>.\n";
      print "-->Here is process steps. If needed update patch01.pl script. Run script.";
      print "Then check into cvs script and hk config files shown in log of patch01.pl script. Do cvs update of files on production.\n";
      print "-->This patch 01 fixes TRAC bug #261\n";
      print "-->Example run of script: patch01.pl  or patch01.pl > LOG-PATCH01\n";
   }
   else
   {
       print "WARNING: Probably could find  msgnum variable value. Check why the warning_flag was not set correctly in check_patch01().\n";
   }

}
