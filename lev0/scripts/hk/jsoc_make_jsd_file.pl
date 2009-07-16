#!/usr/bin/perl
##############################################################################
# Name:        jsoc_make_jsd_file.pl - Make JSD File                         #
# Description: Strip data from HKPDF files and create global and keyword     #
#              section for the JSOC Series Definition file. Uses GTCIDS file #
#              to get packet version number data. Use the <APID>-JSVN-TO-PVN #
#              map files to get jsoc series version number. Run first to make#
#              prelim jsd files. Then run a second time to make final JSDs.  #
#              The prelim files need to be there for all packet versions to  #
#              successfully create final jsd file.                           #
# Process:     (1)Run first to make prelim jsd files for needed APIDs.       #
#              (2)Run jsoc_do_jsvn_map_file.pl to make JSVN-TO-PVN map files #
#                 for each APID.                                             #
#              (3)Run make final jsd file to make the final files to use to  #
#                 create.                                                    #
# Execution:   (1)To run :                                                   #
#              jsoc_make_jsd_files.pl <prelim | final > <File-Version-Number>#
#                                     < processing type>                     #
#              (2)For help: jsoc_make_jsd_file.pl  -h                        #
# Example:     jsoc_make_jsd_file.pl prelim 1.79 <SDO | HK >                 #
# Example:     jsoc_make_jsd_file.pl final  1.79 <SDO | HK >                 #
# Author:      Carl                                                          #
# Date:        Move from EGSE to JSOC software environment on March 21, 2008 #
##############################################################################
# main function                                                              #    
##############################################################################
##get environment variables and initialize variables.
$hm=$ENV{'HOME'};
$scriptname="jsoc_make_jsd_file.pl";
$cvsver="1.5";

#PRODUCTION SETTINGS
$ENV{'HK_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_config_file";
$ENV{'HK_JSD_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_file/prod";
$ENV{'HK_JSD_PRELIM_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsd_prelim_file/prod";
$ENV{'HK_JSVN_MAP_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/hk_jsn_map_file/prod";
$ENV{'HK_SH_CONFIG_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/sdo_hk_config_file";
$ENV{'HK_SHCIDS_DIRECTORY'}="$hm/cvs/TBL_JSOC/lev0/sdo_hk_config_file/";
$ENV{'HK_GTCIDS_FILE'}="gtcids.txt";
$ENV{'HK_SHCIDS_FILE'}="shcids.txt";
$lnjsds=$ENV{'HK_LIST_NEW_JSDS'}="$hm/cvs/JSOC/proj/lev0/scripts/hk/jsoc_new_jsd_files.txt";

#setup Processing Type Constants
use constant NO_PROCESSING_TYPE =>  0;
use constant HK_PROCESSING_TYPE =>  1;
use constant SDO_PROCESSING_TYPE => 2;

# (1)check arguments used are correct and set processing type flag,  
#    prelim or final values flags and file version value.
($ffp,$ptf,$file_version_number)=&check_arguments();

# (2) make directory name to find hkpdf files
if ($ptf == SDO_PROCESSING_TYPE)
{
  $dncf = sprintf("%s/%s",  $ENV{'HK_SH_CONFIG_DIRECTORY'}, $file_version_number);
}
elsif ($ptf == HK_PROCESSING_TYPE)
{
  $dncf = sprintf("%s/%s",  $ENV{'HK_CONFIG_DIRECTORY'}, $file_version_number);
}
else
{
  print "ERROR:jsoc_make_jsd_file.pl: Bad processing type use:$ptf \n";
}

# (3) get all file in folder based on file version number based on filenames there process those apids.
(@fn_list)=&get_all_HKPDF_files($dncf);

# open New JSD File for holding new jsd files created
if($ffp == 1) 
{
  open(NJF, ">$lnjsds") || die "(1)Can't Open:<$lnjsds>: $!\n" ;
}

# (4) loop thru all HKPDF file and create a Prelim or Final JSD files for each apid in HKPDF file
# (4a)get list of apid to do
(@apidslist)=&get_list_apids_to_do($ptf, $dncf, $file_version_number);
print "@apidslist\n";

# (4b)start loop thru each HKPDF file checking which apids to process
foreach $file ( @fn_list)
{
  if ( $file eq "." || $file eq ".." || $file eq "CVS")
  {
     next;#skip doing
  }

  # check apid is one we want to create jsd for.These are packets with VER_NUM keywords
  # this skips doing prelim and final jsd that are not in list
  if(substr($file,0,5) eq "apid-" )
  {
    # check the old file format for HKPDF files is correct
    # get current apid
    $currapid=substr($file,5,3); 

    # check list of apid want to create JSDs for
    ($foundflg)=&check_apid_list($currapid);

    # if apid not in list skip creating jsd
    if ($foundflg == 0)
    {
      print "Warning:skipping $currapid because not in list of packets-probably because no VER_NUM field in packet.\n";
      next;#skip doing
    }
  }
  elsif(substr($file,0,4) eq "HMI-" || substr($file,0,4) eq "AIA-" || substr($file,0,4) eq "SDO-")
  {
    #got new file format for HKPDF filesnames- used for merged series name enhancement
    # get current apid
    $currapid=substr($file,4,3); 


    # handle packet names before checking in apidlist by converting to apid number
    if($currapid eq "ISP" && substr($file,0,4) eq "HMI-")
    {
      $currapid="1BD"; #445 decimal
    }
    elsif($currapid eq "ISP" && substr($file,0,4) eq "AIA-")
    {
      $currapid="211"; # 529 decimal
    }
    elsif($currapid eq "SEQ" && substr($file,0,4) eq "HMI-")
    {
      $currapid="1C3"; # 451 decimal
    }
    elsif($currapid eq "SEQ" && substr($file,0,4) eq "AIA-")
    {
      $currapid="218"; # 536 decimal
    }
    elsif($currapid eq "OBT" && substr($file,0,4) eq "HMI-")
    {
      $currapid="1C0"; #448 decimal
    }
    elsif($currapid eq "OBT" && substr($file,0,4) eq "AIA-")
    {
      $currapid="21C"; # 540 decimal
    }
    elsif($currapid eq "ASD" && substr($file,0,4) eq "SDO-")
    {
      $currapid="081"; # 540 decimal
    }

    # check list of apid want to create JSDs for
    ($foundflg) = &check_apid_list($currapid);

    # if apid not in list skip creating jsd
    if ($foundflg == 0)
    {
      next;#skip doing
    }
  }
  else
  {
      next;#skip doing
  }

# (5) get all lines in HKPDF hk or sdo configuration files for this apid
print "1--before get all lines:$dncf $file\n";
  ($apid_value)=&get_all_lines($dncf, $file);
print "apid returned for get all lines:$apid_value\n";

# (6) get packet version number or packet date based on file version number used for HK config files 
  if ($ptf == SDO_PROCESSING_TYPE)
  {
    # for sdo hk apids get packet date
    ($pkt_date,$found_pkt_date) = &get_pkt_date($currapid,$file_version_number);
  }
  elsif ($ptf == HK_PROCESSING_TYPE)
  {
    # for hk apids get packet version number
    ( $packet_version_number, $found_fvn) = &get_packet_version_number($file_version_number);

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


  }
  else
  {
    print "ERROR:jsoc_make_jsd_file.pl: Bad processing type use:$ptf \n";
  }


  if ($found_fvn == 0 && $ptf == HK_PROCESSING_TYPE)
  {
    # This occurs when there is a STHA file version but is not mapped in GTCIDS file
    # to a packet version number. This are no file versions 1.53,1.66,1.68,1.70.
    print "NOTE:For Apid <$apid_value>:Skipping creating preliminary jsd file.There is no packet version number ";
    print "found corresponding to file version number<$file_version_number> in GTCIDS file.\n";
    next;
  }
  elsif ($found_pkt_date == 0 && $ptf == SDO_PROCESSING_TYPE)
  {
    # This occurs when there is a GODDARD file version but is not mapped in SHCIDS file
    # to a packet date. This are no file versions like 1.1,1.2,1.3,1.4.
    print "NOTE:For Apid <$apid_value>:Skipping creating preliminary jsd file.There is no packet date that maps to file version number. ";
    print "Did not find file version number<$file_version_number> in SHCIDS file.\n";
    next;
  }

# (8) open <SDO|HMI|AIA><APID NUMBER|APID NAME>-JSVN-TO-PVN map files   #
#     and locate packet-version#(hk) or file-version(sdo) and get jsvn  #
  if( ($ffp == 1)  && (check_if_merged_pvn($ptf,$pvn_wn,$pvn_dn) || $ptf == SDO_PROCESSING_TYPE )) 
  {
    #got case of SDO APID or HK APID with merged packet version number-use new map filenames
    $jsvn = &get_jsvn($apid_value); 
  }
  else
  {
    # got case of HK APID with non-merged packet version number- use old map filenames
    &get_jsvn_from_old_files($apid_value );
  } 
#(6) create jsd file based on data in each HKPDF file
  &create_jsd_file($ptf, $jsvn); 
}

#close list in New JSD File 
if ($ffp == 1)
{
  close NJF;
}



###############################################################################
# sub get_jsvn: gets jsoc series version number from map files                #
###############################################################################
sub get_jsvn($)
{
  # setup my variables
  my($apid,$dnm,$mfn_list,$mapfile, $jsvn);
  $apid=$_[0];

  # geti directory name for map file to find jsvn map files
  $dnm = sprintf("%s",  $ENV{'HK_JSVN_MAP_DIRECTORY'});

  # open map file
  opendir(DIR_JSVN_MAP, $dnm) || die ":(2)Can't open directory:$dnm:$!\n"; #open subdirectory

  # read in all map files and put in list 
  @mfn_list = readdir(DIR_JSVN_MAP); #get a list of directory contents

  # initilize jsvn to some value to denote a problem getting jsvn
  $jsvn = "999";#no JSVN Map files exists

  # loop thru map files looking for match of apid to map file- if match-break from loop
  # note special handling for apid packet name in map filename
  # skip mapfiles that don't match the apid values
  foreach $mapfile  (@mfn_list)
  {
    if( hex substr($apid,2,4) == 445 &&  substr($mapfile,0,7) ne  "HMI-ISP")
    {
         next;#skip did not find map file
    }
    elsif( hex substr($apid,2,4) == 529 &&  substr($mapfile,0,7) ne  "AIA-ISP")
    {
         next;#skip did not find map file
    }
    elsif(hex substr($apid,2,4) == 451 &&  substr($mapfile,0,7) ne  "HMI-SEQ")
    {
         next;#skip did not find map file
    }
    elsif(hex substr($apid,2,4) == 536 &&  substr($mapfile,0,7) ne  "AIA-SEQ")
    {
         next;#skip did not find map file
    }
    elsif(hex substr($apid,2,4) == 448 &&  substr($mapfile,0,7) ne  "HMI-OBT")
    {
         next;#skip did not find map file
    }
    elsif(hex substr($apid,2,4) == 540 &&  substr($mapfile,0,7) ne  "AIA-OBT")
    {
         next;#skip did not find map file
    }
    elsif(hex substr($apid,2,4) == 129 &&  substr($mapfile,0,7) ne  "SDO-ASD")
    {
         next;#skip did not find map file
    }
    elsif(int substr($mapfile,4,3) ne  hex substr($apid,2,4) && !check_for_merged_apid(hex $apid))
    {
       next;#if non-merged apid for map and don't match apid then skip
    }

    # got match of map file to apid- so get map file to lookup jsvn
    $jsvn = &get_mapfile_jsvn($mapfile);
    last; #break from loop-found jsvn value!
  }#foreach
  closedir DIR_JSVN_MAP; #close directory
  return($jsvn);
}


  
################################################################################
# sub get_file_version_number: gets file version number based on packet version#
################################################################################
sub get_packet_version_number($)
{
  # setup my variables
  my( $dnm, $dn_fn, $fn, @all_gtcids_lines, $pvn, $foundfvn, $fvnArg);

  # arguments passed
  $fvnArg=$_[0];

  # get directory name to find jsvn map files
  $dnm = sprintf("%s",  $ENV{'HK_JSVN_MAP_DIRECTORY'});

  # open GTCIDS file to look for packet version number
  $dn_fn = sprintf("%s/%s",  $ENV{'HK_CONFIG_DIRECTORY'}, $ENV{'HK_GTCIDS_FILE'});
  open(FILE, "$dn_fn") || die "(4)Can't Open $dn_fn: $!\n";
  $foundfvn = 0;
  while (<FILE>)
  {
      @s_line= split / \s*/,  $_;
      if ( substr($s_line[0],0,1) eq '#')
      {
        next; #skip comment lines
      }
      else
      {
        #if hmi get pvn from field 7
        if( $sys_value eq "hmi" ) 
        {
          if ( $s_line[13] eq $fvnArg)
          {
            $pvn = $s_line[7];
            $stha_ci_date = $s_line[0];
            $stha_ci_time = $s_line[1];
            $stha_gen_time = $s_line[15];
            $foundfvn = 1;
            last;
          }
        }
        elsif( $sys_value eq "aia" ) 
        {
          #if hmi get pvn from field 9
          if ( $s_line[13] eq $fvnArg)
          {
            $pvn = $s_line[9];
            $stha_ci_date = $s_line[0];
            $stha_ci_time = $s_line[1];
            $stha_gen_time = $s_line[15];
            $foundfvn = 1;
            last;
          }
        }
        elsif( $sys_value eq "aiaker" || $sys_value eq "hmiker") 
        {
          #if hmi get ker from field 5
          if ( $s_line[13] eq $fvnArg)
          {
            $pvn = $s_line[5];
            $stha_ci_date = $s_line[0];
            $stha_ci_time = $s_line[1];
            $stha_gen_time = $s_line[15];
            $foundfvn = 1;
            last;
          }
        }
        elsif( $sys_value eq "ssim" ) 
        {
          if ( $s_line[13] eq $fvnArg)
          {
            $pvn = $s_line[7]; #use hmi column as default value
            $stha_ci_date = $s_line[0];
            $stha_ci_time = $s_line[1];
            $stha_gen_time = $s_line[15];
            $foundfvn = 1;
            last;
          }
        }
        else
        {
           print "Warning:jsoc_make_jsd_file.pl: $sys_value is not hmi, aia, sdo, ker or ssim. Check setting.\n";
           $foundfvn = 0;
        }
      }#big else
  }
  close(FILE);
  # return Packet Version number found and set foundfvn flag to 1 if found
  return ( $pvn, $foundfvn);
}



###############################################################################
# sub get_all_hkpdf_files: gets all files in folder specified by file version#
###############################################################################
sub get_all_HKPDF_files($)
{
  # declare local variables
  my (@fnl, $dn_arg);

  # get directory name argument passed
  $dn_arg=$_[0];

  # open directory and read in all HKPDF file
  opendir(DIR, $dn_arg) || die "(5)Can't open:$!\n"; #open subdirectory
  @fnl = readdir(DIR); #get a list of directory contents
  closedir(DIR); #close directory
  
  # return list of all HKPDF files - apid-#-version-# files
  return (@fnl);
} 



###############################################################################
# sub get_all_lines: gets all lines from HKPDF file
###############################################################################
sub get_all_lines($$)
{
    
    my($i,$cdn,$apid,$cfn);
    $cdn = $_[0];   #directory for HKPDF files
    $cfn = $_[1];   #HKPDF config filename
    @all_kwd_lines="";
    @all_dcon_lines="";
    @all_acon_lines="";

    # skip over processing . or CVS, etc files
    if( (substr($cfn,0,4) ne "apid" && substr($cfn,9,7) ne "version"   ) &&
        (substr($cfn,0,3) ne  "HMI"  && substr($cfn,9,7) ne "version"  ) &&
        (substr($cfn,0,3) ne  "AIA"  && substr($cfn,9,7) ne "version"  ) &&
        (substr($cfn,0,3) ne  "SDO"  && substr($cfn,9,7) ne "version"  )  )
    {
        return;
    }

    # open one HKPDF (apid-#-version-# or HMI-#-version-#) file
    open(FILE, "$cdn/$cfn") || die "(6)Can't Open $cdn/$cfn file: $!\n";
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
        # set value from APID line in HKPDF file
        $apid_lines= $_ ;
        $apid = $s_line[1] ;
        $sys_value = lc $s_line[3] ;

        # fix apid value for merged HKPDF files 
        if ($apid eq "ISP" && $s_line[3] eq "AIA")
        {
           $apid = "0x211"; #dec value 529
        }
        elsif ($apid eq "ISP" && $s_line[3] eq "HMI")
        {
           $apid = "0x1BD"; #decimal value 445
        }
        elsif ($apid eq "SEQ" && $s_line[3] eq "HMI")
        {
           $apid = "0x1C3"; #dec value 451
        }
        elsif ($apid eq "SEQ" && $s_line[3] eq "AIA")
        {
           $apid = "0x218"; #dec value 536
        }
        elsif ($apid eq "OBT" && $s_line[3] eq "HMIKER")
        {
           $apid = "0x1C0"; #dec value 448
        }
        elsif ($apid eq "OBT" && $s_line[3] eq "AIAKER")
        {
           $apid = "0x21C"; #dec value 540
        }
        else
        {
          ;#skip fixing
        }

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
    return ($apid);
}



###############################################################################
# sub create_jsd: create jsd based on data in HKPDF file
###############################################################################
sub create_jsd_file($$)
{
  # declare local variable and set passed arguments
  my ($p_ptf, $jsvn, $default_value);
  $p_ptf=$_[0];#processing type flag
  $jsvn=$_[1];

  if ($ffp == 1 )
  {
    if ($jsvn eq "999")
    {
       #print "skip doing jsn is 999\n";
       return;
    }
  }

  # use environment to set or use arguments(-n<namespace> -i<identifer> ) to set
  if( $sys_value eq "hmi" || $sys_value eq "aia" ||  $sys_value eq "ssim" ||  $sys_value eq "sdo") 
  {
    $jsd_namespace = sprintf("%s%s", $sys_value,"");
  }
  elsif($sys_value eq "aiaker" || $sys_value eq "hmiker")
  {
    $jsd_namespace = sprintf("%3.3s%s",  $sys_value,"");
  }
  else
  {
    $jsd_namespace = sprintf("%s%s", $sys_value,"");
  }

  $jsd_identifier = "lev0";
  $author="production";
  $owner="production";
  $unitsize="1";
  $archive= "0";
  $index="PACKET_TIME";

  # set tape ground based on if hmi,aia or ssim 
  if ( $sys_value eq "hmi")
  {
    $tapegroup= "0";
    $retention= "60";
  }
  elsif ( $sys_value eq "aia")
  {
    $tapegroup= "0";
    $retention= "60";
  }
  elsif ( $sys_value eq "sdo")
  {
    $tapegroup= "0";
    $retention= "60";
  }
  elsif ( $sys_value eq "aiaker" || $sys_value eq "hmiker")
  {
    $tapegroup= "0";
    $retention= "60";
  }
  elsif ( $sys_value eq "ssim")
  {
    $tapegroup= "0";
    $retention= "60";
  }
  elsif ( $sys_value eq "unkn")
  {
    $tapegroup= "0";
    $retention= "0";
  }
  else
  {
    print "-Warning:jsoc_make_jsd_file.pl: Do not know sys value <$sys_value>\n";
  }

  if ($p_ptf == HK_PROCESSING_TYPE)
  {

    # format packet version number in MMMmmm decimal format which is used for filename
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
      print "**WARNING:jsoc_make_jsd_file.pl: The Packet Version number for jsd file name is not set correctly.\n";
      print "Check code and updating perl script is probably needed.\n";
    }
  }

  if ($ffp == 1)
  {
    # check if merged series

    if($s = &check_for_merged_apid(hex  substr($apid_value,2,4)) && (&check_if_merged_pvn($p_ptf,$pvn_wn,$pvn_dn) || &check_if_merged_fvn($p_ptf,$file_version_number)))
    {
      # got HK APID with merged packet version number(1.194) or SDO APID
      # get apid "string" for merged series
      $pn= get_series_packet_name(hex substr($apid_value,2,4));

      # set apid value
      $jsd_fn= sprintf("%s.%s_%s_%-04.4d.jsd",$jsd_namespace,$jsd_identifier,lc $pn, $jsvn);
      $jsd_sn= sprintf("%s.%s_%s_%-04.4d",$jsd_namespace,$jsd_identifier, lc $pn, $jsvn);
    }
    else
    {
      # create file name for final jsd files using apid "number"
      $jsd_fn= sprintf("%s.%s_%-04.4d_%-04.4d.jsd",$jsd_namespace,$jsd_identifier,hex  substr($apid_value,2,4), $jsvn);
      $jsd_sn= sprintf("%s.%s_%-04.4d_%-04.4d",$jsd_namespace,$jsd_identifier,hex  substr($apid_value,2,4), $jsvn);
    }

    # get directory to write jsd file
    $directory=$ENV{'HK_JSD_DIRECTORY'};
  }
  else ## ffp==0
  {
    if ($p_ptf == SDO_PROCESSING_TYPE)
    {
      $jsd_fn= sprintf("%s.%s_%-04.4d_%-8.8d",$jsd_namespace,$jsd_identifier,hex  substr($apid_value,2,4), $pkt_date);
    }
    elsif ($p_ptf == HK_PROCESSING_TYPE)
    {
      # create file name for preliminary jsd files used for creating jsn map files.
print  substr($apid_value,2,4);
print "\n";
      $jsd_fn= sprintf("%s.%s_%-04.4d_%-03.3d%-03.3d",$jsd_namespace,$jsd_identifier,hex  substr($apid_value,2,4), $pvn_wn,$pvn_dn);
    }
    # get directory to write jsd file
    $directory=$ENV{'HK_JSD_PRELIM_DIRECTORY'};
  }

  # check if version_number folder exists
  opendir(DIR,$directory) || die "(7)Can't open:$directory  $!\n";

  # check if file exists and put out warning message if does
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
    
    printf(OUTFILE "#============================================================================\n");
    printf(OUTFILE "# Created by:$scriptname  cvs version:$cvsver\n");
    printf(OUTFILE "#============================================================================\n");
    printf(OUTFILE "#======================== Global Series Information =========================\n");
    printf(OUTFILE  "SeriesName:  %-s\n", $jsd_sn);
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
  $print_value="2";
  $default_value = "TSEQ_EPOCH";
  if ($p_ptf == SDO_PROCESSING_TYPE)
  {
    $unit_value="UTC";
  }
  else
  {
    $unit_value="TAI";
  }
  $comment_value = "\"PACKET_TIME\"";

  # Setting based on manually updated latest values in hmi.lev0_60d.jsd
  $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
              "Keyword:", $keyword_name, $type_value,
              "ts_eq", "record", $default_value, 
              $print_value, $unit_value, $comment_value); 
  printf(OUTFILE "%s", $tmp_line);

  # Slotted keywords for PACKET TIME epoch
  $keyword_name = "PACKET_TIME_epoch";
  $type_value = "time";
  $print_value="2";
  $default_value = "TSEQ_EPOCH";
  $unit_value="TAI";
  $comment_value = "\"PACKET_TIME_epoch\"";

  # Setting based on manually updated latest values in hmi.lev0_60d.jsd
  $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
              "Keyword:", $keyword_name, $type_value,
              "constant", "record", $default_value, 
              $print_value, $unit_value, $comment_value); 
  printf(OUTFILE "%s", $tmp_line);

  # Slotted keywords for PACKET TIME step
  $keyword_name = "PACKET_TIME_step";
  $type_value = "float";
  $print_value="%f";
  $default_value = "0.1";
  #if ( &check_if_merged_pvn($p_ptf,$pvn_wn,$pvn_dn) || &check_if_merged_fvn($p_ptf,$file_version_number))
  if ( &check_if_merged_pvn($p_ptf,$pvn_wn,$pvn_dn) || ($ffp == 1))
  {
    $unit_value="\"0.1 sec\"";
  }
  elsif ( &check_if_merged_fvn($p_ptf,$file_version_number))
  {
    # SDO_PROCESSING 
    $unit_value="\"0.1 sec\"";
  }
  else
  {
    #keep old way for older versions less than packet version number 1.196
    $unit_value="0.1 sec";
  }
  $comment_value = "\"PACKET_TIME_step\"";

  # Setting based on manually updated latest values in hmi.lev0_60d.jsd
  $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
              "Keyword:", $keyword_name, $type_value,
              "constant", "record", $default_value, 
              $print_value, $unit_value, $comment_value); 
  printf(OUTFILE "%s", $tmp_line);


  if ($p_ptf == HK_PROCESSING_TYPE)
  {
    # create and print to out file PACKET VERSION NUMBER keyword
    $keyword_name = "PACKET_VERSION_NUMBER";
    $type_value = "string";
    $print_value="%s";
    $default_value ="DRMS_MISSING_VALUE";
    $unit_value="none";
    $comment_value = "\"PACKET_VERSION_NUMBER, decimal value with format:MMM.mmm\"";
    $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
                "Keyword:", $keyword_name, $type_value,
                "variable", "record", $default_value, 
                $print_value, $unit_value, $comment_value); 
    printf(OUTFILE "%s", $tmp_line);
  }
  elsif ($p_ptf == SDO_PROCESSING_TYPE)
  {
    # create and print to out file FILE VERSION NUMBER keyword
    $keyword_name = "FILE_VERSION_NUMBER";
    $type_value = "string";
    $print_value="%s";
    $default_value ="DRMS_MISSING_VALUE";
    $unit_value="none";
    $comment_value = "\"FILE_VERSION_NUMBER, decimal value with format:MMM.mmm\"";
    $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
                "Keyword:", $keyword_name, $type_value,
                "variable", "record", $default_value, 
                $print_value, $unit_value, $comment_value); 
    printf(OUTFILE "%s", $tmp_line);
  }
  else
  {
    print "ERROR:jsoc_make_jsd_file.pl: got unkown processing type:<$p_ptf>\n";
    exit;
  }

  # check if packet version number is for merged series
  if ( &check_if_merged_pvn($p_ptf,$pvn_wn,$pvn_dn) || &check_if_merged_fvn($p_ptf,$file_version_number))
  {
    # create and print to out file HK SOURCE keyword
    $keyword_name = "HK_SOURCE";
    $type_value = "string";
    $print_value="%s";
    $default_value ="DRMS_MISSING_VALUE";
    $unit_value="none";
    $comment_value = "\"HK_SOURCE, string format:<dayfile name>[<packet-data-value>][<packet-apid>][<source>]\"";
    $tmp_line = sprintf("%-s%-s,%-s,%-s,%-s,%-s,%-s,%-s,%-s\n",
                "Keyword:", $keyword_name, $type_value,
                "variable", "record", $default_value, 
                $print_value, $unit_value, $comment_value); 
    printf(OUTFILE "%s", $tmp_line);
  }

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
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[6] eq "SB" && $s_line[7] eq "R") 
    {
      $type_value = "short";
      $print_value=" %hhd";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[6] eq "IU1" && $s_line[7] eq "R") 
    {
      $type_value = "int";
      $print_value=" %d";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[6] eq "IS1" && $s_line[7] eq "R") 
    {
      $type_value = "short";
      $print_value=" %hd";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[6] eq "IL1" && $s_line[7] eq "R") 
    {
      $type_value = "int";
      $print_value=" %d";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[6] eq "UL1" && $s_line[7] eq "R") 
    {
      if ($p_ptf == SDO_PROCESSING_TYPE)
      {
        $print_value=" %lu";
      }
      else
      {
        $print_value=" %ld";
      }
      $type_value = "longlong";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[6] eq "DFP" && $s_line[7] eq "R") 
    {
      $type_value = "double";
      $print_value=" %f";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[6] eq "SFP" && $s_line[7] eq "R") 
    {
      $type_value = "float";
      $print_value=" %f";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[7] eq "D") 
    {
      $type_value = "string";
      $print_value=" %s";
      $default_value ="DRMS_MISSING_VALUE";
    }
    elsif($s_line[7] eq "A") 
    {
      $type_value = "double";
      $print_value=" %lf";
      $default_value ="DRMS_MISSING_VALUE";
    }  
    else
    {
      $type_value = "int";
      $print_value=" %d";
      $default_value ="DRMS_MISSING_VALUE";
      print "**WARNING: Apid: $apid_value Keyword: $s_line[1] type value not found. NOT VALID TYPE! Skipping using keyword.\n";
      next;#skip writing to jsd file
    }

    # use hardcoded value for units
    $unit_value="none";

    # use long name in comments
    $comment_value = sprintf("\"%s\"", $keyword_name);

    # if apid is merged apid series which includes setting TIMECODE and CHECKSUM values with packet name instead of hex number
    if (&check_for_merged_apid(hex  substr($apid_value,2,4)) && &check_for_timecode_checksum($s_line[2]) &&  (&check_if_merged_pvn($p_ptf,$pvn_wn,$pvn_dn) || &check_if_merged_fvn($p_ptf,$file_version_number)))
    {
      ($kwl,$kws) = create_merged_keyword_name($s_line[2], $comment_value, hex substr($apid_value,2,4));
      $s_line[2]=$kwl;
      $comment_value=$kws; 
    }
  
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



###########################################
#  check apid list: check apid in list   #
###########################################
sub check_apid_list ($)
{
  my(@alist,$capid,$ff, $item);
  # get arguments passed
  $capid = $_[0];
  @alist=@apidslist; #set alist to global array apidslist
#print "alist is  @alist\n";

  $ff=0;
  foreach $item (@alist)
  {
    if($capid =~ m/$item/i && $capid ne "")
    {
      $ff = 1;
      last;
    }
    else
    {
      $ff=0;
    }
  }
 return ($ff);
}

###########################################
#  get list of apids to do                #
###########################################
sub get_list_apids_to_do($$$)
{
  # declare locals
  my ($dn, $files, @al,$process_type_flg, $fvn);

  # get arguments passed
  $process_type_flg=$_[0];
  $dn=$_[1];
  $fvn=$_[2];

  if ($process_type_flg == SDO_PROCESSING_TYPE)
  {
      # get files in sdo hk config directory
      $sdofiles=`cd $dn;ls SDO* `;

      # remove HMI text regular expression
      $sdofiles =~ s/(SDO.?-)//g;
     
      # remove everything from after apid number to end of line regular expression
      $sdofiles =~ s/(-.+?\n)/ /g;
    
      # fix for dealing with apid list -substitute ISP for hmi with 1BD and SEQ for hmi with 1E1
      $sdofiles =~ s/ASD/081/g;

      #put sdo list in $files
      $files=sprintf("%s", $sdofiles);
  }
  elsif ($process_type_flg == HK_PROCESSING_TYPE)
  {
    # misses VER_TEMPERATURE
    # $files=`cd /home1/carl/cvs/TBL_JSOC/lev0/hk_config_file/$file_version/; grep VER_NUM apid*`;
    # get all apids that have VER_NUM or VER_TEMPERATURE Keyword
    @s_fvn=split('\.', $fvn );
    if (int $s_fvn[0] == 1  && int $s_fvn[1] >= 160 ) 
    {
      # get hmi files with VER
      $hmifiles =`cd $dn;  egrep '(VER_NUM|VER_TEMPERATURE)' HMI-*`;
      
      # remove HMI text regular expression
      $hmifiles =~ s/(H.+?-)//g;
     
      # remove everything from after apid number to end of line regular expression
      $hmifiles =~ s/(-.+?\n)/ /g;
    
      # fix for dealing with apid list -substitute ISP for hmi with 1BD and SEQ for hmi with 1E1
      $hmifiles =~ s/ISP/1BD/g;
      $hmifiles =~ s/SEQ/1C3/g;
      $hmifiles =~ s/OBT/1C0/g;
   
      # get aia files with VER
      $aiafiles =`cd $dn;  egrep '(VER_NUM|VER_TEMPERATURE)' AIA-*`;

      # remove AIA text regular expression
      $aiafiles =~ s/(AIA-)//g;

      # remove everything from after apid number to end of line regular expression
      $aiafiles =~ s/(-.+?\n)/ /g;

      # fix for dealing with apid list -substitute ISP for hmi with 1BD and SEQ for hmi with 1E1
      $aiafiles =~ s/ISP/211/g;
      $aiafiles =~ s/SEQ/218/g;
      $aiafiles =~ s/OBT/21C/g;

      #get oldfiles format for KER type files-this case occurs only for fvn 1.160 where formats are mixed.
      $oldfiles=`cd $dn;  egrep '(VER_NUM|VER_TEMPERATURE)' apid*`;
      # remove apid text regular expression
      $oldfiles =~ s/(a.+?-)//g;
      # remove everything from after apid number to end of line regular expression
      $oldfiles =~ s/(-.+?\n)/ /g;
      print "oldfiles forgot= $oldfiles\n";

      #put both aia and hmi lists together
      $files=sprintf("%s%s%s", $hmifiles,$aiafiles, $oldfiles);
      print "files = <$files>\n";
    }
    else
    {
      #handle old file formats for HKPDF files--apid-1BD-version-1.159
      $files=`cd $dn;  egrep '(VER_NUM|VER_TEMPERATURE)' apid*`;
      # remove apid text regular expression
      $files =~ s/(a.+?-)//g;
      # remove everything from after apid number to end of line regular expression
      $files =~ s/(-.+?\n)/ /g;
    }
  }

  # split apids into separate field in array
  @al = split(/ /, $files);

  # return list of apids to do
  return ( @al);
}



####################################################
#  check arguments
####################################################
sub check_arguments()
{
  my ($ffp, $pt_type, $fvn);
  if ( ($#ARGV != 2) && ($#ARGV == 0 && $ARGV[0] ne "-h")) 
  {
    die "Usage: $0 < prelim | final | -h > < File Version Number > < HK | SDO >:$!";
  }
  elsif  ($#ARGV == 1 || $#ARGV > 2)
  {
    die "Usage: $0 < prelim | final | -h > < File Version Number > < HK | SDO >:$!";
  }

  if ( $ARGV[0] eq "-h" )
  {
    print "Help Listing\n";
    print "(1)Ways to Execute Perl Script: \n";
    print "(1a)Create Preliminary JSD Files: jsoc_make_jsd_file.pl prelim <File Version Number> <processing type>\n";
    print "(1b)Create Final JSD Files:       jsoc_make_jsd_file.pl final  <File Version Number> <processing type>\n";
    print "    Where file type can be prelim or final.\n";
    print "    Where file version number for HK can be 1.32 to 1.255 and for SDO can be 1.1 to 1.255\n";
    print "    Where processing type is SDO or HK.\n";
    print "(1c)Get Help Information:         jsoc_make_jsd_file.pl -h \n";
    print "(2) Environment variable HK_JSD_PRELIM_DIRECTORY and HK_JSD_DIRECTORY\n"; 
    print "(2a)Used to control where to put prelim or final JSD files\n";
    print "(3) Environment variable HK_GTCIDS_FILE\n";
    print "(3a)Used to store filename of gtcids.txt file which is use to look up packet version number.\n";
    print "(4) Environment variable HK_SHCIDS_FILE\n";
    print "(4a)Used to store filename of shcids.txt file which is use to look up packet date based on file version number.\n";
    print "(5) Environment variable HK_CONFIG_DIRECTORY and HK_SH_CONFIG_DIRECTORY\n";
    print "(5a)Used to store location of hk and sdo config (HKPDF) files for each apid\n";
    print "(5b)The filename format is apid-<apid#>-version-<version#>.\n";
    print "(5c)The HKPDF files are used to get print values, keyword and long names, etc >.\n";
    exit;
  }
  elsif ( $ARGV[0] eq "prelim")
  {
    $fvn = $ARGV[1];
    @s_fvn=split('\.', $fvn );
    $ffp= 0; #flag for final or prelim value;
  }
  elsif ( $ARGV[0] eq "final")
  {
    $fvn = $ARGV[1];
    @s_fvn=split('\.', $fvn );
    $ffp= 1;
  }
  else
  {
    print "ERROR:jsoc_make_jsd_file.pl: use either prelim or final for first argument(i.e.,prelim)\n";
    print "ERROR:jsoc_make_jsd_file.pl: use packet version number for second argument(i.e.,1.132)\n";
  }
  if ( $ARGV[2] eq "SDO" )
  {
    $pt_flag=SDO_PROCESSING_TYPE;
  }
  elsif ( $ARGV[2] eq "HK" )
  {
    $pt_flag=HK_PROCESSING_TYPE;
  }
  else
  {
    $pt_flag=NO_PROCESSING_TYPE;
    print "ERROR:jsoc_make_jsd_file.pl:Bad processing type parameter for arguments:$ARGV[2]\n";
    die "Usage: $0 < prelim | final | -h > < File Version Number > < HK | SDO >:$!";
  }
  return($ffp,$pt_flag,$fvn);
}



################################################################################
# sub get_pkt_date: gets packet date based on file version and apid            #
################################################################################
sub get_pkt_date($$)
{
  my( $dn_fn,$found_pd,$p_apid,$p_fvn,$p_date);
  $p_apid=sprintf("%-04.4d",hex $_[0]);
  $p_fvn=$_[1];

  # open SHCIDS file to look for packet version number
  $dn_fn = sprintf("%s/%s",  $ENV{'HK_SH_CONFIG_DIRECTORY'}, $ENV{'HK_SHCIDS_FILE'});
  open(FILE, "$dn_fn") || die "(6)Can't Open $dn_fn: $!\n";
  $found_pd = 0;
  while (<FILE>)
  {
    #push(@all_gtcids_lines, $_) ;
    @s_line= split / \s*/,  $_;
    if ( substr($s_line[0],0,1) eq '#')
    {
      next; #skip
    }
    else
    {
      if ( $s_line[5] eq $p_fvn and  $s_line[3] eq $p_apid)
      {
        $p_date = $s_line[0];
        $found_pd = 1;
        #reformat to yyyymmdd
        $p_date =~ s/\///g; #regular exp rm /
        break;
      }
    }
  }
  close(FILE);

  if ($found_pd != 1)
  {
    print "ERROR:jsoc_make_jsd_file.pl: Could not find apid and file version number in shcids.txt file. Exiting...\n";
    print " Need to add this line to file to process JSD files for apid:<$p_apid> file version number:<$p_fvn>\n";
    exit();
  }
  return ($p_date,$found_pd);
}

################################################################################
# sub check_for_merged_apid: check if the passed apid is a merged series apid  #
################################################################################
sub check_for_merged_apid($)
{
  my( $apid);
  $apid=$_[0];

  if ($apid == 445 || $apid == 475 || $apid == 29)
  {
     return (1);#HMI ISP's!
  }
  elsif($apid == 529 || $apid == 569 || $apid == 39)
  {
     return (1);#AIA ISP's!
  }
  elsif($apid == 481 || $apid == 451 || $apid == 21)
  {
     return (1);#HMI SEQ's!
  }
  elsif($apid == 536 || $apid == 576 || $apid == 46)
  {
     return (1);#AIA SEQ's!
  }
  elsif($apid == 448 || $apid == 478 || $apid == 18)
  {
     return (1);#HMI OBT's!
  }
  elsif($apid == 540 || $apid == 580 || $apid == 50)
  {
     return (1);#AIA OBT's!
  }
  elsif($apid == 129 )
  {
     #return (1);#SDO-ASD's!
     return (1);
  }
  else
  {
     return (0);
  }

}
################################################################################
# sub check_for_timecode_checksum: check if timecode or checksum keyword       #
################################################################################
sub check_for_timecode_checksum($)
{
  my( $keyword);
  $keyword=$_[0];

  if($keyword =~ m/(TIMECODE)/ && $keyword =~ m/(APID)/)
  {
    return (1);
  }
  elsif($keyword =~ m/(CHECKSUM)/ && $keyword =~ m/(APID)/)
  {
    return (1);
  }
  else 
  {
    return (0);
  }
} 


################################################################################
# sub create_merged_keyword_name                                               #
# Create timecodes and checksum names for merged series.                       #
################################################################################
sub create_merged_keyword_name($$$)
{
  #pass arguments
  my( $keyword, $apid);
  $keywordL =$_[0];
  $keywordS =$_[1];
  $apid=$_[2];

  #local variables
  my($newstrS, $newstrL);

  #get packet name to use(i.e., ISP,SEQ,OBT,ASD, etc)
  $packetname= &get_series_packet_name($apid);

  #fix timecode or checksum short and long keyword names
  if($keywordL =~ m/(TIMECODE)/ && $keywordL =~ m/(SECONDS)/)
  { 
    #fix long and short timecode seconds keyword names
    $newstrL = sprintf("%s%s%s",substr($keywordL,0,4),$packetname,substr($keywordL,7)); 
    $newstrS = sprintf("\"%s%s\"",substr($keywordS,1,4),$packetname); 
    return ($newstrL, $newstrS );
  }
  elsif($keywordL =~ m/(TIMECODE)/ && $keywordL =~ m/(SUBSECS)/)
  {
    $newstrL = sprintf("%s%s%s",substr($keywordL,0,4),$packetname,substr($keywordL,7)); 
    $newstrS = sprintf("\"%s%s\"",substr($keywordS,1,5),$packetname); 
    return ($newstrL, $newstrS );
  }
  elsif($keywordL =~ m/(CHECKSUM)/)
  {
    $newstrL = sprintf("%s%s%s",substr($keywordL,0,4),$packetname,substr($keywordL,7)); 
    $newstrS = sprintf("\"%s%s\"",substr($keywordS,1,5),$packetname); 
    return ($newstrL, $newstrS );
  }
  else
  {
    print "jsoc_make_jsd_file.pl:WARNING: should not get into this cased.\n";
    return("","");
  }
}



################################################################################
# sub check_for_merged_apid: check if the passed apid is a merged series apid  #
################################################################################
sub get_series_packet_name($)
{
  my( $apid);
  $apid=$_[0];

  if ($apid == 445 || $apid == 475 || $apid == 29)
  {
     return ("ISP");#HMI ISP's!
  }
  elsif($apid == 529 || $apid == 569 || $apid == 39)
  {
     return ("ISP");#AIA ISP's!
  }
  elsif($apid == 481 || $apid == 451 || $apid == 21)
  {
     return ("SEQ");#HMI SEQ's!
  }
  elsif($apid == 536 || $apid == 576 || $apid == 46)
  {
     return ("SEQ");#AIA SEQ's!
  }
  elsif ($apid == 448 || $apid == 478 || $apid == 18)
  {
     return ("OBT");#HMI OBT's!
  }
  elsif($apid == 540 || $apid == 580 || $apid == 50)
  {
     return ("OBT");#AIA ISP's!
  }
  elsif($apid == 129 )
  {
     return ("ASD");#SDO-ASD's!
  }
  else
  {
     return ("");
  }
}

################################################################################
# sub check_for_merged_apid: check if the passed apid is a merged series apid  #
################################################################################
sub check_for_merged_apid_in_mapfile($$)
{
  my( $apid, $mp);
  $apid=$_[0];
  $mf=$_[1];

  if(hex substr($apid,2,4) == 445 &&  substr($mf,4,3) eq  "ISP")
  {
    return (1);
  }
  elsif(hex substr($apid,2,4) == 529 &&  substr($mf,4,3) eq  "ISP")
  {
    return (1);
  }
  elsif(hex substr($apid,2,4) == 451 &&  substr($mf,4,3) eq  "SEQ")
  {
    return (1);
  }
  elsif(hex substr($apid,2,4) == 536 &&  substr($mf,4,3) eq  "SEQ")
  {
    return (1);
  }
  elsif(hex substr($apid,2,4) == 448 &&  substr($mf,4,3) eq  "OBT")
  {
    return (1);
  }
  elsif(hex substr($apid,2,4) == 540 &&  substr($mf,4,3) eq  "OBT")
  {
    return (1);
  }
  else
  {
    return (0);
  }
}

################################################################################
# sub get_mapfile_jsvn: get jsoc version number from mapfile                   #
################################################################################
sub get_mapfile_jsvn($)
{
  my( $mapfile,  $dn_fn );
  $mapfile=$_[0];

  #found map file for this APID
  $dn_fn = sprintf("%s/%s",  $ENV{'HK_JSVN_MAP_DIRECTORY'}, $mapfile);

  #open map file
  open(FILE1, "$dn_fn") || die "(3)Can't Open:$dn_fn $!\n";

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
      ;#print "skipping\n";
      ;#print "fvn  is  $s_line[5]\n";
    }
  }#while
  close FILE1;
 return($jsvn);
}

################################################################################
# sub check_if_merge_pvn:check if packet version number is trigger for         #
# writing HK_SOURCE keyword to series                                          #
################################################################################
sub check_if_merged_pvn($$$)
{
  my( $pv_wn, $pv_dn, $ptype);
  $ptype=$_[0];
  $pv_wn=$_[1];
  $pv_dn=$_[2];
  if ($ptype == HK_PROCESSING_TYPE && $pv_wn == 1 &&  $pv_dn >= 194)
  {
    return (1)
  }
  else
  {
    return (0)
  }
}
################################################################################
# sub check_if_merge_fvn: Check if file version number is trigger for writing  #
# HK_SOURCE keyword to series.                                                 #
################################################################################
sub check_if_merged_fvn($,$)
{
  my( $fv, $ptype);
  $ptype=$_[0];
  $fv=$_[1];
  if ($ptype == SDO_PROCESSING_TYPE )
  {
    return (1)
  }
  else
  {
    return (0)
  }
}

###############################################################################
# sub get_jsvn_from_old_files: gets jsoc series version number from old files
###############################################################################
sub get_jsvn_from_old_files($)
{
  # setup my variables
  my($apid,$dn,$fn,$mfn_list,$mapfile,$found_flag);
  $apid=$_[0];

  # get directory name to find jsvn map files
  $dn = sprintf("%s",  $ENV{'HK_JSVN_MAP_DIRECTORY'});

  # open map file
  opendir(DIR_JSVN_MAP, $dn) || die ":(2)Can't open:$!\n"; #open subdirectory

  # read in all map files and put in list 
  @mfn_list = readdir(DIR_JSVN_MAP); #get a list of directory contents
  $found_flag = "0";
  foreach $mapfile  (@mfn_list)
  {
    if(substr( int  $mapfile,0,4) ne  hex substr($apid,2,4))
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
