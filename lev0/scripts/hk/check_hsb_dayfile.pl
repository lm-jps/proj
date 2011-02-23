#!/usr/local/bin/perl
############################################################################ 
#check_hsb_dayfile.pl is use to show whether to save dayfiles with
#hsb or hsb_r. This information is then used to run ingest_dayfile.pl
#to save hsb dayfiles to series.
#created 2/23/2011
#created by:carl
############################################################################ 
$infile="./hsb.txt";
$log=`sort ./hsb.txt > ./sort.hsb.txt`;
$sort_infile="./sort.hsb.txt";
#get apid argument
if (substr($ARGV[0],0,2) eq "-h")
{
    print "USAGE: % check_hsb_dayfile.pl\n";
    print "Requirements:\n";
    print "(1) Create a hsb.txt file from email containing all files to process and save to ./hsb.txt file\n";
    print "(2) When editing email of hsb files remove all text except filenames.\n";
    print "(3) If file not at ./ then won't work unless adjust infile variable at top of script.\n";
    print "(4) The output can be redirected to file (i,e, check_hsb_dayfile.pl > OUTPUT).\n";
    print "(5) The output tells you which files to save as \"hsb\" or \"hsb_r\" when using ingest_dayfile.pl script.\n";
    exit;
}

open(FILE, "$sort_infile") || die "Can't Open: <$sort_infile> file : $!\n";
while (<FILE>)
{
  if(m/hsb_/) 
  {
    #print "file got was $_";
    #check if in series
    # strip time value
    @s_time = split(/_/, $_);
    #print "$s_time[0]\n";
    #print " $s_time[1]\n";
    $s_time[1] = int $s_time[1];
    #print "$s_time[1]\n";
    #print "$s_time[2]\n";
    #print "$s_time[3]\n";
    #print "$s_time[4]\n";
    #print "$s_time[5]\n";
    $time_str=sprintf("%s.%s.%s_00:00:00_UTC",$s_time[2],$s_time[3],$s_time[4]);
    #print "$time_str\n"; 
    if(m/445|475|451|448|481|478/) 
    {
      $out=`show_info ds=hmi.hk_dayfile[$time_str][$s_time[1]][hsb] seg=file`;
      #print "<$out>\n";
    }
    elsif (m/529|536|540|576|580|569/) 
    {
      $out=`show_info ds=aia.hk_dayfile[$time_str][$s_time[1]][hsb] seg=file`;
      #print "<$out>\n";
    }
    #check Result
   if ($out=~/No record/)
   {
      print "THIS DAYFILE DOES NOT HAVE PRIMARY HSB SOURCE RECORD- USE \"hsb\" to this load file===>$_";
   }
   else
   {
      print "THIS DAYFILE HAS A PRIMARY HSB SOURCE RECORD THAT EXISTS!- USE \"hsb_r\" to load this file---------->$_";

   }
  }
}
 
