#!/usr/local/bin/perl
#hsb_0445_2011_02_20_00_00_02_00.hkt 
$infile="./hsb.txt";
$log=`sort ./hsb.txt > ./sort.hsb.txt`;
$sort_infile="./sort.hsb.txt";
open(FILE, "$sort_infile") || die "Can't Open: <$sort_infile> file: $!\n";
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
 
