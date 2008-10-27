#!/usr/bin/perl
#Gets dayfile information from dayfiles series like size
#Script not to automated-need to edit at least items 1 and 2
$ENV{'PATH'}="/usr/local/bin:/bin:/usr/bin:.:/home/carl/cvs/myprod/JSOC/bin/linux_x86_64";

#(1)hardcode end day for month
#for october 2008 use 1-31 days
$i=1;
while ($i < 32)
{
  push(@time_list, $i);
  $i++;
}
$j=1;
while ($j < 64)
{
  push(@apid_list, $j);
  $j++;
}


#create time string for each october day
foreach $apid (@apid_list)
{
  foreach $time (@time_list)
  {
    #(2)hardcode month
    $time_str=sprintf("2008.10.%02s_00:00:00_TAI",$time);
    print "For DATE:$time_str   ";

    print "For APID:$apid\n";
    if($apid < 33)
    {
      ($si,$file)=`show_info \"ds=hmi.hk_dayfile[$time_str][$apid][moc]\" -p seg=file`;
    }
    elsif ($apid < 64 and $apid > 32)
    {
      ($si,$file)=`show_info \"ds=aia.hk_dayfile[$time_str][$apid][moc]\" -p seg=file`;
    }
    else
    {
      $file="";
    }
    if ($file eq "")
    {
      print "no file for this date!\n\n";
    }
    else
    {
      #should update to get file size-if fsize is zero don't print.
      $log=`ls -lt $file`;
      print "do ls -lt on dayfile:$log\n";
    }
  }
}


