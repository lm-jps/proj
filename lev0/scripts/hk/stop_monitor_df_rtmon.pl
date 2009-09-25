#!/usr/bin/perl
# NAME: STOP MONITOR_DF_RTMON.PL Script
# DESCRIPTION:stops monitor_df_rtmon.pl script from running. 
# AUTHOR:carl
# RUN EXAMPLE: stop_monitor_df_rtmon.pl apid=129
# RUN EXAMPLE: stop_monitor_df_rtmon.pl apid=29
# LIMITATION:Run on machine with access to process running 
# monitor_df_rtmon.pl script.
# DATE: 9/25/2009

#get apid argument
if (substr($ARGV[0],0,5) eq "apid=" )
{
  #use apid following -a for apid
  $apid= substr($ARGV[0],5);
  if ($apid eq "")
  {
    print "USAGE: stop_monitor_df_rtmon.pl  apid=129\n";
    exit;
  }
}
else
{
  print "USAGE: stop_monitor_df_rtmon.pl apid=129\n";
  exit;
}

#get process line
$processline=`ps -eaf |grep "\/monitor_df_rtmon.pl" | grep apid=$apid |grep -v grep`;

#get process id
@s_line= split / \s*/,  $processline;
if ($s_line[1])
{
  push(@list,$s_line[1]);
}

#if have process id on list then kill process.
if(@list)
{
 print "--->Killing process @list for command line: $s_line[7] $s_line[8] $s_line[9]\n";  
 kill 9, @list;
 print "--->Killed process @list for command line: $s_line[7] $s_line[8] $s_line[9]\n";  
}
else
{
  print "-->Cannot kill this process, since no such process exists for <monitor_df_rtmon.pl apid $apid>\n";
}
