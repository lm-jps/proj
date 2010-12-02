#!/usr/bin/perl
#NAME:patch01_hmi_temperature.pl
#patch01 for hmi.temperature_summary_300s series to add records in for times shown below.
#approved by Rock on 12/2/2010 to run on series.
#the script creates records and sets only T_START value and NUMPTS values to zero.
my $series= "hmi.temperature_summary_300s";
my $tstart="";
my @tstart_array;
push (@tstart_array, "2010.11.29_19:55:00_UTC");
push (@tstart_array, "2010.11.29_20:00:00_UTC");
push (@tstart_array, "2010.11.29_20:05:00_UTC");
push (@tstart_array, "2010.11.29_20:10:00_UTC");
push (@tstart_array, "2010.11.29_20:15:00_UTC");
push (@tstart_array, "2010.11.29_20:20:00_UTC");
push (@tstart_array, "2010.11.29_20:25:00_UTC");
push (@tstart_array, "2010.11.29_20:30:00_UTC");
push (@tstart_array, "2010.11.29_20:35:00_UTC");
push (@tstart_array, "2010.11.29_20:40:00_UTC");

foreach $tstart (@tstart_array)
{
print "-->Writing to series <$series> for T_START <$tstart>\n";
$log=`/home/carl/cvs/JSOC/bin/linux_x86_64/set_keys -c  ds=$series  T_START=$tstart NUMPTS=0`;
print "-->Running set_keys -c ds=$series T_START=$tstart NUMPTS=0 \n-->Log Errors:<$log>\n";
}


