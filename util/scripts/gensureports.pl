#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

use constant kTmpFile => "/web/jsoc2/htdocs/sureports/tmp.txt";
use constant kDataSec => "__DATA__";
use constant kEndSec => "__END__";

my($res);
my($cmd);
my($line);
my($datasec);


# All Data
$cmd = "/home/jsoc/cvs/Development/JSOC/proj/util/scripts/sumstapestat.pl jsoc_sums hmidb 5434 production agg group all > " . kTmpFile;

if (open(REPFILE, "<" . kTmpFile))
{
   $datasec = 0;
   while (defined($line = <REPFILE>))
   {
      chomp($line);

      if ($line =~ /kDataSec/)
      {
         $datasec = 1;
      }
      elsif ($line =~ /kEndSec/)
      {
         $datasec = 0;
         last;
      }

      if ($datasec)
      {
         print "$line\n";
      }
   }
}


# AP Data
# `/home/jsoc/cvs/Development/JSOC/proj/util/scripts/sumstapestat.pl jsoc_sums hmidb 5434 production agg group ap > /web/jsoc2/htdocs/sureports/tmp.txt`;
