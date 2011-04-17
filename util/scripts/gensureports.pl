#!/home/jsoc/bin/linux_x86_64/perl5.12.2 -w

use constant kWdPath => "/web/jsoc2/htdocs/sureports";
use constant kTmpFile => "tmp.txt";
use constant kDataSec => "__DATA__";
use constant kEndSec => "__END__";
use constant kScrPath => "/home/jsoc/cvs/Development/JSOC/proj/util/scripts";

my($res);
my($cmd);
my($line);
my($datasec);


# All Data
$cmd = kScrPath . "/sumstapestat.pl jsoc_sums hmidb 5434 production agg group all > " . kWdPath . "/" . kTmpFile;
`$cmd`;

if (open(TMPFILE, "<" . kWdPath . "/" . kTmpFile) && open(REPFILE, ">" . kWdPath . "/sureportA.txt"))
{
   FilterData(\*TMPFILE, \*REPFILE);

   close(TMPFILE);
   close(REPFILE);
}

# AP Data
$cmd = kScrPath . "/sumstapestat.pl jsoc_sums hmidb 5434 production agg group ap > " . kWdPath . "/" . kTmpFile;
`$cmd`;

if (open(TMPFILE, "<" . kWdPath . "/" . kTmpFile) && open(REPFILE, ">" . kWdPath . "/sureportB.txt"))
{
   FilterData(\*TMPFILE, \*REPFILE);

   close(TMPFILE);
   close(REPFILE);
}

exit(0);

# subroutines
sub FilterData
{
   my($tmpfile) = $_[0];
   my($repfile) = $_[1];
   my($datas) = kDataSec;
   my($ends) = kEndSec;
 
   $datasec = 0;
   while (defined($line = <$tmpfile>))
   {
      chomp($line);

      if ($line =~ /$datas/)
      {
         $datasec = 1;
         next;
      }
      elsif ($line =~ /$ends/)
      {
         $datasec = 0;
         last;
      }

      if ($datasec)
      {
         print $repfile "$line\n";
      }
   }     
}
