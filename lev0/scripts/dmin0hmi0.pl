#!/usr/bin/perl -w
$to_list = join ",", 'jps@lmsal.com', 'boerner@lmsal.com', 'green@lmsal.com',
        'wolfson@lmsal.com', 'zoe@lmsal.com', 'jeneen@sun.stanford.edu',
        'rock@sun.stanford.edu', 'thailand@sun.stanford.edu',
        'couvidat@stanford.edu',
        '6509965043@txt.att.net',
        '6504503716@txt.att.net',
        '6502248075@txt.att.net',
        '4084317110@txt.att.net',
        '5103252489@vtext.com',
        '5178969022@vtext.com',
        '6503878335@txt.att.net',
        '6507436500@tmomail.net',
        '6509968590@txt.att.net';
$msg_file = "$ENV{HOME}/hmi_cam_anomaly.txt";
exit if -e $msg_file;
$cname[1] = "1=vector=side";
$cname[2] = "2=Doppler=front";
if ($t = shift @ARGV) { $threshold = $t; } else { $threshold = 100; }
$cmd = "/home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/show_info";
for ($cam=1; $cam<3; $cam++) {
  $fsn0 = `$cmd -q key=fsn 'hmi.lev0a[:#\$]'` - 399;
  $n = `$cmd -qc 'hmi.lev0a[$fsn0/400][?datamin=0?][?camera=$cam?]'` + 0;
#printf "fsn0: %d, n: %d, Thresh: %d\n", $fsn0, $n, $threshold;
  if ($n > $threshold) {
    $c = $cname[$cam];
    $subj = "HMI camera $c anomaly";
    open MSG, ">$msg_file" or die "Can not write mail";
    print MSG "\n", `date -u`,`date`, "$subj detected.\n";
    close MSG or die "Can not close mail message file";
    `mail -s "$subj" $to_list < $msg_file`;
  }
}
