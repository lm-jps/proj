#!/usr/bin/perl -w
# Count number of AIA images per camera with AIHIS192 > 9 for 4000 most
# recent recnums (about 1000 per camera). If count > threshold (default 5)
# send email to list. Body of email is written to $ENV{HOME}/bit_flip_his.txt
# which also serves as a flag that camera anomaly was detected. To prevent
# flood of mail at cron cadence, $ENV{HOME}/bit_flip_his.txt must be
# deleted after anomaly is fixed to re-enable automatic anomaly detection.
if ($t = shift @ARGV) { $threshold = $t; } else { $threshold = 5; }
$test = shift @ARGV or $test = 0;
$to_list = join ",", 'jps@lmsal.com',
        '6509965043@txt.att.net',
        '6504503716@txt.att.net',
        '4084317110@txt.att.net',
        '5103252489@vtext.com',
        '5178969022@vtext.com',
        '6507436500@tmomail.net',
        '6509968590@txt.att.net';
$to_list = join ",", $to_list, 'boerner@lmsal.com', 'green@lmsal.com',
        'zoe@lmsal.com', 'jeneen@sun.stanford.edu',
        'rock@sun.stanford.edu', 'thailand@sun.stanford.edu',
        'wolfsonjake@gmail.com', 'mbobra@sun.stanford.edu',
        'aanorton@stanford.edu', 'baldner@sun.stanford.edu',
        'william.d.pesnell@nasa.gov', 'cheung@lmsal.com'
        unless $test;
$msg_file = "$ENV{HOME}/bit_flip_his.txt";
if (-e $msg_file) {
  $subj = "Stale AIA $msg_file during camera anomaly test";
  `echo "$subj" | mail -s "$subj" $to_list` if $test;
  exit;
}
$cmd = "/home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/show_info";
for ($cam=1; $cam<5; $cam++) {
  $fsn0 = `$cmd -q key=fsn 'aia.lev0[:#\$]'` - 3999;
  $n = `$cmd -c 'aia.lev0[$fsn0/4000][?aihis192>9?][?camera=$cam?]'`;
  @w = split /\s+/, $n;
  if ($w[0] > $threshold) {
    $subj = "AIA camera $cam histogram anomaly";
    $subj = join " ", $test, $subj if $test;
    open MSG, ">$msg_file" or die "Can not write mail";
    print MSG "\n", `date -u`,`date`, "$subj detected.\n";
    close MSG or die "Can not close mail message file";
    `mail -s "$subj" $to_list < $msg_file`;
  }
}
unlink $msg_file if $test;
