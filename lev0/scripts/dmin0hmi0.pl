#!/usr/bin/perl -w
# Count number of HMI images per camera with DATAMIN = 0 for 400 most
# recent recnums (about 200 per camera). If count > threshold (default
# 100) send email to list. Body of email is written to
# $ENV{HOME}/hmi_cam_anomaly.txt which also serves as a flag that camera
# anomaly was detected. To prevent flood of mail at cron cadence
# $ENV{HOME}/hmi_cam_anomaly.txt must be deleted after anomaly is fixed
# to re-enable automatic anomaly detection.
if ($t = shift @ARGV) { $threshold = $t; } else { $threshold = 100; }
$test = shift @ARGV or $test = 0;
$to_list = join ",", 'jps@lmsal.com',
          '6509965043@txt.att.net',
          '6504503716@txt.att.net',
          '5178969022@txt.att.net',
          '6507436500@tmomail.net',
          '4088964682@messaging.sprintpcs.com',
          'mark.cheung@gmail.com';

$to_list = join ",", $to_list, 'green@lmsal.com',
        '4087208116@txt.att.net', 'weiliu@lmsal.com', 'wei.liu2004@gmail.com',
        'jeneen@sun.stanford.edu',
        'rock@sun.stanford.edu', 'thailand@sun.stanford.edu',
        'wolfsonjake@gmail.com', 'mbobra@sun.stanford.edu',
        'aanorton@stanford.edu', 'baldner@sun.stanford.edu',
        'william.d.pesnell@nasa.gov', 'cheung@lmsal.com',
        '4086125678@vtext.com', 
        '4089338703@txt.att.net'
        unless $test;

$msg_file = "$ENV{HOME}/hmi_cam_anomaly.txt";
if (-e $msg_file) {
  $subj = "Stale HMI $msg_file during camera anomaly test";
  `echo "$subj" | mail -s "$subj" $to_list` if $test;
  exit;
}
$t = time - 86400;
my ($sc, $mn, $hr, $da, $mo, $yr) = gmtime($t);
$ts = sprintf '$(%d.%2.2d.%2.2d_00:00)', $yr+1900, $mo+1, $da;
$cname[1] = "1=vector=side";
$cname[2] = "2=Doppler=front";
$cmd = "/home/jsoc/cvs/Development/JSOC/bin/linux_avx/show_info";
for ($cam=1; $cam<3; $cam++) {
  $fsn0 = `$cmd -q key=fsn 'hmi.lev0a[:#\$]'` - 399;
  $qs = "hmi.lev0a[$fsn0/400][?datamin=0?][?camera=$cam?][?T_OBS>$ts?]";
  $n = `$cmd -qc '$qs'` + 0;
  if ($n > $threshold) {
    $c = $cname[$cam];
    $subj = "HMI camera $c anomaly";
    $subj = join " ", $test, $subj if $test;
    open MSG, ">$msg_file" or die "Can not write mail";
    print MSG "\n", `date -u`,`date`, "$subj detected.\n";
    close MSG or die "Can not close mail message file";
    `mail -s "$subj" $to_list < $msg_file`;
  }
}
unlink $msg_file if $test;
