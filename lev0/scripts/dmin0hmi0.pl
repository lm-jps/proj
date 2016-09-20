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
          '5178969022@vtext.com',
          '6507436500@tmomail.net',
          '6509968590@txt.att.net',
          '4088964682@messaging.sprintpcs.com',
          '7192333079@vmobl.com',
          'mark.cheung@gmail.com';

$to_list = join ",", $to_list, 'green@lmsal.com',
        'shelbe@lmsal.com', 'weiliu@lmsal.com', 'wei.liu2004@gmail.com',
        'zoe@lmsal.com', 'jeneen@sun.stanford.edu',
        'rock@sun.stanford.edu', 'thailand@sun.stanford.edu',
        'wolfsonjake@gmail.com', 'mbobra@sun.stanford.edu',
        'aanorton@stanford.edu', 'baldner@sun.stanford.edu',
        'william.d.pesnell@nasa.gov', 'cheung@lmsal.com',
        '4084317110@txt.att.net'
        unless $test;

$msg_file = "$ENV{HOME}/hmi_cam_anomaly.txt";
if (-e $msg_file) {
  $subj = "Stale HMI $msg_file during camera anomaly test";
  `echo "$subj" | mail -s "$subj" $to_list` if $test;
  exit;
}
$cname[1] = "1=vector=side";
$cname[2] = "2=Doppler=front";
$cmd = "/home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/show_info";
for ($cam=1; $cam<3; $cam++) {
  $fsn0 = `$cmd -q key=fsn 'hmi.lev0a[:#\$]'` - 399;
  $n = `$cmd -qc 'hmi.lev0a[$fsn0/400][?datamin=0?][?camera=$cam?]'` + 0;
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
