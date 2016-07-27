#!/usr/bin/perl -w
# Count number of HMI images per camera with DATAMIN = 0 for 400 most
# recent recnums (about 200 per camera). If count > threshold (default
# 100) send email to list. Body of email is written to
# $ENV{HOME}/hmi_cam_anomaly.txt which also serves as a flag that camera
# anomaly was detected. To prevent flood of mail at cron cadence
# $ENV{HOME}/hmi_cam_anomaly.txt must be deleted after anomaly is fixed
# to re-enable automatic anomaly detection.
$to_list = join ",", 'jps@lmsal.com', 'boerner@lmsal.com', 'green@lmsal.com',
        'zoe@lmsal.com', 'jeneen@sun.stanford.edu',
        'rock@sun.stanford.edu', 'thailand@sun.stanford.edu',
        'wolfsonjake@gmail.com', 'mbobra@sun.stanford.edu',
        'aanorton@stanford.edu', 'baldner@sun.stanford.edu',
        'william.d.pesnell@nasa.gov', 'cheung@lmsal.com',
        '6509965043@txt.att.net',
        '6504503716@txt.att.net',
        '4084317110@txt.att.net',
        '5103252489@vtext.com',
        '5178969022@vtext.com',
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
  if ($n > $threshold) {
    $c = $cname[$cam];
    $subj = "HMI camera $c anomaly";
    open MSG, ">$msg_file" or die "Can not write mail";
    print MSG "\n", `date -u`,`date`, "$subj detected.\n";
    close MSG or die "Can not close mail message file";
    `mail -s "$subj" $to_list < $msg_file`;
  }
}
