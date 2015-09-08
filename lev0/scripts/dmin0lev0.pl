#!/usr/bin/perl -w
#$to_list = join ",", 'jps@lmsal.com', 'zoe@lmsal.com', 'wolfson@lmsal.com',
#        'emma.m.lehman@lmco.com', 'linford@lmsal.com',
#        '6505213455@txt.att.net', '6502242826@txt.att.net',
#        '6504503716@txt.att.net', '4084317110@txt.att.net',
#        '6502248075@txt.att.net', '6503878335@txt.att.net',
#        '4153141403@vtext.com',   '6034963597@vtext.com';
$to_list = join ",", 'jps@lmsal.com', 'wolfson@lmsal.com',
                      'zoe@lmsal.com';
$msg_file = "$ENV{HOME}/bit_flip.txt";
exit if -e $msg_file;
if ($t = shift @ARGV) { $threshold = $t; } else { $threshold = 200; }
$cmd = "/home/jsoc/cvs/Development/JSOC/bin/linux_x86_64/show_info";
for ($cam=1; $cam<5; $cam++) {
  $fsn0 = `$cmd -q key=fsn 'aia.lev0[:#\$]'` - 3999;
  $n = `$cmd -qc 'aia.lev0[$fsn0/4000][?datamin=0?][?camera=$cam?]'` + 0;
  if ($n > $threshold) {
    $subj = "AIA camera $cam anomaly";
    open MSG, ">$msg_file" or die "Can not write mail";
    print MSG "$subj at approximately\n",`date -u`,`date`;
    close MSG or die "Can not close mail message file";
    `mail -s "$subj" $to_list < $msg_file`;
  }
}
