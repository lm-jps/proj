#! /bin/csh -f
rm -f /tmp/ptt*
date > /tmp/ptt
./jsoc_info  > /tmp/pttinfo

echo "Content-type: text/html" >>/tmp/ptt
echo "" >>/tmp/ptt
echo "<html><body><pre>" >>/tmp/ptt
echo "$status before ########################" >>/tmp/ptt
cat  /tmp/pttinfo >>/tmp/ptt
echo "$status after ########################" >>/tmp/ptt
env >>/tmp/ptt
date >>/tmp/ptt
echo "</pre></body></html>" >>/tmp/ptt

cat /tmp/pttinfo
#cat /tmp/ptt
