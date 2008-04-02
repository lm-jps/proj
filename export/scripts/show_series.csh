#! /bin/csh -f
rm -f /tmp/ptt /tmp/pttresp
./show_series  > /tmp/pttresp

echo "Content-type: text/html" >>/tmp/ptt
echo "" >>/tmp/ptt
echo "<html><body><pre>" >>/tmp/ptt
echo "$status before ########################" >>/tmp/ptt
cat /tmp/pttresp >>/tmp/ptt
echo "$status after ########################" >>/tmp/ptt
date >>/tmp/ptt
env >>/tmp/ptt
echo "</pre></body></html>" >>/tmp/ptt

cat /tmp/pttresp
