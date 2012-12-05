#!/usr/bin/perl
#
#This is called as a cron job, typically after midnight,
#to restart the socdc script that runs the datacapture front 
#end on dcs0 and dcs1.
#
`touch /usr/local/logs/soc/RESTARTNOW`;

