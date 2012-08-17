#!/usr/bin/perl
use CGI; 
while ($x=<>) { print CGI::escape($x); }
