#! /bin/csh -f

set TMP = $$.tmp
echo -n "escaped_html='" > $TMP
./escape_file.pl < export_request_form.html >> $TMP
echo "';" >> $TMP

mv $TMP export_request_form.htmlesc
