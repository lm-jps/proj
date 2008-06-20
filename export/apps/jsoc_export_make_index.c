#define DEBUG 1
#define DEBUG 0

/*
 *  jsoc_export_make_index - Generates index.XXX files for dataset export.
 *  Should be run in the directory containing a jsoc export index.txt file.
 *  Will read the index.txt and generate index.json and index.html
 *
*/

#include <stdio.h>
#include "json.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <strings.h>

// #include "json.h"

static char x2c (char *what) {
  register char digit;

  digit = (what[0] >= 'A' ? ((what[0] & 0xdf) - 'A')+10 : (what[0] - '0'));
  digit *= 16;
  digit += (what[1] >= 'A' ? ((what[1] & 0xdf) - 'A')+10 : (what[1] - '0'));
  return (digit);
}

static void CGI_unescape_url (char *url) {
  register int x, y;

  for (x = 0, y = 0; url[y]; ++x, ++y) {
    if ((url[x] = url[y]) == '%') {
      url[x] = x2c (&url[y+1]);
      y += 2;
    }
  }
  url[x] = '\0';
}

char * string_to_json(char *in)
  {
  char *new, *c = in;
  wchar_t *tmp, *work;
  tmp = work = (wchar_t *)malloc(sizeof(wchar_t)*(strlen(in)+1));
  while (c && *c) *tmp++ = *c++;
  *tmp = NULL;
  new = json_escape(work);
  free(work);
  return(new);
  }

#define HTML_INTRO "<HTML><HEAD TITLE=\"JSOC Data Export Request Results Index\">\n" \
		"</HEAD><BODY>\n" \
		"<H2>JSOC Data Request Summary</H2>\n" \
		"<TABLE>\n"
		

/* Module main function. */
int main(int argc, char **argv)
  {
  FILE *index_txt, *index_html, *index_json;
  char buf[1000];
  char dir[1000];
  json_t *jroot;
  json_t *recinfo;
  json_t *fileinfo;
  char *json, *final_json;
  int state = 0;

  if (argc && strcmp(*argv, "-h") == 0)
    {
    printf ("Usage:\njsoc_export_make_index {-h}\n");
    return(0);
    }

  index_txt = fopen("index.txt","r");
  if (!index_txt)
    {
    fprintf(stderr,"XX jsoc_export_make_index - failed to open index.txt.\n");
    return(1);
    }

  index_html = fopen("index.html","w");
  if (!index_html)
    {
    fprintf(stderr,"XX jsoc_export_make_index - failed to create index.html.\n");
    return(1);
    }
  fprintf(index_html, HTML_INTRO);

  index_json = fopen("index.json","w");
  if (!index_json)
    {
    fprintf(stderr,"XX jsoc_export_make_index - failed to create index.json.\n");
    return(1);
    }
  jroot = json_new_object();
  recinfo = json_new_array();

  while (fgets(buf, 1000, index_txt))
    {
    char *p = buf + strlen(buf) - 1;
    if (p >= buf && *p == '\n')
      *p = '\0'; 
    p = buf;
    switch (state)
      {
      char *name, *val;
      char *namestr, *valstr;
      case 0: // Initial read expect standard header line
	if (strncmp(buf, "# JSOC ",7) != 0)
	  {
	  fprintf(stderr, "XX jsoc_export_make_index - incorrect index.txt file.\n");
	  return(1);
	  }
	state = 1;
	break;
      case 1:  // In header section, take name=val pairs.
	if (strncmp(buf, "# DATA",6) == 0) // done with header ?
	  {
	  state = 2;
	  fprintf(index_html, "</TABLE><P><H2><B>Selected Data</B></H2><P><TABLE>\n");
	  break;
	  }
	if (*p == '#' || !*p) // skip blank and comment lines
	  break;
	if ((val=index(p, '='))==NULL)
	  {
	  fprintf(stderr, "XX jsoc_export_make_index - ignore unexpected line in header, \n    %s\n", buf);
	  break;
	  }
	p = val++;
	name = buf;
	while (isblank(*name))
	  name++;
        *p-- = '\0';
	while (isblank(*p) && p >= buf)
          *p-- = '\0';
        while (isblank(*val))
          val++;
	p = val + strlen(val);
        p--;
	while (isblank(*p) && p >= val)
          *p-- = '\0';

	// save dir for use in data section
	if (strcmp(name, "dir") == 0)
	  strncpy(dir, val, 1000);
	// put name=value pair into index.json
	namestr = string_to_json(name);
        valstr = string_to_json(val);
	json_insert_pair_into_object(jroot, namestr, json_new_string(valstr));
	free(namestr);
	free(valstr);
	// put name=value pair into index.html
	fprintf(index_html, "<TR><TD><B>%s</B></TD<TD>%s</TD></TR>\n", name, val);
	break;
      case 2: // Data section contains pairs of record query and filenames
	if (*p == '#' || !*p) // skip blank and comment lines
	  break;
	name = buf;
	while (isblank(*name)) // skip leading blanks
	  name++;
        p = name;
	while (*p && !isblank(*p)) // skip past query
	  p++;
	if (*p)
	  *p++ = '\0'; // mark end of query
	val = p;
	while (isblank(*val)) // skip leading blanks
	  val++;
	p = val + strlen(val);
        p--;
	while (isblank(*p) && p >= val) // omit trailing blanks
          *p-- = '\0';
	// put query : filename pair into index.json
	fileinfo = json_new_object();
	namestr = string_to_json(name);
	json_insert_pair_into_object(fileinfo, "record", json_new_string(namestr));
	free(namestr);
        valstr = string_to_json(val);
	json_insert_pair_into_object(fileinfo, "filename", json_new_string(valstr));
	free(valstr);
	json_insert_child(recinfo, fileinfo);
	// put name=value pair into index.html
	fprintf(index_html, "<TR><TD>%s</TD<TD><A HREF=\"http://jsoc.stanford.edu/%s/%s\">%s</A></TD></TR>\n", name, dir, val, val);
	break;
      default:
      }
    }
  fclose(index_txt);

  // finish json
  json_insert_pair_into_object(jroot, "data", recinfo);
  json_tree_to_string(jroot,&json);
  final_json = json_format_string(json);
  fprintf(index_json, "%s\n",final_json);
  fclose(index_json);

  // finish html
  fprintf(index_html, "</TABLE></BODY></HTML>\n");
  fclose(index_html);

  return(0);
  }
