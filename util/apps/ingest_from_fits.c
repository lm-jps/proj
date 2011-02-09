/* ingest_from_fits - read FITS file and print JSD, Keyword map file, or ingest file */
/*
 * ingest_from_fits [-j] [in=]<fitsfile> [map=<mapfile>] [ds=<series>]
 *   -j means print JSD for the fitsfile.
 *   -c means create the series if ds= is given and 
 *      a primekey is given and
 *      if the series does not exist already.
 *  in= is optional, but the fitsfile name is required.
 *  ds= gives series name.  If present a new record will be added
 *      containing the fitsfile keywords and data.
 *  map= Thea optional mapfile contains a line for each
 *       keyword that needs special action.
 *
 *  This program is intended to be called at least twice, first
 *  to get a sample JSD file then, once that JSD file has been
 *  modified and create_series has been called with the JSD file,
 *  this program can be called again to ingest the fitsfile.
 *  

 *  Several things need to be changed in the JSD template file:
 *    The target seriesname must be specified.
 *    The PrimeKeys and DBindex lines must be fixed or removed.
 *    The mapfile lines at the end should be extracted if needed
 *       but certainly removed from the JSD file.
 *
 *  If any illegal FITS keywords were encountered they will be
 *  added to the sample mapfile lines.  The mapfile structure is
 *    <drmskeyname> <fitskeyname> <action>
 *
 *  The code only recognizes the action "copy" but the place
 *  is marked where other actions can be added as needed.
 *
 *  NOTE:  You can change any desired attributes of keywords
 *  simply by changing the JSD entry.  If you change the name
 *  you will need to add a line in a mapfile to let the program
 *  match the correct fits keyword values with the new name.
 *
 *  NOTE:  The present code that reads the FITS file takes its own
 *  action on illegal DRMS keywords.  Keywords with trailing "_mh" have
 *  had hyphens changed to underscore.  If desired for later export
 *  from DRMS, these lines in the JSD should be changed to use the
 *  export name mapping rules.  Move the _mh name to the comment
 *  field in [] and change the key name to have double underscore
 *  instead of the "_mh" form.  Then make a map file to be used
 *  when ingesting fitsfiles into this series.
 *
 *  Example, this line:
 *    Keyword: DATE_OBS_mh, string, variable, record, "", "%s", "none", ""
 *  should be changed to:
 *    Keyword: DATE__OBS, string, variable, record, "", "%s", "none", "[DATE-OBS]"
 *
 *  And a map file should contain the line:
 *    DATE__OBS DATE_OBS_mh copy
 *
 *  Until the fitsrw_read code allows the original keyword to be put into the
 *  HContainer it returns with keywords, there is no way to do this in the
 *  code automatically.
 *
 *  This code makes only a single segment.  If you want to
 *  ingest multiple segments into a single record you can use
 *  fits_into_drms or use this program to make a JSD for each
 *  kind of fitsfile, then edit into a single JSD with multiple
 *  segments, create the series, then use fits_into_drms.
 *
 *  If the -c flag is given and a series and primekey are specified and the series does not
 *  already exist then JSD that is created will
 *  be used to create a series.  It will also be printed if the -j flag is present.
 */  

/**
\defgroup ingest_from_fits ingest_from_fits
@ingroup su_migration

\par Synopsis:
\code
ingest_from_fits [-j] [-c] {in=}<fitsfile> [ds=<seriesname>] [primekey=<primekeys>] [map=<mapfile>]
\endcode

\details

\b Ingest_from_fits provides tools to aid ingesting FITS files into DRMS.
It can help by making a draft JSD, and optionally by creating a new series from that JSD, and/or
ingesting a fitsfile into the specified series.

\par Options:

The program can make a JSD and exit, make and use a JSD to create a series, and/or
ingest a fits file into a (possible new) series.

\li \c -c: Create a new series using the <seriesname> given in the ds= argument.
\li \c -j: Create a JSD from a given <fitsfile> and print it to stdout.
\li \c {in=}<fitsfile> - Specifies the FITS file to be ingested and/or used to generate a JSD.
\li \c [ds=<seriesname>] - specifies the target DRMS series to be used in the generated JSD and/or to create and/or to insert data into.
\li \c [primekey=<primekeys>] - specifies one or more keywords to use as PrimeKey and DBindex in the JSD.
\li \c [map=<mapfile>] - specifies an optional keyword mapfile to use while ingesting the fitsfile.

\par Usage:

In the base mode with only the -j flag and a fitsfile provided, ingest_from_fits will
print a draft JSD file that can be captured and editted by the user.  The properly
editted JSD file can then be used with \ref create_series to generate a new
DRMS dataseries that is appropriate to use when ingesting fits files like the
sample used.  To then actually ingest the fits file into the new series,
call ingest_from_fits with the fitsfile and the seriesname passed in the ds= argument.

The draft JSD file generated does not contain a seriesname nor PrimeKey or DBindex lists.
However, if the primekey= argument is provided the given <primekeys> string will be
put into both the PrimeKeys and DBindex fields in the JSD.  If the ds= argument is
given, then the priovided <seriesname> will be put into the Seriesname field in the JSD.

Thus if both the ds= and primekey= argument are given, a complete JSD is created.

If the -c flag is given as well as the ds= and primekey= arguments, then a new
series will be created with the given seriesname.  If it already exists an error
message is printed and the program quits.

\b NOTE the default JSD has no archiving and a retention time of 10 days.

If the -j flag is NOT given but the ds= argument is given then the fitsfile will be ingested
into that series.  If the -c flag is given then the new record will be the first record
in the new series.

If in the process of generating a JSD from the fitsfile, some illegal DRMS names are
found among the FITS keywords, then two lines will be printed for each such keyword.
The first line will be a comment with the original FITS name and the auto-generated
substitute name.  Next a sample mapfile line will be provided which can be editted
and included in a mapfile if desired.  This second line contains first the desired
DRMS name, then the name to be found in the input file (this needs to be the auto-converted
name), then an action.  The defualt action is "copy".
If you do not want the auto-generated substitute keyword name, change the first column
to the desired name AND change the matching line in the draft JSD to also have
the desired name.
This
is the keyword mapping format also used by \ref ingest_dsds_a and can be captured from
the ingest_from_fits stdout into a mapfile.  That mapfile, after possible editting,
can be given to subsequent calls of ingest_from_fits to be used when ingesting fits files.
The original FITS keyword will be placed into the JSD in the "note" section of the Keyword
line.  Then upon export via e.g. \ref jsoc_export_as_fits the keyword will be mapped
back into the original name.

When a mappped keyword is encountered in the ingest process, an action is taken depending
on the value of the "action" field in the mapfile.  In the present code, only the "copy"
action is implemented but the place for the user to add special code for other user
defined actions is marked in the code.  See \ref ingest_dsds_a for examples.
The keyword list in an ingested fitsfile is inspected for illegal names even in the
case where the -j flag is not given and data is simply ingested into a series.
In this case the <mapfile> and any newly found bad keywords are merged with
the mapfile taking precedence.

\par Output:

Stderr:  Some output may be generated in the internal call of \ref fitsrw_read which scans the
<fitsfile> to make a list of keywords.  Multiple instances of a given keyword for instance will
generate information lines.  Other diagnostics are also directed to stderr.

Stdout: The normal output stream is reserved for the generated JSD information and self-generated
<mapfile> entries if some of the keywords need to be mapped to DRMS compliant keyrord names.
If the stdout is captured into a file to be used as a JSD, then the <mapfile> lines at the end
should be extracted to a separate mapfile for later use.

\par Examples:

\b Example 1:
To print a draft JSD file appropriate for ingesting e.g. a MDI magnetogram with the "coffee-cup sunspot":
\code
  ingest_from_fits -j /mag//fd_M_96m_01d.001994/fd_M_96m_01d.1994.0010.fits
\endcode

\b Example 2:
To make the same JSD but with specified Seriesname and Primekeys then make a series manually:
\code
  ingest_from_fits -j /mag/fd_M_96m_01d.001994/fd_M_96m_01d.1994.0010.fits ds=su_phil.test primekey=T_REC >pt.jsd
  create_series pt.jsd
\endcode

\b Example 3:
To ingest several fits files into the series created in example 2:
\code
  cd /mag/fd_M_96m_01d.001994
  foreach fitsfile ( *[0-9].fits )
    ingest_from_fits ds=su_phil.test $fitsfile
  end
\endcode

\b Example 4:
To ingest a single fitsfile into a not-yet created series, all in one command:
\code
  ingest_from_fits -c /mag/fd_M_96m_01d.001994/fd_M_96m_01d.1994.0010.fits ds=su_phil.test primekey=T_REC
\endcode

\bug

*/


#include "jsoc_main.h"

char *module_name = "ingest_from_fits";

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT_SPECIFIED",  "Input FITS file."},
     {ARG_STRING, "ds", "NOT_SPECIFIED",  "Target DRMS data series."},
     {ARG_STRING, "map", "NOT_SPECIFIED",  "Map file for newname from oldname."},
     {ARG_STRING, "primekey", "NOT_SPECIFIED",  "keyword name to use as a prime key in the JSD created"},
     {ARG_FLAG, "c", "0",  "Use generated jsd to create a series. Requires both ds and primekey args."},
     {ARG_FLAG, "j", "0",  "Print jsd."},
     {ARG_END}
};

# define MAXJSDLEN 100000
# define MAXMAPLEN 10000

int DoIt(void)
  {
  const char *in = params_get_str(&cmdparams, "in");
  const char *ds = params_get_str(&cmdparams, "ds");
  const char *map = params_get_str(&cmdparams, "map");
  const char *primekey = params_get_str(&cmdparams, "primekey");
  int printjsd = params_isflagset(&cmdparams, "j");
  int wantcreate = params_isflagset(&cmdparams, "c");
  int haveseriesname =  strcmp(ds, "NOT_SPECIFIED") != 0;
  int haveprime =  strcmp(primekey, "NOT_SPECIFIED") != 0;
  int havemap = strcmp(map, "NOT_SPECIFIED") != 0;
  int wantjsd = printjsd || wantcreate;
  int insertrec = haveseriesname  && !printjsd;
  int status = DRMS_SUCCESS;
  DRMS_Array_t *data = NULL;
  DRMS_Keyword_t *key=NULL;
  HContainer_t *keywords = NULL;
  HIterator_t hit;
  int readraw = 1;
  char jsd[MAXJSDLEN];
  char *newnames[MAXMAPLEN];
  char *oldnames[MAXMAPLEN];
  char *actions[MAXMAPLEN];
  int imap, nmap = 0;

  if (strcmp(in, "NOT_SPECIFIED") == 0)
    {
    if (cmdparams_numargs(&cmdparams) < 1 || !(in = cmdparams_getarg(&cmdparams, 1)))
      DIE("No input data found");
    }

  data = drms_fitsrw_read(drms_env, in, readraw, &keywords, &status);
  if (status || !keywords)
  {
     DIE("No keywords found");
  }

  if (wantjsd)
    {
    char *pjsd = jsd;
    char keyname[DRMS_MAXNAMELEN];
    char *colon;
    DRMS_Type_t datatype;
    int iaxis, naxis, dims[10];
    double bzero, bscale;

    // build jsd in internal string
    pjsd += sprintf(pjsd, "#=====General Series Information=====\n");
    pjsd += sprintf(pjsd, "Seriesname:  %s\n", (haveseriesname ? ds : "<NAME HERE>"));
    pjsd += sprintf(pjsd, "Author:      %s\n", getenv("USER"));
    pjsd += sprintf(pjsd, "Owner:       nobody_yet\n");
    pjsd += sprintf(pjsd, "Unitsize:    1\n");
    pjsd += sprintf(pjsd, "Archive:     0\n");
    pjsd += sprintf(pjsd, "Retention:   10\n");
    pjsd += sprintf(pjsd, "Tapegroup:   0\n");
    pjsd += sprintf(pjsd, "PrimeKeys:   %s\n", (haveprime ? primekey : "<PRIME KEYS HERE OR DELETE LINE>"));
    pjsd += sprintf(pjsd, "DBIndex:     %s\n", (haveprime ? primekey : "<PRIME KEYS HERE OR DELETE LINE>"));
    pjsd += sprintf(pjsd, "Description: \"From: %s\"\n", in);
  
    pjsd += sprintf(pjsd, "#===== Keywords\n");
  
    // drms_fitsrw_read() does not place reserved fits keywords in the keywords container.
    // The BITPIX, NAXIS, BLANK, BZERO, BSCALE, SIMPLE, EXTEND values are copied or
    // set in various fields in the in the DRMS_Array_t struct returned. END is dropped.
    // Another function, fitsrw_read(), WILL put every FITS keyword into the keywords
    // container, but it does not convert their names into DRMS-compatible keyword names, 
    // unlike drms_fitsrw_read().
    datatype = data->type;
    naxis = data->naxis;
    memcpy(dims, data->axis, sizeof(int) * naxis);
    bzero = data->bzero;
    bscale = data->bscale;

    hiter_new (&hit, keywords);
    while ( key = (DRMS_Keyword_t *)hiter_getnext(&hit) )
      {
      strcpy(keyname, key->info->name);
  
      colon = index(key->info->description, ':');
      // check for lllegal or reserved DRMS names 
      // In this case the FITS Keyword structure note section will contain the
      // original FITS keyword.
      if (*(key->info->description) == '[')
        {
        char *c;
	char originalname[80];
        strcpy(originalname, key->info->description+1);
        c = index(originalname, ':');
        if (c)
          *c = '\0';
        c = index(originalname, ']');
        if (c)
          *c = '\0';
        newnames[nmap] = strdup(keyname);
        oldnames[nmap] = strdup(originalname);
        actions[nmap] = strdup("copy");
        nmap++; 
        }
      
      pjsd += sprintf(pjsd, "Keyword: %s, ", keyname);
      // make all but note section of jsd.
      switch (key->info->type)
        {
        case DRMS_TYPE_CHAR:
          if (colon) // probably type logical, leave as DRMS CHAR
  	    pjsd += sprintf(pjsd, "char, variable, record, DRMS_MISSING_VALUE, \"%%d\", \"none\", ");
          else
  	    pjsd += sprintf(pjsd, "int, variable, record, DRMS_MISSING_VALUE, \"%%d\", \"none\", ");
  	  break;
        case DRMS_TYPE_SHORT:
        case DRMS_TYPE_INT:
  	  pjsd += sprintf(pjsd, "int, variable, record, DRMS_MISSING_VALUE, \"%%d\", \"none\", ");
  	  break;
        case DRMS_TYPE_LONGLONG:
  	  pjsd += sprintf(pjsd, "longlong, variable, record, DRMS_MISSING_VALUE, \"%%lld\", \"none\", ");
  	  break;
        case DRMS_TYPE_FLOAT:
        case DRMS_TYPE_DOUBLE:
  	  pjsd += sprintf(pjsd, "double, variable, record, DRMS_MISSING_VALUE, \"%%f\", \"none\", ");
  	  break;
        case DRMS_TYPE_TIME:
  	  pjsd += sprintf(pjsd, "time, variable, record, DRMS_MISSING_VALUE, 0, \"UTC\", ");
  	  break;
        case DRMS_TYPE_STRING:
  	  pjsd += sprintf(pjsd, "string, variable, record, \"\", \"%%s\", \"none\", ");
  	  break;
        default:
  	  DIE("bad key type");
        }
      pjsd += sprintf(pjsd, "\"%s\"\n", key->info->description);
      }

    hiter_free(&hit);
  
    pjsd += sprintf(pjsd, "#======= Segments =======\n");
    pjsd += sprintf(pjsd, "Data: array, variable, %s, %d, ", drms_type2str(datatype), naxis);

    for (iaxis = 0; iaxis < naxis; iaxis++)
      pjsd += sprintf(pjsd, "%d, ", dims[iaxis]);
    pjsd += sprintf(pjsd, "\"\", fits, \"%s\", %f, %f, \"%s\"\n", 
                    (datatype != DRMS_TYPE_FLOAT && 
                     datatype != DRMS_TYPE_DOUBLE && 
                     datatype != DRMS_TYPE_LONGLONG) ? "compress Rice" : "",
                    bzero, bscale, in);
    pjsd += sprintf(pjsd, "#======= End JSD =======\n");

    if (printjsd)
      {
      printf("%s\n", jsd);
      // print keymap info
      printf("#====== BEGIN KEYNAME MAP =======\n");
      printf("# REMOVE these keyname map lines from the JSD\n");
      printf("# place the keyname map into a file for later use\n");
      printf("# mapfile has structure: wantedDRMSname namefromFITSfile action\n");
      printf("# but the second column needs to be the auto-converted name to provoke substitution\n");
      printf("# use \"copy\" for default action - without quotes\n");
      for (imap=0; imap<nmap; imap++)
          printf("# FITS name %s is converted to %s on input.\n%s\t%s\t%s\n",
            oldnames[imap], newnames[imap], newnames[imap], newnames[imap], actions[imap]);
      printf("#======END KEYNAME MAP =======\n");
      }
    }

  // if keyname mapfile is given, append to or replace names in map list
  if (havemap)
    {
    FILE *mapfile = fopen(map, "r");
    char line[10000];
    char newname[DRMS_MAXNAMELEN];
    char oldname[DRMS_MAXNAMELEN];
    char action[100];
    while (fgets(line, 1000, mapfile))
      {
      if (*line == '#')
        continue;
      if (sscanf(line,"%s%s%s", newname, oldname, action) != 3)
      {
         DIE("A mapfile line does not contain 3 words\n");
      }
      for (imap=0; imap<nmap; imap++)
        if (strcmp(oldname, oldnames[imap]) == 0)
          {
          if (newnames[imap]) free(newnames[imap]);
          newnames[imap] = strdup(newname);
          if (actions[imap]) free(actions[imap]);
          actions[imap] = strdup(action);
          break;
          }
      if (imap == nmap) // new name set found
        {
        if (nmap >= MAXMAPLEN)
        {
           DIE("Too many mapped keywords, increase MAXMAPLEN\n");
        }
        newnames[nmap] = strdup(newname);
        oldnames[nmap] = strdup(oldname);
        actions[nmap] = strdup(action);
        nmap++;
        }
      }
    fclose(mapfile);
    }

  if (wantcreate)
    {
    DRMS_Record_t *template;
    if (!haveseriesname || !haveprime)
    {
       DIE("Cant create series without ds and primekey args.\n");
    }
    if (drms_series_exists(drms_env , ds, &status))
    {
       DIE("Cant create existing series\n");
    }
    template = drms_parse_description(drms_env, jsd);
    if (template==NULL)
    {
       DIE("Failed to parse\n");
    }
    if (drms_create_series(template, 0))
    {
       DIE("Failed to create series.\n");
    }
    drms_free_record_struct(template);
    free(template);
    fprintf(stderr,"Series %s created.\n", ds);
    }

  if (insertrec)
    {
    char *usename;
    char *action;
    DRMS_RecordSet_t *rs;
    DRMS_Record_t *rec;
    if (!drms_series_exists(drms_env , ds, &status))
    {
       DIE("Series does not exist, cant insert record\n");
    }
    rs = drms_create_records(drms_env, 1, ds, DRMS_PERMANENT, &status);
    if (status)
    {
       DIE("Could not create new records in series");
    }
    rec = rs->records[0];
    hiter_new (&hit, keywords);
    while ( key = (DRMS_Keyword_t *)hiter_getnext(&hit) )
      {
      DRMS_Keyword_t *outkey;
      usename = key->info->name;
      action = "copy";
      for (imap=0; imap<nmap; imap++)
        if (strcmp(usename, oldnames[imap]) == 0)
          {
          usename = newnames[imap];
          action = actions[imap];
          break;
          }
      outkey = drms_keyword_lookup(rec, usename, 0);
      if (outkey)
        {
        if (strcmp(action, "copy") == 0)
          drms_setkey(rec, usename, key->info->type, &key->value);
//      else if ##### this is where you add new actions on keyword mapping
        else
          {
          fprintf(stderr, "old keyword %s has no action to make new key %s\n",
            key->info->name, usename);
          DIE("no action found for keyword\n");
          }
        }
      }

    hiter_free(&hit);
    status = drms_segment_write(drms_segment_lookupnum(rec,0), data, 0);

    if (status)
    {
      DIE("Could not write record");
    }
    drms_close_records(rs, DRMS_INSERT_RECORD);
    }

  drms_free_array(data);
  hcon_destroy(&keywords);
  return (DRMS_SUCCESS);
  }
