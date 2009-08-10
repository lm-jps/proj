/* ingest_from_fits - read FITS file and print JSD, Keyword map file, or ingest file */
/*
 * ingest_from_fits [-j] [in=]<fitsfile> [map=<mapfile>] [ds=<series>]
 *   -j means print JSD for the fitsfile.
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
 */  

#include "jsoc_main.h"

char *module_name = "ingest_from_fits";

#define DIE(msg) {fflush(stdout);fprintf(stderr,"%s, status=%d\n",msg,status); return(status);}

ModuleArgs_t module_args[] =
{
     {ARG_STRING, "in", "NOT_SPECIFIED",  "Input FITS file."},
     {ARG_STRING, "ds", "NOT_SPECIFIED",  "Target DRMS data series."},
     {ARG_STRING, "map", "NOT_SPECIFIED",  "Map file for newname from oldname."},
     {ARG_FLAG, "j", "0",  "Print jsd and exit."},
     {ARG_END}
};

# define MAXJSDLEN 100000
# define MAXMAPLEN 10000

int DoIt(void)
  {
  char *in = params_get_str(&cmdparams, "in");
  char *ds = params_get_str(&cmdparams, "ds");
  char *map = params_get_str(&cmdparams, "map");
  int insertrec = strcmp(ds, "NOT_SPECIFIED") != 0;
  int wantjsd = params_isflagset(&cmdparams, "j") || !insertrec;
  int havemap = strcmp(map, "NOT_SPECIFIED") != 0;
  int status = DRMS_SUCCESS;
  DRMS_Array_t *data = NULL;
  DRMS_Keyword_t *key=NULL;
  HContainer_t *keywords = NULL;
  HIterator_t hit;
  int readraw = 1;
  char *newnames[MAXMAPLEN];
  char *oldnames[MAXMAPLEN];
  char *actions[MAXMAPLEN];
  int imap, nmap = 0;

  if (strcmp(in, "NOT_SPECIFIED") == 0)
    {
    if (cmdparams_numargs(&cmdparams) < 1 || !(in = cmdparams_getarg(&cmdparams, 1)))
      DIE("No input data found");
    }

  data = drms_fitsrw_read(in, readraw, &keywords, &status);
  if (status || !keywords)
    DIE("No keywords found");

  if (wantjsd)
    {
    char jsd[MAXJSDLEN], *pjsd = jsd;
    char keyname[DRMS_MAXNAMELEN];
    char *dash;
    int bitpix, iaxis, naxis, dims[10];
    double bzero, bscale;

    // build jsd in internal string
    pjsd += sprintf(pjsd, "#=====General Series Information=====\n");
    pjsd += sprintf(pjsd, "Seriesname:  <NAME HERE>\n");
    pjsd += sprintf(pjsd, "Author:      %s\n", getenv("USER"));
    pjsd += sprintf(pjsd, "Owner:       nobody_yet\n");
    pjsd += sprintf(pjsd, "Unitsize:    1\n");
    pjsd += sprintf(pjsd, "Archive:     0\n");
    pjsd += sprintf(pjsd, "Retention:   10\n");
    pjsd += sprintf(pjsd, "Tapegroup:   0\n");
    pjsd += sprintf(pjsd, "PrimeKeys:   <PRIME KEYS HERE OR DELETE LINE>\n");
    pjsd += sprintf(pjsd, "DBIndex:     <INDEX KEYS HERE OR DELETE LINE>\n");
    pjsd += sprintf(pjsd, "Description: \"From: %s\"\n", in);
  
    pjsd += sprintf(pjsd, "#===== Keywords\n");
  
    hiter_new (&hit, keywords);
    while ( key = (DRMS_Keyword_t *)hiter_getnext(&hit) )
      {
      strcpy(keyname, key->info->name);
      // Special actions for FITS basic structure
      if (strcmp(keyname, "SIMPLE") == 0) continue;
      if (strcmp(keyname, "EXTEND") == 0) continue;
      if (strcmp(keyname, "END") == 0) continue;
  
      // Special actions for external data representation
      if (strcmp(keyname, "BITPIX") == 0)
        {
        bitpix = drms2int(key->info->type, &key->value, NULL);
        continue;
        }
      if (strcmp(keyname, "NAXIS") == 0)
        {
        naxis = drms2int(key->info->type, &key->value, NULL);
        continue;
        }
      if (strncmp(keyname, "NAXIS", 5) == 0 && isdigit(keyname[5]))
        {
        dims[keyname[5] - '0' - 1] = drms2int(key->info->type, &key->value, NULL);
        continue;
        }
      if (strcmp(keyname, "BLANK") == 0) continue;
      if (strcmp(keyname, "BZERO") == 0)
        {
        bzero = drms2double(key->info->type, &key->value, NULL);
        continue;
        }
      if (strcmp(keyname, "BSCALE") == 0)
        {
        bscale = drms2double(key->info->type, &key->value, NULL);
        continue;
        }

      // check for illegal DRMS names 
      dash = index(keyname, '-');
      if (dash)
        {
        char *c = key->info->name;
        char *k = keyname;
        if (nmap >= MAXMAPLEN)
          DIE("Too many mapped keywords, increase MAXMAPLEN\n");
        while (*c)
          {
          if (*c == '-')
            {
            *k++ = '_';
            *k++ = '_';
            c++;
            }
          else
            *k++ = *c++;
          }
        while (k - keyname < 8)
          *k++ = '_';
        *k = '\0';
        newnames[nmap] = strdup(keyname);
        oldnames[nmap] = strdup(key->info->name);
        actions[nmap] = strdup("copy");
        nmap++; 
        }
      
      pjsd += sprintf(pjsd, "Keyword: %s, ", keyname);
      // make all but note section of jsd.
      switch (key->info->type)
        {
        case DRMS_TYPE_CHAR:
        case DRMS_TYPE_SHORT:
        case DRMS_TYPE_INT:
  	pjsd += sprintf(pjsd, "int, variable, record, DRMS_MISSING_VALUE, \"%%d\","
  	" \"none\", ");
  	break;
        case DRMS_TYPE_LONGLONG:
  	pjsd += sprintf(pjsd, "longlong, variable, record, DRMS_MISSING_VALUE, \"%%lld\","
  	" \"none\", ");
  	break;
        case DRMS_TYPE_FLOAT:
        case DRMS_TYPE_DOUBLE:
  	pjsd += sprintf(pjsd, "double, variable, record, DRMS_MISSING_VALUE, \"%%f\","
  	" \"none\", ");
  	break;
        case DRMS_TYPE_TIME:
  	pjsd += sprintf(pjsd, "time, variable, record, DRMS_MISSING_VALUE, 0,"
  	" \"UTC\", ");
  	break;
        case DRMS_TYPE_STRING:
  	pjsd += sprintf(pjsd, "string, variable, record, \"\", \"%%s\","
  	" \"none\", ");
  	break;
        default:
  	DIE("bad key type");
        }
      if (dash)
        pjsd += sprintf(pjsd, "\"[%s]\"\n", key->info->name);
      else
        pjsd += sprintf(pjsd, "\"\"\n");
      }
  
    pjsd += sprintf(pjsd, "#======= Segments =======\n");
    pjsd += sprintf(pjsd, "Data: array, variable, %s, %d, ",
      (bitpix == -64 ? "double" : 
      (bitpix == -32 ? "float"  :
      (bitpix == 8 ? "char" :
      (bitpix == 16 ? "short" :
      (bitpix == 32 ? "int" :
      (bitpix == 64 ? "longlong" : "UNKNOWN")))))), naxis);
    for (iaxis = 0; iaxis < naxis; iaxis++)
      pjsd += sprintf(pjsd, "%d, ", dims[iaxis]);
    pjsd += sprintf(pjsd, "\"\", fits, \"%s\", %f, %f, \"%s\"\n", 
      (bitpix > 0 ? "compress Rice" : ""), bzero, bscale, in);
    pjsd += sprintf(pjsd, "#======= End JSD =======\n");

    printf("%s\n", jsd);
    // print keymap info
    printf("#====== BEGIN KEYNAME MAP =======\n");
    printf("# REMOVE the keyname map lines from the JSD\n");
    printf("# place the keyname map into a file for later use\n");
    printf("# mapfile has structure: drmsname fitsname action\n");
    printf("# use \"copy\" for default action - without quotes\n");
    for (imap=0; imap<nmap; imap++)
        printf("%s\t%s\t%s\n", newnames[imap], oldnames[imap], actions[imap]);
    printf("#======END KEYNAME MAP =======\n");
    }

fprintf(stderr,"done with jsd part, do keymap part\n");
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
        DIE("A mapfile line does not contain 3 words\n");
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
          DIE("Too many mapped keywords, increase MAXMAPLEN\n");
        newnames[nmap] = strdup(newname);
        oldnames[nmap] = strdup(oldname);
        actions[nmap] = strdup(action);
        nmap++;
        }
      }
    fclose(mapfile);
    }

fprintf(stderr, "keymap info, nmap=%d\n",nmap);

  if (insertrec)
    {
    char *newname;
    char *action;
    DRMS_RecordSet_t *rs = drms_create_records(drms_env, 1, ds, DRMS_PERMANENT, &status);
    DRMS_Record_t *rec = rs->records[0];
    if (status)
      DIE("Could not create new records in series");
    hiter_new (&hit, keywords);
    while ( key = (DRMS_Keyword_t *)hiter_getnext(&hit) )
      {
      DRMS_Keyword_t *outkey;
      newname = key->info->name;
      action = "copy";
      for (imap=0; imap<nmap; imap++)
        if (strcmp(newname, oldnames[imap]) == 0)
          {
          newname = newnames[imap];
          action = actions[imap];
          }
      outkey = drms_keyword_lookup(rec, newname, 0);
      if (outkey)
        {
        if (strcmp(action, "copy") == 0)
          drms_setkey(rec, newname, key->info->type, &key->value);
//      else if ##### this is where you add new actions on keyword mapping
        else
          {
          fprintf(stderr, "old keyword %s has no action to make new key %s\n",
            key->info->name, newname);
          DIE("no action found for keyword\n");
          }
        }
      }
    status = drms_segment_write(drms_segment_lookupnum(rec,0), data, 0);
    if (status)
      DIE("Could not write record");
    drms_close_records(rs, DRMS_INSERT_RECORD);
    }
  return (DRMS_SUCCESS);
  }
