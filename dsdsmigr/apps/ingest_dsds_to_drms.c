/* ingest_dsds_to_drms.c */

#include "jsoc_main.h"
#include "drms_types.h"
#include "dsdsmigr.h"
#include "printk.h"
#include <time.h>
#include <math.h>

char *module_name = "opendsrecs";

#define kRecSetIn	"in"
#define kRecSetOut	"out"
#define kNameList	"map"
#define kNOT_SPEC	"Not Specified"

#define DIE(msg)	return(fprintf(stderr,"%s",msg),1)
#define DIE_status(msg)	return(fprintf(stderr,"%s, status of death=%d",msg,status),1)

#define FDSCALE_DEFAULT (1.97784)
#define CANT_OPEN_FILE  (730)
#define BEYOND_TABLE_RANGE      (1)
#define EPHEMERIS_GAP   (2)
#define TRUNCATED_FIT   (3)

ModuleArgs_t module_args[] =
{
     {ARG_STRING, kRecSetIn, kNOT_SPEC, "Input data series."},
     {ARG_STRING, kRecSetOut, kNOT_SPEC, "Output data series."},
     {ARG_STRING, kNameList, kNOT_SPEC, "Name conversion list."},
     {ARG_STRING, "SCALE_CORRECTIONS", "", "scale correction recordset (1 record)"},
     {ARG_FLAG, "M", "0", "SkipMissingFiles - no records if DATAFILE is blank."},
     {ARG_FLAG, "v", "0", "verbose - more diagnostics"},
     {ARG_END}
};

double soho_obs_dist (TIME t) {
  return 1.0;
}

static int fgetline (char *line, FILE *tbl) {
/*
 *  A local utility for reading correction (or other) tables with the
 *    convention of one \n line delimited line per information record,
 *    comments with # as first character, and continuations flagged by
 *    \ at the end of the line
 */
  int c, b, cloc, len = 0;

  while ((c = fgetc (tbl)) != EOF) {
    if (c == '\n') {
      line[len++] = 0;
      if (line[0] == '#') line[len = 0] = 0;
      else if (b != '\\') return len;
      else {
	line[len=cloc] = ' ';
	line[++len] = 0;
      }
    } else {
      if (c != ' ' && c != '\t') {
	b = c;
	cloc = len;
      }
      line[len++] = c;
    }
  }
  return 0;
}

int MDI_plate_scale_adjustment (char *filepath, double *fdscl, TIME t, int *focus, int *qual,
    double obs_dist) {
  static struct plist {
    int qual;
    int focus;
    TIME t0;
    TIME t1;
    TIME tmid;
    double a0;
    double a1;
    double a2;
    double a3;
    int status;
  } *fitp;
  FILE *tbl;
  TIME t0, t1, t2, tmid;
  double a0, a1, a2, a3, rsun_adj;
  static int allocd = 32, file_read = 0;
  int status = 0;
  static int nct;
  int c, fields, m, n;
  char line[4096], tstr0[128], tstr1[128], tstr2[128];

  if (!file_read) {
     //    if (!(tbl = fopen ("/home/soi/CM/tables/calib/geom/scale_corrections", "r"))) {
     if (!(tbl = fopen (filepath, "r"))) {
      *fdscl = FDSCALE_DEFAULT;
      *qual = 0;
      *focus = 0;
      return CANT_OPEN_FILE;
    }
    fitp = (struct plist *)malloc (allocd * sizeof (struct plist));
    nct = 0;
/*
 * Reads a file in the format of multiple lines of the form
 *   Qual Focus T1  T2  T0  a0  a1  a2 a3
 * e.g.
 *   2 1 1995.12.21 1997.11.03 1997.02.09_14:31 1.005 -0.005  5.7e-12 -1.4e-19
 * and returns a set of n parameter sets
 */
    while (fgetline (line, tbl)) {
      fields = sscanf (line, "%d %d %s %s %s %lg %lg %lg %lg",
          qual, focus, tstr1, tstr2, tstr0, &a0, &a1, &a2, &a3);
      if (fields < 4) continue;
				/*  ignore lines without start & stop times  */
      fitp[nct].qual = *qual;
      fitp[nct].focus = *focus;
      fitp[nct].status = (fields < 9) ? TRUNCATED_FIT : 0;
      fitp[nct].t0 = sscan_time (tstr1);
      fitp[nct].t1 = sscan_time (tstr2);
      fitp[nct].tmid = (fields > 4) ? sscan_time (tstr0) : 0.0;
      fitp[nct].a0 = (fields > 5) ? a0 : 1.0;
      fitp[nct].a1 = (fields > 6) ? a1 : 0.0;
      fitp[nct].a2 = (fields > 7) ? a2 : 0.0;
      fitp[nct].a3 = (fields > 8) ? a3 : 0.0;
      if (++nct >= allocd) {
        allocd += 32;
	fitp = (struct plist *)realloc (fitp, allocd * sizeof (struct plist));
      }
    }
    fclose (tbl);
    file_read = 1;
  }

  if (obs_dist == 0.0) obs_dist = soho_obs_dist (t);
  if (obs_dist == 0.0) status = EPHEMERIS_GAP;
  n = 0;
  while (n < nct) {
    if (t >= fitp[n].t0 && t <= fitp[n].t1) break;
    n++;
  }
  if (n >= nct) {
    *fdscl = FDSCALE_DEFAULT;
    *qual = 0;
    *focus = 0;
    return BEYOND_TABLE_RANGE;
  }
  t -= fitp[n].tmid;
  rsun_adj = fitp[n].a0;
  rsun_adj += obs_dist * (fitp[n].a1 + t * (fitp[n].a2 + t * fitp[n].a3));
/*
  *fdscl *= FDSCALE_DEFAULT;
  	WRONG:  FDSCALE_DEFAULT is based on nominal rsun0 and corrected rsun
	is rsun0 * rsun_adj
*/
  *fdscl = FDSCALE_DEFAULT / rsun_adj;
  *qual = fitp[n].qual;
  *focus = fitp[n].focus;
  if (!status) {
    status = fitp[n].status;
  }
  return status;
}

// #include <sys/time.h>
TIME time_now()
  {
  TIME now;
  TIME UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
  struct timeval tp;
  gettimeofday(&tp, NULL);
  now = (double)tp.tv_sec + (double)tp.tv_usec/1.0e6;
  return(now +  UNIX_epoch);
  }

/* Name check code.  This code uses an external table to drive the mapping of
dsds keywords into drms keywords.  It should contain a row for each keyword in
the target drms series.  Each row should contain 3 fields, drms_name, dsds_name, action.
The action will be a string which matches the action table.  The action table
converts the string to an integer.
Sample file lines like:
  cadence CADENCE copy
  date DATE copy
Lines begging with "#" are comments and will be ignored.
Blank lines are not allowed.
*/

typedef struct NameListLookup_struct
  {
  char *drms_name;
  char *dsds_name;
  int action;
  struct NameListLookup_struct *next;
  } NameListLookup_t;

typedef enum {
  ACT_NOP, ACT_COPY, ACT_ANGLE, ACT_CENTER, ACT_TIME, ACT_AU, ACT_RESCALE, ACT_NOT_FOUND
  } Action_t;

Action_t actions[] = {
  ACT_NOP, ACT_COPY, ACT_ANGLE, ACT_CENTER, ACT_TIME, ACT_AU, ACT_RESCALE, ACT_NOT_FOUND
  };

char *action_names[] = {
  "nop",   "copy",   "pangle",   "center",   "time",  "au",   "scale", "done"
  };

int keyNameCheck(char *name, char **fromname)
  {
  static int first_call = 1;
  static NameListLookup_t *actionlist;
  NameListLookup_t *this;
  if (first_call)
    {
    char *actionlistname;
    FILE *flist;
    char drms_name[100], dsds_name[100], action_name[100];
    char line[1024];
    // int action;
    NameListLookup_t *last;

    first_call = 0;
    actionlistname = strdup(cmdparams_get_str(&cmdparams, kNameList, NULL));
    if (strcmp(actionlistname, kNOT_SPEC) == 0)
      {
      fprintf(stderr, "Name mapping list must be specified\n");
      exit(1);
      }
    flist = fopen(actionlistname, "r");
    if (!flist)
      {
      fprintf(stderr, "Name mapping file not found\n");
      exit(1);
      }
    last = actionlist = (NameListLookup_t *)malloc(sizeof(struct NameListLookup_struct));
    last->next = NULL;
    while (fgets(line, 1024, flist))
      {
      if (line[0] == '#')
        continue;
      else if (sscanf(line,"%s %s %s\n", drms_name, dsds_name, action_name)==3)
        {
        int iname;
	this = last;
        last = this->next = (NameListLookup_t *)malloc(sizeof(struct NameListLookup_struct));
        last->next = NULL;
        this->drms_name = strdup(drms_name);
        this->dsds_name = strdup(dsds_name);
        this->action = ACT_NOT_FOUND;
        for (iname=0; actions[iname] != ACT_NOT_FOUND; iname++)
          if (strcmp(action_name, action_names[iname]) == 0)
            {
            this->action = actions[iname];
            break;
            }
        }
      else fprintf(stderr, "Name map read error on line containing %s, ignoring this line\n",line);
      }
    }
  /* lookup name and return action */
  for (this=actionlist; this->next; this = this->next)
    {
    if (strcmp(this->drms_name, name) == 0)
      {
      *fromname = this->dsds_name;
      return(this->action);
      }
    }
  *fromname = name;
  return(ACT_NOT_FOUND);
  }

int DoIt(void) 
   {
   int status = DRMS_SUCCESS;
   int SkipMissingFiles;
   int verbose;
   int qualnodata=0x80000000;
   int nRecs, iRec;
   int qualstat;
   double val;
# define AU_m (149597870691.0)
   char *inRecQuery, *outRecQuery;
   DRMS_RecordSet_t *inRecSet, *outRecSet; 
   unsigned int quality;
   inRecQuery = strdup(cmdparams_get_str(&cmdparams, kRecSetIn, NULL));
   outRecQuery = strdup(cmdparams_get_str(&cmdparams, kRecSetOut, NULL));
   SkipMissingFiles = cmdparams_get_int(&cmdparams, "M", NULL) != 0;
   verbose = cmdparams_get_int(&cmdparams, "v", NULL) != 0;
   char scfilepath[DRMS_MAXPATHLEN];
   char infilepath[DRMS_MAXPATHLEN];

   if (strcmp(inRecQuery, kNOT_SPEC) == 0 || strcmp(outRecQuery, kNOT_SPEC) == 0)
      DIE("Both the "kRecSetIn" and "kRecSetOut" dataseries must be specified.\n");

   inRecSet = drms_open_records(drms_env, inRecQuery, &status);
   if (!inRecSet)
      DIE_status("Input dataseries not found\n");
   if ((nRecs = inRecSet->n) == 0)
      DIE("No input records found\n");
   printf("%d input records found\n", nRecs);
   char *screcquery = (char *)cmdparams_get_str(&cmdparams, "SCALE_CORRECTIONS", NULL);
   DRMS_RecordSet_t *screcset = drms_open_records(drms_env, screcquery, &status);
   if (status != DRMS_SUCCESS || screcset == NULL)
   {
     fprintf(stderr, "ERROR: problem reading scale corrections recordset: query = %s, status = %d\n", screcquery, status);
     return 1;
   }

   if (screcset->n != 1)
   {
     fprintf(stderr, "ERROR: scale corrections recordset contains more than one record: query = %s, nrecs = %d\n", screcquery, screcset->n);
     return 1;
   }

   DRMS_Segment_t *scseg = drms_segment_lookupnum(screcset->records[0], 0);
   drms_segment_filename(scseg,scfilepath);
   printf("scale corrections filepath = %s\n", scfilepath);


   DRMS_Segment_t *inseg = drms_segment_lookupnum(inRecSet->records[0], 0);
   drms_segment_filename(inseg,infilepath);
   printf("incoming filepath = %s\n", infilepath);


   for (iRec=0; iRec<nRecs; iRec++)
      {
      char *DataFile;
      char timebuf[1024];
      float UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
      float now;
      float result;
      //double result;
      int Record_OK = 1;
      DRMS_Record_t *inRec, *outRec;
      DRMS_Keyword_t *outKey;
      DRMS_Segment_t *inSeg, *outSeg;
      HIterator_t *outKey_last = NULL;
//      DRMS_Link_t *outLink;

      /* create output series rec prototype */
      inRec = inRecSet->records[iRec];
      outRecSet = drms_create_records(drms_env, 1, outRecQuery, DRMS_PERMANENT, &status);
      if (!outRecSet || outRecSet->n != 1)
         DIE_status("Output dataseries not found or can't create records\n");

      outRec = outRecSet->records[0];

      /* assume only one segment */
	DataFile = drms_getkey_string(inRec,"DATAFILE",&status);
	if (status && verbose)fprintf(stderr,"*** Segment Read DATAFILE status=%d\n",status);
	char filepath[DRMS_MAXPATHLEN];
	inSeg = drms_segment_lookupnum(inRec, 0);
        if (inSeg)
	  drms_segment_filename(inSeg, filepath);
        else 
          filepath[0] = '\0';
        //printf("filepath=%s\n",filepath);             
        //printf("ss=%d\n",access(filepath, R_OK | F_OK));
	if (*DataFile && access(filepath, R_OK | F_OK) == 0)
	  {
          outSeg = drms_segment_lookupnum(outRec, 0);
          if (inSeg && outSeg)
            {
        //    printf("ss=%d\n",iRec);
            DRMS_Array_t *data;
            /* read the data ad doubles so allow rescaling on output */
            data = drms_segment_read(inSeg, DRMS_TYPE_DOUBLE, &status);
            if (!data)
                  {
                     fprintf(stderr, "Bad data record %lld, status=%d\n",inRec->recnum, status);
                  DIE_status("giveup\n");
                  }
            /* use the zero and offset values in the JSD for the new record segment */
            data->bscale = outSeg->bscale;
            data->bzero = outSeg->bzero;
            drms_segment_write(outSeg, data, 0);
            drms_free_array(data);
            Record_OK = 1;    
	    quality = drms_getkey_int(inRec, "QUALITY", &qualstat);
            quality = quality & (~qualnodata);
            drms_setkey_int(outRec,"QUALITY",quality);
        //    printf("QUALITY=%08x\n",quality);
        //  drms_keyword_fprint(stderr, drms_keyword_lookup(iRec, "QUALITY",0));
        //    printf("QUALITY=%08x, DRMS_MISSING_INT=%08x\n",quality,DRMS_MISSING_INT);
	//    printf("qualstat=%d\n",qualstat);
        //    if (quality & 0x40000000) {
	//      fprintf (stderr, "Error: found bit 30 set in quality (%08x) for record %d\n",
	//	       quality, iRec);
	//      return 1;
	//       }
	//	if (quality & 0x8a000000) {
	//	  quality |= 0x40000000;
	//	  quality &= 0x7fffffff;
	//	}
	    }
          else
            DIE("Bad data segment lookup, in or out\n");
	  }
        else
	  { /* record is missing, copy t_rec and soho ephemeris keywords anyway*/
          drms_copykey(outRec, inRec, "T_REC");
	  drms_copykey(outRec, inRec, "OBS_VW");
	  drms_copykey(outRec, inRec, "OBS_VR");
	  drms_copykey(outRec, inRec, "OBS_VN");
	  sprint_time(timebuf, (double)time(NULL) + UNIX_epoch, "ISO", 0);
	  drms_setkey_string(outRec, "DATE", timebuf);
	  val = drms_getkey_double(inRec, "OBS_DIST", &status);
	  drms_setkey_double(outRec, "DSUN_OBS", val*AU_m);
	  val = drms_getkey_double(inRec, "OBS_B0", &status);
	  drms_setkey_double(outRec, "CRLT_OBS", val);	
	  val = drms_getkey_double(inRec, "OBS_CR", &status);
	  drms_setkey_double(outRec, "CAR_ROT", val);	
	  val = drms_getkey_double(inRec, "OBS_R0", &status);
	  drms_setkey_double(outRec, "RSUN_OBS", val);
  	  val = drms_getkey_double(inRec, "OBS_L0", &status);
	  drms_setkey_double(outRec, "CRLN_OBS", val); 
	  drms_setkey_int(outRec, "QUALITY", 0X80000000);
  	  val = drms_getkey_int(inRec, "DPC", &status);
	  drms_setkey_int(outRec, "DPC", val);
          /* commenting out the qualstat procedure
	  qualstat = 0;
          int quality = drms_getkey_int(outRec, "QUALITY", &qualstat);
	  if (!qualstat)
	    drms_setkey_int(outRec, "QUALITY", 0X80000000 | quality); 
          end commenting out the qualstat procedure*/
	  if (drms_keyword_lookup(outRec, "DATAVALS", 0))
	    drms_setkey_int(outRec, "DATAVALS", 0);
          if (SkipMissingFiles)
             {
             Record_OK = 0;
             if (verbose) 
               fprintf(stderr,"DSDS Record %d has no datafile, T_REC=%s, set missing.\n", iRec, drms_getkey_string(outRec,"T_REC",NULL));
             }
          else
             Record_OK = 1;
          drms_close_records(outRecSet,(Record_OK ? DRMS_INSERT_RECORD : DRMS_FREE_RECORD));
          continue;
	  }
      /* loop through all target keywords */
      outKey_last = NULL;
      while (outKey = drms_record_nextkey(outRec, &outKey_last, 1))
	{
	char *wantKey, *keyName = outKey->info->name;
        int action = keyNameCheck(keyName, &wantKey);
        if (!drms_keyword_inclass(outKey, kDRMS_KeyClass_Explicit))
	    continue;  // skip implicit keywords.
        switch (action)
          {
	  case ACT_NOP:
		break;
	  case ACT_COPY:
		{
		DRMS_Value_t inValue = {DRMS_TYPE_STRING, NULL};
		inValue = drms_getkey_p(inRec, wantKey, &status);
		if (status == DRMS_ERROR_UNKNOWNKEYWORD)
			break;
			if (status && verbose)fprintf(stderr,"*** ACT_COPY drms_getkey_p %s status=%d\n",wantKey,status);
		drms_setkey_p(outRec, keyName, &inValue);
		if ((inValue.type == DRMS_TYPE_STRING) && inValue.value.string_val)
		  free(inValue.value.string_val);
		inValue.value.string_val = NULL;
		break;
                }
	  case ACT_ANGLE:
		{ /* on CROTA2 set CROTA2, SAT_ROT, INST_ROT */
		double pangle, sat_rot;
		pangle = drms_getkey_double(inRec, "SOLAR_P", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_ANGLE drms_getkey_double SOLAR_P status=%d, pangle=%f\n",status,pangle);
		sat_rot = -pangle;
		drms_setkey_double(outRec, "CROTA2", sat_rot);
		drms_setkey_double(outRec, "SAT_ROT", sat_rot);
		drms_setkey_double(outRec, "INST_ROT", 0.0);
		break;
		}
	  case ACT_RESCALE:
	        { /* on X_SCALE, Y_SCALE, MAGNIFY, FD_SCALE set CDELT1, CDELT2 */
                double X_SCALE, Y_SCALE, fdscl, MAGNIFY, rsunobs;
		int qual, focus;
		double obs_dist;
		TIME t_obs;
                drms_setkey_string(outRec, "BLD_VERS", jsoc_version);
		obs_dist = drms_getkey_double(inRec, "OBS_DIST", &status);
                X_SCALE = drms_getkey_double(inRec, "X_SCALE", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_RESCALE drms_getkey_double X_SCALE status=%d, x0=%f\n",status,X_SCALE);
		Y_SCALE = drms_getkey_double(inRec, "Y_SCALE", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_RESCALE drms_getkey_double Y_SCALE status=%d, x0=%f\n",status,Y_SCALE);
		MAGNIFY = drms_getkey_double(inRec, "MAGNIFY", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_RESCALE drms_getkey_double MAGNIFY status=%d, x0=%f\n",status,MAGNIFY);
		/* it is critical that the status of the t_obs is what is fed into the two subsequent if statements */
                t_obs = drms_getkey_time(inRec, "T_OBS", &status);
                sprint_time(timebuf, t_obs, "ISO",0);
                //printf("t_obs=%s\n",timebuf);
                if (status && verbose)fprintf(stderr,"*** ACT_RESCALE drms_getkey_time t_obs status=%d, x0=%f\n",status,t_obs);
                //printf("qual=%d,focus=%d,obs_dist=%f,fdscl=%f,X_SCALE=%f,MAGNIFY=%f\n",qual,focus,obs_dist,fdscl,X_SCALE,MAGNIFY);
                if ( status == DRMS_SUCCESS)
		  {
                     status = MDI_plate_scale_adjustment(scfilepath, &fdscl, t_obs, &focus, &qual, obs_dist);
                        if (status != 0)
                          {
                             fprintf(stderr,"status of plate scale adjustment=%d\n",status);
                             fprintf(stderr,"incoming filepath with error = %s\n", infilepath);
                             return status;
                          }
		  }	
		drms_setkey_double(outRec, "CDELT1", X_SCALE*(fdscl/MAGNIFY));
		drms_setkey_double(outRec, "CDELT2", Y_SCALE*(fdscl/MAGNIFY));
                rsunobs=drms_getkey_double(inRec,"OBS_R0",&status);
                result=(rsunobs)/(X_SCALE*(fdscl/MAGNIFY));
                status=drms_setkey_float(outRec, "R_SUN", result);  
                //printf("result=%f,status=%d,qual=%d,focus=%d,obs_dist=%f,fdscl=%f,rsunobs=%f,X_SCALE=%f,MAGNIFY=%f\n",result,status,qual,focus,obs_dist,fdscl,rsunobs,X_SCALE,MAGNIFY);
	        }
	  case ACT_CENTER:
		{ /* on CRPIX1 set CRPIX1, CRPIX2, CRVAL1, CRVAL2 */
		double x0, y0;
		x0 = drms_getkey_double(inRec, "X0", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_CENTER drms_getkey_double X0 status=%d, x0=%f\n",status,x0);
		y0 = drms_getkey_double(inRec, "Y0", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_CENTER drms_getkey_double Y0 status=%d, y0=%f\n",status,y0);
		drms_setkey_double(outRec, "CRPIX1", x0+1.0);
		drms_setkey_double(outRec, "CRPIX2", y0+1.0);
		drms_setkey_double(outRec, "CRVAL1", 0.0);
		drms_setkey_double(outRec, "CRVAL2", 0.0);
		break;
		}
	  case ACT_TIME:
		{ /* on T_OBS set T_OBS, DATE-OBS, EXPTIME, CADENCE, TIME, MJD */
		char timebuf[1024];
		TIME MJD_epoch = -3727641600.000; /* 1858.11.17_00:00:00_UT  */
		TIME UNIX_epoch = -220924792.000; /* 1970.01.01_00:00:00_UTC */
                TIME t_obs, t_rec, date__obs, mjd, now; 
		double t_step;
		double exptime, mjd_day, mjd_time;
		t_rec = drms_getkey_time(inRec, "T_REC", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_TIME drms_getkey_time T_REC status=%d, t_rec=%f\n",status,t_rec);
		t_step = drms_getkey_double(outRec, "T_REC_step", &status); /* note from outRec */
		if (status && verbose)fprintf(stderr,"*** ACT_TIME drms_getkey_double T_REC_step status=%d, t_step=%f\n",status,t_step);
		// exptime = t_step; /* note - for lev1.5 */
                exptime = drms_getkey_double(inRec, "INTERVAL", &status);
		t_obs = drms_getkey_time(inRec, "T_OBS", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_TIME drms_getkey_time T_OBS status=%d, t_obs=%f\n",status,t_obs);
                if (status == DRMS_ERROR_UNKNOWNKEYWORD) // T_OBS is not present, skip this record.
		date__obs = t_obs - exptime/2.0;
                mjd = date__obs - MJD_epoch; /* sign error corrected by tplarson 2008.05.29 */
		mjd_day = floor(mjd / 86400.0);
		mjd_time = mjd - 86400.0 * mjd_day;
		now = (double)time(NULL) + UNIX_epoch;
		drms_setkey_time(outRec, "T_REC", t_rec);
		drms_setkey_time(outRec, "T_OBS", t_obs);
		drms_setkey_double(outRec, "EXPTIME", exptime);
		drms_setkey_double(outRec, "CADENCE", t_step);
		drms_setkey_double(outRec, "MJD", mjd_day);
		drms_setkey_double(outRec, "TIME", mjd_time);
                // allow either string or time types for DATE and DATE_OBS
                if (drms_keyword_type(drms_keyword_lookup(outRec, "DATE__OBS", 1)) == DRMS_TYPE_STRING)
                  {
		  sprint_time(timebuf, date__obs, "ISO", 0);
		  drms_setkey_string(outRec, "DATE__OBS", timebuf);
                  }
                else {
                  date__obs = t_obs - exptime/2.0;
		  drms_setkey_time(outRec, "DATE__OBS", date__obs);}
		if (drms_keyword_type(drms_keyword_lookup(outRec, "DATE", 1)) == DRMS_TYPE_STRING)
                  {
		  sprint_time(timebuf, now, "ISO", 0);
		  drms_setkey_string(outRec, "DATE", timebuf);
		  }
                else
                  drms_setkey_time(outRec, "DATE", now);

		break;
		}
	  case ACT_AU:
		{
# define AU_m (149597870691.0)
//#define AU_m (1.49597892e11) bad
		double au;
		au = drms_getkey_double(inRec, "OBS_DIST", &status);
		if (status && verbose)fprintf(stderr,"*** ACT_AU drms_getkey_double OBS_DIST status=%d,au=%f\n",status,au);
		if (status != DRMS_ERROR_UNKNOWNKEYWORD)
		    drms_setkey_double(outRec, "DSUN_OBS", au * AU_m);
		break;
		}
          case ACT_NOT_FOUND:
          default:
                /* name not in table, just take same name from input series */
                {
		DRMS_Value_t inValue = {DRMS_TYPE_STRING, NULL};
                DRMS_Keyword_t *outKey = drms_keyword_lookup(outRec, keyName, 0);
                if (drms_keyword_isconstant(outKey))
			break;
                inValue = drms_getkey_p(inRec, keyName, &status);
		if (status == DRMS_ERROR_UNKNOWNKEYWORD)
			break;
		if (status && verbose)fprintf(stderr,"*** DEFAULT drms_getkey_p %s status=%d\n",keyName, status);
                drms_setkey_p(outRec, keyName, &inValue);
		if ((inValue.type == DRMS_TYPE_STRING) && inValue.value.string_val)
		  free(inValue.value.string_val);
		inValue.value.string_val = NULL;
                break;
                }
          }
        }

/* populate the ROLL_TBL keyword */

char *rollrecstring;
char *rolltbls=drms_getkey_string(inRec, "ROLL_TBL", &status);

if ( status == DRMS_SUCCESS)
	{
        if (!strcmp("", rolltbls)) 
  		rollrecstring="";
 	if (!strcmp("not used", rolltbls)) 
  		rollrecstring="not used";
	if (!strcmp("/CM/tables/ephemeris/mdi_roll.smooth", rolltbls)) 
  		rollrecstring="mdi.roll_table";
	if (!strcmp("/CM/tables/ephemeris/mdi_roll.ascii", rolltbls)) 
  		rollrecstring="mdi.roll_table_ascii";
        drms_setkey_string(outRec, "ROLL_TBL", rollrecstring);
	}			



/* populate the CALTBLS keyword */

char *calrecstring;
char *caltbls=drms_getkey_string(inRec, "CALTBLS", &status);

if ( status == DRMS_SUCCESS)
{
if (!strcmp("", caltbls)) 
  calrecstring="";
if (!strcmp("/home/soi/CM/tables/calib/flat/fd/vers_0/", caltbls)) 
  calrecstring="mdi.caltables_intensity[fd_vers_0]";
if (!strcmp("/home/soi/CM/tables/calib/flat/fd/vers_1/", caltbls))
  calrecstring="mdi.caltables_intensity[fd_vers_1]";
if (!strcmp("/home/soi/CM/tables/calib/flat/hr/", caltbls))
  calrecstring="mdi.caltables_intensity[hr_vers_1]";
if (!strcmp("/home/soi/CM/tables/calib/flat/hr/vers_1", caltbls))
  calrecstring="mdi.caltables_intensity[hr_vers_1]";
if (!strcmp("/home/soi/CM/tables/calib/obflat/vers_0/", caltbls))
  calrecstring="mdi.caltables_intensity[obflat_vers_0]";
if (!strcmp("/home/soi/CM/tables/calib/obflat/vers_1/", caltbls))
  calrecstring="mdi.caltables_intensity[obflat_vers_1]";
if (!strcmp("/home/soi/CM/tables/calib/obflat/vers_2/", caltbls))
  calrecstring="mdi.caltables_intensity[obflat_vers_2]";
if (!strcmp("/home/soi/CM/tables/calib/dop/fd/tune_0/", caltbls))
  calrecstring="mdi.caltables_doppler[fd_tune_0]";
if (!strcmp("/home/soi/CM/tables/calib/dop/fd/tune_1/", caltbls))
  calrecstring="mdi.caltables_doppler[fd_tune_1]";
if (!strcmp("/home/soi/CM/tables/calib/dop/fd/tune_2/", caltbls))
  calrecstring="mdi.caltables_doppler[fd_tune_2]";
if (!strcmp("/home/soi/CM/tables/calib/dop/fd/tune_4/", caltbls))
  calrecstring="mdi.caltables_doppler[fd_tune_4]";
if (!strcmp("/home/soi/CM/tables/calib/dop/fd/tune_5/", caltbls))
  calrecstring="mdi.caltables_doppler[fd_tune_5]";
if (!strcmp("/home/soi/CM/tables/calib/dop/fd/tune_6/", caltbls))
  calrecstring="mdi.caltables_doppler[fd_tune_6]";
if (!strcmp("/home/soi/CM/tables/calib/dop/fd_spec/tune_0/", caltbls))
  calrecstring="mdi.caltables_doppler[fd_spec_tune_0]";
if (!strcmp("/home/soi/CM/tables/calib/dop/hr/tune_1/", caltbls))
  calrecstring="mdi.caltables_doppler[hr_tune_1]";
if (!strcmp("/home/soi/CM/tables/calib/dop/hr/tune_2/", caltbls))
  calrecstring="mdi.caltables_doppler[hr_tune_2]";
if (!strcmp("/home/soi/CM/tables/calib/dop/hr/tune_3/", caltbls))
  calrecstring="mdi.caltables_doppler[hr_tune_3]";
if (!strcmp("/home/soi/CM/tables/calib/dop/hr/tune_4/", caltbls))
  calrecstring="mdi.caltables_doppler[hr_tune_4]";
if (!strcmp("/home/soi/CM/tables/calib/dop/hr/tune_5/", caltbls))
  calrecstring="mdi.caltables_doppler[hr_tune_5]";
if (!strcmp("/home/soi/CM/tables/calib/dop/hr_spec/tune_0/", caltbls))
  calrecstring="mdi.caltables_doppler[hr_spec_tune_0]";
if (!strcmp("/home/soi/CM/tables/calib/dop/orig/tune_1/", caltbls))
  calrecstring="mdi.caltables_doppler_orig[orig_tune_1]";
if (!strcmp("/home/soi/CM/tables/calib/dop/orig/tune_2/", caltbls))
  calrecstring="mdi.caltables_doppler_orig[orig_tune_2]";
if (!strcmp("/home/soi/CM/tables/calib/dop/orig/tune_4/", caltbls))
  calrecstring="mdi.caltables_doppler_orig[orig_tune_4]";
if (!strcmp("/home/soi/CM/tables/calib/dop/orig/tune_5/", caltbls))
  calrecstring="mdi.caltables_doppler_orig[orig_tune_5]";
if (!strcmp("/home/soi/CM/tables/calib/dop/orig/tune_6/", caltbls))
  calrecstring="mdi.caltables_doppler_orig[orig_tune_6]";
drms_setkey_string(outRec, "CALTBLS", calrecstring); 
}
		
      /* loop through all target links */
      //while ( (outLink=(DRMS_Keyword_t *)hiter_getnext(&link_list)) )
        //{
        ///* assume no links - do nothing here now. */
        //}
      drms_close_records(outRecSet,(Record_OK ? DRMS_INSERT_RECORD : DRMS_FREE_RECORD));
      }
   drms_close_records(inRecSet,DRMS_FREE_RECORD);
   return(DRMS_SUCCESS);
   }
