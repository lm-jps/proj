#include "jsoc_main.h"
#include "drms_types.h"

char *module_name = "export";

typedef enum
{
   kMymodErr_Success,
   kMymodErr_MissingArg,
   kMymodErr_DRMS
} MymodError_t;

#define kArg_recid     "reqid"
#define kArg_clname    "kmclass"
#define kArg_file      "kmfile"
#define kNotSpecified  "NOT SPECIFIED"
#define kArg_expSeries "expseries"
#define kDef_expSeries "jsoc.exports"

ModuleArgs_t module_args[] =
{
     {ARG_STRING, kArg_recid,  "", 
        "Export series primary key value that identifies the output record."},
     {ARG_STRING, kArg_expSeries, kDef_expSeries, "Series to which exported data are saved."},
     {ARG_STRING, kArg_clname,    kNotSpecified,  "Export key map class."},
     {ARG_STRING, kArg_file,      kNotSpecified,  "Export key map file."},
     {ARG_END}
};

ExportStrings_t gExpStr[] =
{
   {kExport_ReqID, "RequestID"},
   {kExport_Request, "DataSet"},
   {kExport_SegList, "Seg"},
   {kExport_Requestor, "Requestor"},
   {kExport_Notification, "Notify"},
   {kExport_ReqTime, "ReqTime"},
   {kExport_ExpTime, "ExpTime"},
   {kExport_DataSize, "Size"},
   {kExport_Format, "Format"},
   {kExport_FileNameFormat, "FilenameFmt"},
   {kExport_Status, "Status"},
   {(DRMS_ExportKeyword_t)-99, ""}
};

static int CallExportToFile(DRMS_Segment_t *segout, 
			    DRMS_Segment_t *segin, 
			    const char *clname,
			    const char *mapfile,
			    unsigned long long *szout)
{
   int status = DRMS_SUCCESS;
   char fileout[DRMS_MAXPATHLEN];
   char filein[DRMS_MAXPATHLEN];
   char *basename = NULL;
   unsigned long long size = 0;
   struct stat filestat;

   if (segout)
   {
      drms_segment_filename(segin, filein); /* input seg file name */
      if (!stat(filein, &filestat))
      {
	 size = filestat.st_size;
	 basename = rindex(filein, '/');
	 if (basename) 
	 {
	    basename++;
	 }
	 else 
	 {
	    basename = filein;
	 }

	 CHECKSNPRINTF(snprintf(segout->filename, DRMS_MAXSEGFILENAME, "%s", basename), DRMS_MAXSEGFILENAME);
	 drms_segment_filename(segout, fileout);

	 status = drms_segment_mapexport_tofile(segin, clname, mapfile, fileout);
      }
      else
      {
	 fprintf(stderr, "Unable to open source file '%s'.\n", filein);
	 status = DRMS_ERROR_EXPORT;
      }
   }

   *szout = 0;
   if (status == DRMS_SUCCESS)
   {
      *szout = size;
   }

   return status;
}

/* recout is the export series' record to which the exported data will be saved. */
static int MapexportRecord(DRMS_Record_t *recout,
                           DRMS_Record_t *recin,
                           const char *classname, 
                           const char *mapfile,
                           int *status)
{
   int stat = DRMS_SUCCESS;
   HIterator_t *hit = NULL;
   DRMS_Segment_t *segout = NULL;
   DRMS_Segment_t *segin = NULL;
   unsigned long long size = 0;
   unsigned long long tsize = 0;
   char dir[DRMS_MAXPATHLEN];

   segout = drms_segment_lookupnum(recout, 0);

   drms_record_directory(recin, dir, 1); /* This fetches the input data from SUMS. */

   if (segout)
   {
      /* The input rs query can specify a subset of all the series' segemnts - 
       * this is encapsulated in recin. */
      hit = hiter_create(&(recin->segments));
      while ((segin = hiter_getnext(hit)) != NULL)
      {
	 size = 0;
	 stat = CallExportToFile(segout, segin, classname, mapfile, &size);
	 if (stat != DRMS_SUCCESS)
	 {
	    fprintf(stderr, "Failure exporting segment '%s'.\n", segin->info->name);
	    break;
	 }
	 else
	 {
	    tsize += size;
	 }
      }

      hiter_destroy(&hit);
   }
   else
   {
      fprintf(stderr, "Export series contains no segment!\n");
      stat = DRMS_ERROR_EXPORT;
   }

   if (status)
   {
      *status = stat;
   }

   return tsize;
}		

/* recout is the export series' record to which the exported data will be saved. 
 * reqid is the primary key in the export series 
 */
static int Mapexport(DRMS_Env_t *env,
                     const char *reqid,
                     const char *classname, 
                     const char *mapfile,
                     const char *expseries,
                     int *status)
{
   int stat;
   int nSets = 0;
   int iSet = 0;
   int nRecs = 0;
   int iRec = 0;
   unsigned long long tsize = 0;
   DRMS_Record_t *recout = NULL;
   DRMS_Record_t *recin = NULL;
   DRMS_RecordSet_t *rs = NULL;
   DRMS_RecordSet_t *rsout = NULL;
   DRMS_RecordSet_t *rsin = NULL;
   char rsoutquery[DRMS_MAXQUERYLEN];
   char *rsinquery = NULL;

   snprintf(rsoutquery, sizeof(rsoutquery), "%s[%s]", expseries, reqid);

   rs = drms_open_records(env, rsoutquery, &stat);

   if (rs)
   {
      /* XXX Change to DRMS_REPLACE_SEGMENTS (which isn't implemented yet). */
      rsout = drms_clone_records(rs, DRMS_PERMANENT, DRMS_COPY_SEGMENTS, &stat);
      drms_close_records(rs, DRMS_FREE_RECORD);
   }
   
   if (rsout && rsout->n == 1)
   {
      recout = rsout->records[0];
      rsinquery = drms_getkey_string(recout, gExpStr[kExport_Request].str, &stat);
      
      if (rsinquery && (rsin = drms_open_records(env, rsinquery, &stat)))
      {
	 nSets = rsin->ss_n;

	 for (iSet = 0; stat == DRMS_SUCCESS && iSet < nSets; iSet++)
	 {
	    /* Perhaps we will need this in the future?
	     * request = rsin->ss_queries[iSet];
	     */
	    nRecs = drms_recordset_getssnrecs(rsin, iSet, &stat);

	    for (iRec = 0; stat == DRMS_SUCCESS && iRec < nRecs; iRec++)
	    {
	       recin = rsin->records[(rsin->ss_starts)[iSet] + iRec];
	       tsize += MapexportRecord(recout, 
                                        recin,
                                        classname, 
                                        mapfile, 
                                        &stat);
	    }
	 }
	 
	 if (stat != DRMS_SUCCESS)
	 {
	    fprintf(stderr, "Export halted due to DRMS failure.\n");
	 }
      }
      else
      {
	 fprintf(stderr, 
		 "Export series keyword '%s' did not contain a valid recordset query.\n", 
		 gExpStr[kExport_Request].str);
	 stat = DRMS_ERROR_INVALIDDATA; 
      }
   }
   else
   {
      fprintf(stderr, "Could not open export destination record set with query '%s'.\n", reqid);
      stat = DRMS_ERROR_INVALIDDATA;
   }

   /* Set output record keywords. */
   if (stat == DRMS_SUCCESS)
   {
      drms_setkey_time(recout, gExpStr[kExport_ExpTime].str, CURRENT_SYSTEM_TIME);
      drms_setkey_int(recout, gExpStr[kExport_DataSize].str, tsize);
      drms_setkey_int(recout, gExpStr[kExport_Status].str, 1);
   }

   if (rsout)
   {
      drms_close_records(rsout, DRMS_INSERT_RECORD);
   }

   if (rsinquery)
   {
      free(rsinquery);
   }

   if (status)
   {
      *status = stat;
   }
   
   return tsize;
}

int DoIt(void) 
{
   MymodError_t err = kMymodErr_Success;
   int drmsstat = DRMS_SUCCESS;
   long long tsize = 0;

   const char *expseries = NULL;
   const char *clname = NULL;
   const char *mapfile = NULL;
   const char *reqid = 0;

   reqid = cmdparams_get_str(&cmdparams, kArg_recid, &drmsstat);
   if (drmsstat != DRMS_SUCCESS)
   {
      fprintf(stderr, "Missing required argument '%s'.\n", kArg_recid);  
   }
   else
   {
      expseries = cmdparams_get_str(&cmdparams, kArg_expSeries, &drmsstat);

      clname = cmdparams_get_str(&cmdparams, kArg_clname, &drmsstat);
      if (drmsstat != DRMS_SUCCESS || !strcmp(clname, kNotSpecified))
      {
	 clname = NULL;
      }

      mapfile = cmdparams_get_str(&cmdparams, kArg_file, &drmsstat);
      if (drmsstat != DRMS_SUCCESS || !strcmp(mapfile, kNotSpecified))
      {
	 mapfile = NULL;
      }

      tsize = Mapexport(drms_env, reqid, clname, mapfile, expseries, &drmsstat);

      if (drmsstat != DRMS_SUCCESS)
      {
	 fprintf(stderr, "Failure occurred while processing export Request ID '%s'.\n", reqid);
      }
      else
      {
	 fprintf(stdout, "Successfully processed export Request ID '%s'.\n", reqid);
      }
   }

   if (drmsstat != DRMS_SUCCESS)
   {
      fprintf(stderr, "DRMS error '%d'.\n", drmsstat);
   }

   fprintf(stdout, "'%lld' bytes exported.\n", tsize);

   return err;
}
