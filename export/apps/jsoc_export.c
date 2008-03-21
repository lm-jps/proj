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

ModuleArgs_t module_args[] =
{
     {ARG_STRING, kArg_recid,  "", "Export series primary key value that identifies the output record."},
     {ARG_STRING, kArg_clname, kNotSpecified, "Export key map class."},
     {ARG_STRING, kArg_file,   kNotSpecified, "Export key map file."},
     {ARG_END}
};

int DoIt(void) 
{
   MymodError_t err = kMymodErr_Success;
   int drmsstat = DRMS_SUCCESS;
   long long tsize = 0;

   const char *clname = NULL;
   const char *mapfile = NULL;
   const char *reqid = NULL;

   reqid = cmdparams_get_str(&cmdparams, kArg_recid, &drmsstat);
   if (drmsstat != DRMS_SUCCESS)
   {
      fprintf(stderr, "Missing required argument '%s'.\n", kArg_recid);  
   }
   else
   {
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

      tsize = drms_recordset_mapexport(drms_env, reqid, clname, mapfile, &drmsstat);

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
